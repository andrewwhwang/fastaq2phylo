#!/bin/bash
usage ()
{
  echo 'Usage : main.sh FASTQ/A/5 [-r #READS_PER_SAMPLE] [-s #BOOTSTRAP_SAMPLES] [-t THRESHOLD] [-db nt/viruses] [-f filterdb]'
  exit
}
if [ "$#" -gt 10 ] || [ "$#" -lt 1 ] ; then
  usage
fi
while [ "$1" != "" ]; do
    case $1 in
        -r )            shift
        READS=$1
        ;;
        -s )            shift
        SAMPLES=$1
        ;;
        -t )            shift
        THRES=$1
        ;;
        -db )           shift
        DB=$1
        ;;
        -f )           shift
        FILTER=$1
        ;;
        * )             QUERY=$1
    esac
    shift
done
if [ "$READS" = "" ] ; then
    READS="0"
fi
if [ "$SAMPLES" = "" ] ; then
    SAMPLES="1"
fi
if [ "$THRES" = "" ] ; then
    THRES="0"
fi
if [ "$QUERY" = "" ] ; then
    usage
fi
if [ "$FILTER" = "" ] ; then
    FILTER='human_genomic'
fi
if [ ! -e "$QUERY" ] ; then
    echo "can't find $QUERY"
    usage
fi
if [ "$DB" = "" ] ; then
    DB="nt"
#elif [ "$DB" != 'nt' ] && [ "$DB" != 'viruses' ]; then
#    echo "invalid database"
#    usage
fi

if [ ! -d output ]; then
    mkdir output
fi
if [ ! -d "output/newick" ]; then
    mkdir output/newick
fi
if [ ! -d "output/pngs" ]; then
    mkdir output/pngs
fi
export BLASTDB="C:/Users/Andrew.Hwang/Desktop/fastaq2phylo/db/"
para=4  # paralllelize x4
getLineage ()
{
#    awk '/^>/ {seqtotal+=seqlen;seqlen=0;seq+=1;next;}{seqlen=seqlen+length($0)}END{print seq" sequences, total length "seqtotal+seqlen": average length = "(seqtotal+seqlen)/seq}' output/$1.fasta
    # seqs=$(awk '/^>/ {seq+=1}END{print seq}' output/$1.fasta)
    if [ $DB = 'nt' ] ; then
        TASK="blastn"
        ID="staxid"
    elif [ $DB = 'mir' ] ; then
        TASK="blastn-short"
        ID="sgi"
        ALTNAME="sseqid"
    elif [ $DB = 'hiv' ] ; then
        TASK="blastn"
        ID="sgi"
        ALTNAME="sseqid"
    elif [ $DB = 'blood' ] ; then
        TASK="blastn"
        ID="sgi"
    elif [ $DB = 'viruses' ] ; then
        TASK="blastn"
        ID="sgi"
    fi
    echo "##############################BLASTING REMAINING SEQUENCES TO $DB##############################"
    echo "spliting fasta into $para files"
    python scripts/fastaSplit.py -file "output/filtered_fasta.$1.fasta" -num $para -total $totalSeqs -filenum $1
    echo "blasting fasta sequences"
    for i in $(eval echo {0..$(expr $para - 1)}) ; do
        (blastn -task $TASK -query output/$1.$i.fasta -max_hsps 1 -max_target_seqs 1 -out output/blastout.$i.txt -db db/$DB \
        -outfmt "6 qseqid $ID sstart send slen $ALTNAME" -num_threads $(nproc) >/dev/null 2>&1; echo "part $i done") &
    done
    wait
    cat output/blastout.*.txt > output/blastout.txt
    echo "getting lineage from hits"
    python scripts/lineage.py -file 'output/blastout.txt' -dbType $DB -filenum $1 #> output/lineage.$1.txt
}

filterOut ()
{
    echo "##############################FILTERING OUT $FILTER##############################"
    echo "spliting filter fasta into $para files"
    python scripts/fastaSplit.py -file "output/$1.fasta" -num $para -total $totalSeqs -filenum $1
    echo "blast filtering fasta sequences"
    for i in $(eval echo {0..$(expr $para - 1)}) ; do
        (blastn -task blastn -query output/$1.$i.fasta -out output/filtered.$1.$i.out -evalue 0.1 -max_hsps 1 -num_alignments 0 -num_descriptions 1 -db db/$FILTER \
        -word_size 11 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 -num_threads $(nproc) >/dev/null 2>&1; echo "part $i done") &
    done
    wait
    > output/nohits.$1.ids
    for i in $(eval echo {0..$(expr $para - 1)}) ; do
        grep -B10 "***** No hits" output/filtered.$1.$i.out | grep '^Query=' | sed 's/^Query= //' >> output/nohits.$1.ids
    done
    python scripts/filter.py -blastout "output/nohits.$1.ids" -fa "output/$1.fasta" -num $1
}

worker ()
{
    filterOut $1
    if [ -s output/filtered_fasta.$1.fasta ] ; then
        totalSeqs=$(awk '/^>/ {seq+=1}END{print seq}' output/filtered_fasta.$1.fasta)
    else
        echo "No sequences left after filtering"
        exit
    fi
    getLineage $1
}

if [ -f $QUERY ] ; then
    filename=$(basename "$QUERY")
    extension="${filename##*.}"
    name="${filename%.*}"
    echo "##############################STARTING $name.$extension##############################"
    if [ "$extension" = "fastq" ] || [ "$extension" = "fq" ] ; then
        echo 'converting fastq into fasta'
        cat $QUERY | paste - - - - | cut -f1-2 | sed 's/^@/>/g' | tr '\t' '\n' > output/result.fasta
    elif [ "$extension" = "fasta" ] || [ "$extension" = "fa" ] || [ "$extension" = "fas" ] ; then
        cp $QUERY output/result.fasta
    elif [ "$extension" = "fast5" ] ; then
        echo 'converting fast5 into fasta'
        poretools fasta $QUERY > output/result.fasta
    else
        echo "file format must be either fastq or fasta"
        usage
    fi
elif [ -d $QUERY ] ; then
    echo 'converting fast5 into fasta'
    poretools fasta $QUERY > output/result.fasta
fi
totalSeqs=$(awk '/^>/ {seq+=1}END{print seq}' output/result.fasta)
echo "total number of sequences: $totalSeqs"
if [ $READS -ne '0' ] ; then
    if [ "$READS" -lt "$totalSeqs" ]; then
        totalSeqs=$READS
    fi
    for j in $(eval echo {1..$SAMPLES}) ; do
        echo "------------------------------PROCESSING SAMPLE $j------------------------------"
        echo "selecting $READS sequences at random"
        python scripts/randomFasta.py -file output/result.fasta -num $READS -total $totalSeqs -sampleNum $j #> output/$j.fasta
        worker $j
    done
else
    mv output/result.fasta output/0.fasta
    worker 0
fi
cat output/lineage.*.txt > output/lineage.txt
cat output/filtered_fasta.*.fasta > output/filtered_fasta.fasta
echo 'creating taxonomy tree'
python scripts/makeTree.py -file output/lineage.txt -thres $THRES -samples $SAMPLES -param "$name.$DB.$FILTER.$READS.$THRES"
find output/ -maxdepth 1 ! -name 'readme.txt' -and ! -name 'lineage.txt' -and ! -name 'blastout.txt' -and ! -name 'filtered_fasta.fasta' -type f -exec rm {} +
echo 'done!'
