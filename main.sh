#!/bin/bash
usage ()
{
  echo 'Usage : main.sh FASTQ/A [-r #READS_PER_SAMPLE] [-s #BOOTSTRAP_SAMPLES] [-t THRESHOLD] [-db nt/viruses]'
  exit
}
if [ "$#" -gt 9 ] || [ "$#" -lt 1 ] ; then
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
if [ ! -e "$QUERY" ] ; then
    echo "can't find $QUERY"
    usage
fi
if [ "$DB" = "" ] ; then
    DB="nt"
elif [ "$DB" != 'nt' ] && [ "$DB" != 'viruses' ]; then
    echo "invalid database"
    usage
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

getLineage ()
{
#    awk '/^>/ {seqtotal+=seqlen;seqlen=0;seq+=1;next;}{seqlen=seqlen+length($0)}END{print seq" sequences, total length "seqtotal+seqlen": average length = "(seqtotal+seqlen)/seq}' output/$1.fasta
    # seqs=$(awk '/^>/ {seq+=1}END{print seq}' output/$1.fasta)
    if [ $DB = 'viruses' ] ; then
        format="sgi"
    elif [ $DB = 'nt' ] ; then
        format="staxids"
    fi
    para=4  # paralllelize x4
    echo "spliting fasta into $para files"
    python scripts/fastaSplit.py -file "output/$1.fasta" -num $para -total $totalSeqs -filenum $1
    echo "blasting fasta sequences"
    for i in $(eval echo {0..$(expr $para - 1)}) ; do
        (blastn -query output/$1.$i.fasta -max_hsps 1 -max_target_seqs 1 -out output/blastout.$i.txt -db db/$DB -outfmt "10 qseqid $format sstart send slen" -num_threads $(nproc) >/dev/null 2>&1; echo "part $i done") & 
    done
    wait
    cat output/blastout.*.txt > output/blastout.txt
    echo "getting lineage from hits"
    python scripts/lineage.py -file 'output/blastout.txt' -dbType $DB -filenum $1 #> output/lineage.$1.txt
}

filename=$(basename "$QUERY")
extension="${filename##*.}"
name="${filename%.*}"
echo "##############################STARTING $name.$extension##############################"
if [ "$extension" = "fastq" ] || [ "$extension" = "fq" ] ; then
    echo 'converting fastq into fasta'
    cat $QUERY | paste - - - - | cut -f1-2 | sed 's/^@/>/g' | tr '\t' '\n' > output/result.fasta
elif [ "$extension" = "fasta" ] || [ "$extension" = "fa" ] || [ "$extension" = "fas" ] ; then 
    cp $QUERY output/result.fasta
else
    echo "file format must be either fastq or fasta"
    usage
fi

totalSeqs=$(awk '/^>/ {seq+=1}END{print seq}' output/result.fasta)
echo "total number of sequences: $totalSeqs"
if [ $READS -ne '0' ] ; then
    if [ "$READS" -lt "$totalSeqs" ]; then
        totalSeqs=$READS
    fi
    for j in $(eval echo {1..$SAMPLES}) ; do
        echo "------------------------------PROCESSING SAMPLE $j------------------------------"
        echo "selecting $READ sequences at random"
        python scripts/randomFasta.py -file output/result.fasta -num $READS -total $totalSeqs -sampleNum $j #> output/$j.fasta
        getLineage $j 
    done
else
    mv output/result.fasta output/0.fasta
    getLineage 0
fi
cat output/lineage.*.txt > output/lineage.txt
echo 'creating taxonomy tree'
python scripts/makeTree.py -file output/lineage.txt -thres $THRES -samples $SAMPLES -param "$name.$DB.$READS.$THRES"
find output/ -maxdepth 1 ! -name 'readme.txt' -and ! -name 'lineage.txt' -and ! -name 'blastout.txt' -type f -exec rm {} +
echo 'done!'

