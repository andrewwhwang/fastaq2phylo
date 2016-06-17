#!/bin/bash
BASEDIR=$(dirname "$0")
echo $BASEDIR
usage ()
{
  echo 'Usage : filter.sh FASTQ/A [-r #READS_PER_SAMPLE] [-s #BOOTSTRAP_SAMPLES] [-t THRESHOLD] [-db DATABASE]'
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
#elif [ "$DB" != 'nt' ] && [ "$DB" != 'viruses' ]; then
#    echo "invalid database"
#    usage
fi

if [ ! -d $BASEDIR/output ]; then
    mkdir $BASEDIR/output
fi
if [ ! -d "$BASEDIR/output/newick" ]; then
    mkdir $BASEDIR/output/newick
fi
if [ ! -d "$BASEDIR/output/pngs" ]; then
    mkdir $BASEDIR/output/pngs
fi
getLineage ()
{
#    awk '/^>/ {seqtotal+=seqlen;seqlen=0;seq+=1;next;}{seqlen=seqlen+length($0)}END{print seq" sequences, total length "seqtotal+seqlen": average length = "(seqtotal+seqlen)/seq}' output/$1.fasta
    # seqs=$(awk '/^>/ {seq+=1}END{print seq}' output/$1.fasta)
    >$BASEDIR/output/blastout.txt
    if [ $DB = 'viruses' ] ; then
        format="sgi"
    elif [ $DB = 'nt' ] ; then
        format="staxids"
    fi
    para=8  # paralllelize x8
    echo "$(tput setaf 1)spliting fasta$(tput sgr0)"
    ./$BASEDIR/scripts/split.sh "$BASEDIR/output/$1.fasta" $BASEDIR
    echo "$(tput setaf 1)blasting fasta sequences$(tput sgr0)"
    parallel --eta -P $para "blastn -query {} -max_hsps 1 -max_target_seqs 1 -out {.}.tmp -db $BASEDIR/db/$DB -outfmt \"10 qseqid $format sstart send slen score\" -num_threads $(nproc); cat {.}.tmp >> $BASEDIR/output/blastout.txt && rm {.}.tmp" ::: $BASEDIR/output/*.fa
    echo "$(tput setaf 1)filtering out human reads$(tput sgr0)"
    python $BASEDIR/scripts/filter.py -blastout "$BASEDIR/output/blastout.txt" -fa $BASEDIR/output/$1.fasta
    numLinesBefore=$(wc -l < $BASEDIR/output/$1.fasta)
    numLinesAfter=$(wc -l < $BASEDIR/output/filtered_fasta.fasta)
    lineDiff=$((numLinesBefore - numLinesAfter))
    echo "got rid of $((lineDiff / 2)) sequences"
}

filename=$(basename "$QUERY")
extension="${filename##*.}"
name="${filename%.*}"
echo "$(tput setaf 2)##############################STARTING $name.$extension##############################$(tput sgr0)"
if [ "$extension" = "fastq" ] || [ "$extension" = "fq" ] ; then
    echo 'converting fastq into fasta'
    cat $QUERY | paste - - - - | cut -f1-2 | sed 's/^@/>/g' | tr '\t' '\n' > $BASEDIR/output/tmp.fasta
elif [ "$extension" = "fasta" ] || [ "$extension" = "fa" ] || [ "$extension" = "fas" ] ; then 
    echo 'reformating fasta'
    cp  $QUERY $BASEDIR/output/tmp.fasta
else
    echo "file format must be either fastq or fasta"
    usage
fi

python $BASEDIR/scripts/reformat_fasta.py -file $BASEDIR/output/tmp.fasta > $BASEDIR/output/result.fasta
rm $BASEDIR/output/tmp.fasta

totalSeqs=$(awk '/^>/ {seq+=1}END{print seq}' $BASEDIR/output/result.fasta)
echo "total number of sequences: $totalSeqs"
if [ $READS -ne '0' ] ; then
    if [ "$READS" -gt "$totalSeqs" ]; then
        READS=$totalSeqs
    fi
    for j in $(eval echo {1..$SAMPLES}) ; do
        echo "------------------------------PROCESSING SAMPLE $j------------------------------"
        echo "selecting $READS sequences at random"
        python $BASEDIR/scripts/randomFasta.py -file $BASEDIR/output/result.fasta -num $READS -total $totalSeqs -sampleNum $j
        rm $BASEDIR/output/result.fasta
        getLineage $j 
    done
else
    mv $BASEDIR/output/result.fasta $BASEDIR/output/0.fasta
    getLineage 0
fi
find $BASEDIR/output/ -maxdepth 1 ! -name 'filtered_fasta.fasta' -type f -exec rm {} +
echo 'done!'

