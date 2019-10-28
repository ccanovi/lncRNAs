#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 20
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00

## stop on error and undefined vars
## be verbose
set -eux

## load the modules
module load bioinfo-tools diamond

## options
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids"
PROC=20
TNODE=
TMAP=
TNAME=
EVALUE=0.001
MAXTGTSEQ=20
OPTIONS=""

## usage
usage(){
echo >&2 \
"
	Usage: $0 [options] <diamond command> <fasta file> <db> <out dir>

	Options:
		-a report all titles (sets --salltitles)
		-e the e-value (default to $EVALUE)
		-f the output format (default to $FMT)
		-k the max number of alignments (default to $MAXTGTSEQ)
		-m the path to the NCBI Taxonomy prot.accession2taxid.gz file
		-n the path to the NCBI Taxonomy names.dmp file
		-o the path to the NCBI Taxonomy nodes.dmp file
    -p number of threads to use (default $PROC)

    Note:
	     Providing the -f option is broken at the moment
"
	exit 1
}

## get the options
while getopts ae:f:k:m:n:o:p: option
do
    case "$option" in
      a) OPTIONS="$OPTIONS --salltitles";;
      e) EVALUE=$OPTARG;;
	    f) FMT=$OPTARG;;
	    k) $MAXTGTSEQ=$OPTARG;;
	    m) TMAP=$OPTARG;;
	    n) TNAME=$OPTARG;;
	    o) TNODE=$OPTARG;;
	    p) PROC=$OPTARG;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# != 4 ]; then
    echo "This function takes one diamond command, one file, one index and one directory as arguments"
    usage
fi

BLAST=
case "$1" in
    blastp)
	;;
    blastx)
	;;
    *)
	echo "Unknown diamond command. Not one of blastp, blastx. Aborting."
	usage
esac
BLAST=$1

if [ ! -f $2 ]; then
    echo "The second argument needs to be the fasta file"
    usage
fi

if [ ! -f $3 ] ; then
    echo "The third argument needs to be the path of the DIAMOND index"
    usage
fi

if [ ! -d $4 ]; then
    echo "The forth argument needs to be the output directory"
    usage
fi

if [ ! -z $TMAP ] && [ ! -f $TMAP ]; then
    echo "The -m argument expect an existing file"
    usage
else
  OPTIONS="$OPTIONS --taxonmap $TMAP"
fi

if [ ! -z $TNAME ] && [ ! -f $TNAME ]; then
    echo "The -n argument expect an existing file"
    usage
else
  OPTIONS="$OPTIONS --taxonnames $TNAME"
fi

if [ ! -z $TNODE ] && [ ! -f $TNODE ]; then
    echo "The -o argument expect an existing file"
    usage
else
  OPTIONS="$OPTIONS --taxonnodes $TNODE"
fi

## run BLAST
fnam=$(basename $3)_$(basename ${2//.f*a*/.blt})
diamond $1 --db $3 --query $2 --out $4/$fnam --evalue $EVALUE -k $MAXTGTSEQ \
--compress 1 --more-sensitive $OPTIONS --threads $PROC --outfmt $FMT

