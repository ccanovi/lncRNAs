#!/bin/bash -l
#SBATCH -p core
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00

# stop on error, be verbose and expand the commands
set -eux

## usage
USAGETXT=\
"
	Usage: runLncTar.sh <lncRNAs.fasta> <target_genes.fasta> <out dir>
	
	Options:
	            -l            lncRNAs fasta file
	            -m            mRNAs fasta file
	            -o            out file
"

# Check
if [ $# -ne 3 ]; then
    echo "This function needs 3 arguments"
    usage
fi

# array ID
aID=$(printf '%04d' $SLURM_ARRAY_TASK_ID)

if [ ! -f $1.$aID ]; then
  abort "The first argument needs to be the lncRNAs fasta filepath"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be the mRNA fasta filepath"
fi

if [ ! -d $3 ]; then
    abort "The third argument (output dir) needs to be an existing directory"
fi

out=$3/linc$aID
[[ ! -d $out ]] && mkdir -p $out

# run lnctar
cd $out

# output filename
singularity exec --bind /mnt:/mnt \
/mnt/picea/projects/singularity/lnctar.sif \
LncTar.pl \
-p 1 -l $1.$aID -m $2 -d -0.09 -s T -o $out/all.txt
