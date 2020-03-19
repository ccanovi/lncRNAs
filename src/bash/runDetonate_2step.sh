#!/bin/bash -l
#SBATCH -p core -n 12 
#SBATCH --mem=16GB 
#SBATCH -t 1-00:00:00 
#SBATCH --mail-type=END,FAIL

# stop on error
set -ex

## load modules
module load bioinfo-tools bowtie2

USAGETXT="
runDetonate_2step.sh <fwd fastq> <rev fastq> <fasta file> <out dir> <average length of reads> <param>
"

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# test if we get 6 params
if [ $# -ne 6 ]; then
    echo "This function needs 6 arguments"
    usage
fi
# test that $1 (left) exists
if [ ! -f $1 ]; then
  abort "The second argument needs to be the left filepath"
fi
# test that $2 (right) exists
if [ ! -f $2 ]; then
  abort "The first argument needs to be the right filepath"
fi
# test that $3 (fasta) exists
if [ ! -f $3 ]; then
  abort "The third argument needs to be the fasta filepath"
fi
# test that $4 (out dir) exists
if [ ! -d $4 ]; then
    abort "The fourth argument (output dir) needs to be an existing directory"
fi
# test that $5 (length) exists
if [ $5 -lt 1 ]; then
  abort "The fifth argument needs to be the average length of reads"
fi
# test that $6 (param) exists
if [ ! -f $6 ]; then
  abort "The sixth argument needs to be the param filepath"
fi
# submit the job
rsem-eval-calculate-score --bowtie2 -p 12 \
  --transcript-length-parameters $6 --paired-end \
  $1 $2 $3 $4/final_detonate $5
