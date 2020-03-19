#!/bin/bash -l
#SBATCH -p core -n 1 
#SBATCH --mem=16GB 
#SBATCH -t 1-00:00:00 
#SBATCH --mail-type=END,FAIL

# stop on error
set -e

USAGETXT="
runDetonate.sh <fasta file> <out dir> 
"

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

 # test if we get 2 params
if [ $# -ne 2 ]; then
    echo "This function needs 2 arguments"
    usage
fi
# test that $1 (fasta) exists
if [ ! -f $1 ]; then
  abort "The first argument needs to be the fasta filepath"
fi

# test that $2 (out dir) exists
if [ ! -d $2 ]; then
    abort "The fourth argument (output dir) needs to be an existing directory"
fi

# submit the job
rsem-eval-estimate-transcript-length-distribution $1 $2/detonate.params

