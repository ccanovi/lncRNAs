#!/bin/bash -l
#SBATCH -p core -n 1 
#SBATCH --mem=16GB 
#SBATCH -t 1-00:00:00 
#SBATCH --mail-type=END,FAIL

# stop on error
set -ex

USAGETXT="
runTransrate.sh <assembly> <left file> <right file> <threads>
"

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# test if we get 4 params
if [ $# -ne 4 ]; then
    echo "This function needs 4 arguments"
    usage
fi
# test that $1 (assembly) exists
if [ ! -f $1 ]; then
  abort "The first argument needs to be the assembly filepath"
fi
# test that $2 (left) exists
if [ ! -f $2 ]; then
  abort "The second argument needs to be the left filepath"
fi
# test that $3 (right) exists
if [ ! -f $3 ]; then
  abort "The third argument needs to be the right filepath"
fi
# test that $4 (threads) exists
if [ $4 -lt 1 ]; then
  abort "The fourth argument needs to be the threads"
fi
# submit the job
transrate --assembly $1 --left $2 --right $3 --threads $4

