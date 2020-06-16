#!/bin/bash -l
#SBATCH -p core -n 1 
#SBATCH --mem=16GB 
#SBATCH -t 2-00:00:00 
#SBATCH --mail-type=END,FAIL

# stop on error
set -e

USAGETXT="
runtrmap.sh <ref_gff> <query_gff> <out dir>
"

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

export KOGIA_DOCKER_OPTS="-v /mnt/picea:/mnt/picea"

 # test if we get 3 params
if [ $# -ne 3 ]; then
    echo "This function needs 3 arguments"
    usage
fi
# test that $1 (fasta) exists
if [ ! -f $1 ]; then
  abort "The first argument needs to be the ref filepath"
fi

# test that $2 (query file) exists
if [ ! -f $2 ]; then
    abort "The second argument needs to be the query file"
fi

# test that $3 (out dir) exists
if [ ! -d $3 ]; then
    abort "The third argument needs to be the output directory"
fi

# submit the job
~/Git/kogia/scripts/trmap $1 $2 -o $3/trmap.fasta 

