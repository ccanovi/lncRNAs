#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem 40GB
#SBATCH --mail-type FAIL
#SBATCH --mail-user=camilla.canovi@umu.se

# stop on error, be verbose and expand the commands
set -e -x

## usage
USAGETXT=\
"
	Usage: runMatrixPreparation.sh <go_file> <annotation> <network_file>
	
	Options:
	            -G            GO file
	            -A            annotation file
	            -N            network file
	            -O           out file
"

# Check
if [ $# -ne 4 ]; then
    echo "This function needs 4 arguments"
    usage
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the GO filepath"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be the annotation filepath"
fi

if [ ! -f $3 ]; then
    abort "The third argument needs to be the network filepath"
fi

if [ ! -d $4 ]; then
    abort "The third argument (output dir) needs to be an existing directory"
fi

# run matrixPreparation

cd $4

singularity exec --bind /mnt:/mnt \
/mnt/picea/projects/singularity/gni-predictors.sif GNI_predictors \
matrixPreparation \
-G $1 -A $2 -N $3 -O $4
