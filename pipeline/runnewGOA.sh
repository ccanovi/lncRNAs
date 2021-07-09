#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 24:00:00
#SBATCH --mem 32GB
#SBATCH --mail-type FAIL
#SBATCH --mail-user=camilla.canovi@umu.se

# stop on error, be verbose and expand the commands
set -e -x

## usage
USAGETXT=\
"
	Usage: runnewGOA.sh <input_Folder> <go_file> 
	
	Options:
	            --inputFolder results dir from matrixPreparation
	            -G            GO file
"

# Check
if [ $# -ne 2 ]; then
    echo "This function needs 2 arguments"
    usage
fi

if [ ! -d $1 ]; then
  abort "The first argument needs to be the results dir from matrixPreparation"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be the GO filepath"
fi

# run matrixPreparation

singularity exec --bind /mnt:/mnt \
/mnt/picea/projects/singularity/gni-predictors.sif GNI_predictors \
newGOA \
--inputFolder $1 -H -F 5 -P 5 -G $2