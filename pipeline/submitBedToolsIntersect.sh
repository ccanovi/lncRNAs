#!/bin/bash -l

## be verbose and print
set -eux

## source functions
source $UPSCb/src/bash/functions.sh

## variables

mail=camilla.canovi@umu.se
in=~/Git/lncRNAs/data/GMAP/GMAP_all.gff3
out=~/Git/lncRNAs/data/GMAP/BedToolsIntersect
ref=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene-gene-only.gff3


## load the modules
module load bioinfo-tools
module load BEDTools


if [ ! -d $out ]; then
    mkdir -p $out
fi


## submit the job
sbatch --mail-user=$mail -e $out/BedToolsIntersect.err -o $out/BedToolsIntersect.out $UPSCb-common/pipeline/runBedToolsIntersect.sh $in $ref $out -s  -wo