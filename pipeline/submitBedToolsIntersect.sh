#!/bin/bash -l

## be verbose and print
set -eux

## source functions
source $UPSCb/src/bash/functions.sh

## variables

mail=camilla.canovi@umu.se
in=~/Git/lncRNAs/data/GMAP/GMAP_all.gff3
out=~/Git/lncRNAs/data/GMAP/BedToolsIntersect2
ref=../reference/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene-gene-only.gff3
proj=u2015037

## load the modules
module load bioinfo-tools
module load BEDTools

TRINITY_DN761_c0_g1_i10 TRINITY_DN603_c0_g2_i2 TRINITY_DN532_c0_g1_i2 TRINITY_DN357_c0_g1_i6 TRINITY_DN28_c0_g1_i3
TRINITY_DN178927_c0_g1_i1 TRINITY_DN178912_c0_g1_i1 TRINITY_DN178936_c0_g1_i1 TRINITY_DN56_c1_g1_i5 TRINITY_DN53_c1_g1_i2

if [ ! -d $out ]; then
    mkdir -p $out
fi


## submit the job
sbatch --mail-user=$mail -e $out/BedToolsIntersect.err -o $out/BedToolsIntersect.out -A $proj ../UPSCb-common/pipeline/runBedToolsIntersect.sh $in $ref $out -wo
