#!/bin/bash -l

email=camilla.canovi@umu.se

## be verbose and print
set -ex

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
}

## source functions
source $UPSCb/src/bash/functions.sh

## process the argument
in=~/Git/lncRNAs/data/trimmomatic
ref=~/Git/lncRNAs/data/SalmonIndex
out=~/Git/lncRNAs/data/Salmon
bind=/mnt:/mnt
img=/mnt/picea/projects/singularity/salmon.simg
proj=u2015037

## check vars
if [ -z $UPSCb ]; then
    abort "The UPSCb var needs to be set."
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## for every file
cd $out
for f in $(find $in -name "*sortmerna_trimmomatic_1.fq.gz"); do
  fnam=$(basename ${f/_1.fq.gz/})
  
  ## execute
  sbatch --mail-user=$email \
  -e $out/$fnam.err -o $out/$fnam.out \
  -A $proj ../UPSCb-common/pipeline/runSalmon.sh -b $bind \
  -i $img $ref $f $in/${fnam}_2.fq.gz $out

done
