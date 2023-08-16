#!/bin/bash -l

## be verbose and print
set -eux

## source functions
source ../UPSCb-common/src/bash/functions.sh

## variables

mail=camilla.canovi@umu.se
in=~/Git/lncRNAs/data/GMAP/genes_all_sorted.gff3
out=~/Git/lncRNAs/data/GMAP/BedToolsClosest
ref=~/Git/lncRNAs/data/GMAP/ref_genes_all.gff3
proj=u2015037
## load the modules
module load bioinfo-tools
module load BEDTools

if [ ! -d $out ]; then
    mkdir -p $out
fi


## submit the job
sbatch --mail-user=$mail -e $out/BedToolsClosest.err -o $out/BedToolsClosest.out -A $proj runBedToolsClosest.sh $in $ref $out -d 

