#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## check for the UPSCb env. var.
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## default args
mail="camilla.canovi@umu.se"
in=~/Git/lncRNAs/data/blastn/tmp/Trinity.fasta
out=~/Git/lncRNAs/data/blastn
inx=/mnt/picea/storage/reference/NCBI/20190507/nt
cpu=8

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch --mail-user $mail -p core --array=1-390 -c $cpu -e $out/Trinity.fa_%a.err -o $out/Trinity.fasta_%a.out \
$UPSCb-common/pipeline/runBlastPlus2.sh -p $cpu blastn $in $inx $out 

