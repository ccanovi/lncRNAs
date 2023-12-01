#!/bin/bash -l

## error and verbose
set -ex

# modules
#module load bioinfo-tools blast

## default args
mail="camilla.canovi@umu.se"
in=~/Git/lncRNAs/precursors/Pabies_SE_miRNA.precursor.fa
out=$(realpath ../precursors)
inx=$(realpath ../indices/new_blast+/linc_network.fasta)
singularity=$(realpath ../singularity/kogia/ncbi-blast_2.11.0+.sif)
cpu=8
proj=u2019016

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

export SINGULARITY_BINDPATH=/mnt:/mnt

## prepare
sbatch --mail-user $mail -p core -c $cpu \
-e $out/miRNA.err -o $out/miRNA.out -A $proj \
$(realpath ../UPSCb-common/pipeline/runBlastPlus.sh) -p $cpu $singularity blastn $in $inx $out 

