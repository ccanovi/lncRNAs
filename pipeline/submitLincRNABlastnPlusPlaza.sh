#!/bin/bash -l

## error and verbose
set -ex

## default args
mail="camilla.canovi@umu.se"
in=$(realpath ../data/analysis/DE/linc_network.fasta)
out=$(realpath ../data/blastn_plaza)
inx=/mnt/picea/storage/reference/PlantGenIE/genomic_blast_indices/aliases/plaza
singularity=$(realpath ../singularity/kogia/ncbi-blast_2.11.0+.sif)
cpu=12
proj=u2019016

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

export SINGULARITY_BINDPATH=/mnt:/mnt

## prepare
sbatch --mail-user $mail -p core -c $cpu \
-e $out/plaza.err -o $out/plaza.out -A $proj \
$(realpath runBlastPlus.sh) -p $cpu $singularity blastn $in $inx $out 

