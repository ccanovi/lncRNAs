#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## default args
mail="nicolas.delhomme@umu.se"
in=$(realpath ../data/analysis/DE/linc_network.fasta)
out=$(realpath ../data/blastn)
inx=$(realpath ../PlantGenIE/genomic_blast_indices/aliases)/plaza
cpu=12

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch --mail-user $mail -p core -c $cpu \
-e $out/TE.err -o $out/TE.out \
$(realpath ../UPSCb-common/pipeline/runBlastPlus.sh) -p $cpu blastn $in $inx $out 

