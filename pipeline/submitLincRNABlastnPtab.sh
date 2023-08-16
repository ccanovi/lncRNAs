#!/bin/bash -l

## error and verbose
set -ex

## default args
mail="camilla.canovi@umu.se"
in=$(realpath ../data/analysis/DE/first.fasta)
out=$(realpath ../data/blastn_plaza)
inx=$(realpath ../indices/Pinus-tabulliformis/P.tabuliformis_V1.0.fa.gz)
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
-e $out/Pinus_tabulliformis_first.err -o $out/Pinus_tabulliformis_first.out -A $proj \
$(realpath ../UPSCb-common/pipeline/runBlastPlus.sh) -p $cpu $singularity blastn $in $inx $out 

