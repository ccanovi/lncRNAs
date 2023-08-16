#!/bin/bash -l

## error and verbose
set -ex

## default args
mail="camilla.canovi@umu.se"
in=/mnt/picea/storage/reference/Pinus-tabuliformis/fasta/P.tabuliformis_V1.0.fa.gz
out=$(realpath ../indices)/Pinus-tabulliformis
proj=u2019016
singularity=$(realpath ../singularity/kogia/ncbi-blast_2.10.1+.sif)
## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

export SINGULARITY_BINDPATH=/mnt:/mnt

## prepare
sbatch --mail-user $mail -e $out/index.err -o $out/index.out -A $proj \
runBlastPlusMakeblastdb.sh -p nucl -t Pinus-tabulliformis $singularity $in $out
