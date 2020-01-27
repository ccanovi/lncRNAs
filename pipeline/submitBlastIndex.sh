#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## default args
mail="camilla.canovi@umu.se"
in=../swissprot/2019_11/fasta/uniprot_sprot.fasta
out=../swissprot/2019_11/indices/blast+

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch --mail-user $mail -e $out/index.err -o $out/index.out \
../UPSCb-common/pipeline/runBlastPlusMakeblastdb.sh -p prot -t swissprot $in $out

