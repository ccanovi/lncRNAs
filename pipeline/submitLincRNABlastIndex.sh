#!/bin/bash -l

## error and verbose
set -ex

## default args
mail="nicolas.delhomme@umu.se"
in=$(realpath ../data/analysis/DE/linc_network.fasta)
out=$(realpath ../data)/indices/blast+

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# KOGIA is not exactly working that way, problem with linked dirs
export PATH=$PATH:$(realpath ../UPSCb-common/kogia/scripts)
alias makeblastdb='ncbi-blast makeblastdb'

## prepare
sbatch --mail-user $mail -e $out/index.err -o $out/index.out \
$(realpath ../UPSCb-common/pipeline/runBlastPlusMakeblastdb.sh) -p nucl -t lincRNA $in $out
