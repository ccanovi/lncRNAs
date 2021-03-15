#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## default args
mail="nicolas.delhomme@umu.se"
in=$(realpath ../sRNA/ShortStack_genome/Pabies_SE_miRNA.precursor.fa)
out=$(realpath ../data/blastn)
inx=$(realpath ../data/indices/blast+/linc_network.fasta)
cpu=8

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch --mail-user $mail -p core -c $cpu \
-e $out/miRNA.err -o $out/miRNA.out \
$(realpath ../UPSCb-common/pipeline/runBlastPlus.sh) -p $cpu blastn $in $inx $out 

