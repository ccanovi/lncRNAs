#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## default args
mail="nicolas.delhomme@umu.se"
in=$(realpath ../reference/fasta/Repeats/LTR-TE_CopiaGypsy_whole-length.fa.bz2)
out=$(realpath ../data/blastn)
inx=$(realpath ../data/indices/blast+/linc_network.fasta)
cpu=8

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

bunzip2 -c $in > $out/TE.fa

## prepare
sbatch --mail-user $mail -p core -c $cpu \
-e $out/TE.err -o $out/TE.out \
$(realpath ../UPSCb-common/pipeline/runBlastPlus.sh) -p $cpu blastn $out/TE.fa $inx $out 

