#!/bin/bash -l
#SBATCH -t 48:00:00

mail=camilla.canovi@umu.se
proj=u2019016
in=$(realpath ../data/analysis/DE/linc_network.fasta)
out=$(realpath ../data)/cdhit_new
singularity=$(realpath ../singularity/kogia/cd-hit-est_4.8.1.sif)

[[ ! -d $out ]] && mkdir -p $out

export SINGULARITY_BINDPATH=/mnt:/mnt

sbatch -A $proj --mail-user=$mail --mail-type=END,FAIL \
-J cdhit -e $out/cdhit90.err -o $out/cdhit90.out \
runCdHitEst.sh -c 0.9 $singularity $in $out/$(basename ${in/.fasta/_id90.fasta})