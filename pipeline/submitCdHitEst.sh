#!/bin/bash -l

mail=nicolas.delhomme@umu.se
proj=u2015037
in=$(realpath ../data/analysis/DE/linc_network.fasta)
out=$(realpath ../data)/cdhit

[[ ! -d $out ]] && mkdir -p $out

export PATH=$PATH:$(realpath ../UPSCb-common/kogia/scripts)

sbatch -A $proj --mail-user=$mail --mail-type=END,FAIL \
-J cdhit -e $out/cdhit.err -o $out/cdhit.out \
$(realpath ../UPSCb-common/pipeline/runCdHitEst.sh) -c 0.9 $in $out/$(basename ${in/.fasta/_id90.fasta})