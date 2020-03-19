#!/bin/bash -l

email=camilla.canovi@umu.se
fasta=~/Git/lncRNAs/data/trinity/Trinity.fasta
proj=u2015037
out=../data/detonate

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/detonate_step1.err -o $out/detonate_step1.out --mem=16G \
--mail-user=$email -A $proj ../src/bash/runDetonate.sh $fasta $out
