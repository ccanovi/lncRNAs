#!/bin/bash -l

email=camilla.canovi@umu.se
fasta=~/Git/lncRNAs/data/trinity/Trinity.fasta
out=~/Git/lncRNAs/data/PLEK
proj=u2015037

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -p core -n 12 -e $out/plek.err -o $out/plek.out --mem=16G \
--mail-user=$email -A $proj ../UPSCb-common/pipeline/runPLEK.sh $fasta $out
