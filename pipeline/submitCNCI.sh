#!/bin/bash -l

email=camilla.canovi@umu.se

fasta=~/Git/lncRNAs/data/trinity/Trinity.fasta
out=~/Git/lncRNAs/data/CNCI

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -p core -n 12 -e $out/cnci.err -o $out/cnci.out --mem=16G \
--mail-user=$email $UPSCb-common/pipeline/runCNCI.sh $fasta $out
