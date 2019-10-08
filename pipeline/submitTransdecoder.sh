#!/bin/bash -l

email=camilla.canovi@umu.se

out=~/Git/lncRNAs/data/Transdecoder
in=~/Git/lncRNAs/data/trinity/Trinity.fasta

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/transdecoder.err -o $out/transdecoder.out --mem=16G \
--mail-user=$email $UPSCb/pipeline/runTrinityTransDecoder.sh $in $out
