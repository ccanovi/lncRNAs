#!/bin/bash -l

email=camilla.canovi@umu.se
out=~/Git/lncRNAs/data/Transdecoder
in=~/Git/lncRNAs/data/trinity/Trinity.fasta
proj=u2015037

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/transdecoder.err -o $out/transdecoder.out --mem=16G \
--mail-user=$email -A $proj ../UPSCb-common/pipeline/runTrinityTransDecoder.sh $in $out
