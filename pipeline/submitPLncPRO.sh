#!/bin/bash -l

email=camilla.canovi@umu.se
in=~/Git/lncRNAs/data/trinity/Trinity.fasta
out=~/Git/lncRNAs/data/PLncPRO
model=/mnt/picea/storage/reference/PLncPRO/dicot.model
proj=u2015037

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/plncpro.err -o $out/plncpro.out --mem=16G \
--mail-user=$email -A $proj ../UPSCb-common/pipeline/runPLncPRO.sh $in $model $bla $out
