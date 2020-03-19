#!/bin/bash -l

email=camilla.canovi@umu.se
in=~/Git/lncRNAs/data/trinity/Trinity.fasta
out=~/Git/lncRNAs/data/PLncPRO
model=/mnt/picea/storage/reference/PLncPRO/dicot.model
db=/mnt/picea/storage/reference/UniRef90/201908/indices/diamond/uniref90.dmnd 
proj=u2015037

## load modules
module load bioinfo-tools diamond

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/plncpro.err -o $out/plncpro.out --mem=16G \
--mail-user=$email -A $proj ../UPSCb-common/pipeline/runPLncPRO.sh $in $model $db $out
