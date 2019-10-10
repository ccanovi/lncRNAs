#!/bin/bash -l

email=camilla.canovi@umu.se

out=~/Git/lncRNAs/data/GMAP
inxDir=~/Git/lncRNAs/reference/indices/GMAP
inxName=Pabies1.0
in=~/Git/lncRNAs/data/trinity/Trinity.fasta
proj=u2015037

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -p core -n 12 -e $out/gmap.err -o $out/gmap.out --mem=100G \
--mail-user=$email -A $proj ../UPSCb-common/pipeline/runGmapl.sh -i 70000 -p 12 $in $inxDir $inxName $out
