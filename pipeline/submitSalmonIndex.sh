#!/bin/bash -l

# fail on ERROR
set -e

# load helpers
source $UPSCb/src/bash/functions.sh

# vars
email=camilla.canovi@umu.se
in=~/Git/lncRNAs/data/trinity/Trinity.fasta
out=~/Git/lncRNAs/data/SalmonIndex

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi


# exec
sbatch -e $outdir/SalmonIndex.err -o $outdir/SalmonIndex.out --mail-user $email $UPSCb-common/pipeline/runSalmonIndex.sh -t 10 $in $out
