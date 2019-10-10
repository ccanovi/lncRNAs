#!/bin/bash -l

email=camilla.canovi@umu.se

fasta=~/Git/lncRNAs/data/trinity/Trinity.fasta
out=~/Git/lncRNAs/data/CPC2
outfile=$out/results.txt
proj=u2015037

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -p core -n 12 -w picea -e $out/cpc2.err -o $out/cpc2.out --mem=16G \
--mail-user=$email -A $proj ../UPSCb-common/pipeline/runCPC2.sh $fasta $outfile
