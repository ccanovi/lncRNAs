#!/bin/bash -l

email=camilla.canovi@umu.se

outdir=~/Git/lncRNAs/data/FrameDP
infile=~/Git/lncRNAs/data/trinity/Trinity.fasta

# Nico to check - how to regenerate config and trained set
CFG=/mnt/picea/storage/reference/UniRef90/201406/framedp/cfg/FrameDP.cfg
ref_folder=/mnt/picea/storage/reference/UniRef90/201406/framedp/TRAIN_MAT 
options=--no_train

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/framedp.err -o $out/framedp.out --mem=16G \
--mail-user=$email $UPSCb-common/pipeline/runFrameDP.sh $infile $outdir $CFG $ref_folder $options
