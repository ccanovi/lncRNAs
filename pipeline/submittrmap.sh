#!/bin/bash -l

email=camilla.canovi@umu.se
proj=u2015037
out=~/Git/lncRNAs/data/trmap/CuffMerge
ref=$(realpath ../reference/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3)
query=$(realpath ../data/Cuffmerge/transcripts.gff3)

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/trmap.err -o $out/trmap.out --mem=16G \
--mail-user=$email -A $proj ../src/bash/runtrmap.sh $ref $query $out
