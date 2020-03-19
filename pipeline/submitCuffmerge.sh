#!/bin/bash -l

## global vars
mail="camilla.canovi@umu.se"
fasta=/mnt/picea/storage/reference/Picea-abies/v1.1/fasta/Pabies01-genome-collapsed-for-STAR.fa
gtf=/mnt/picea/storage/reference/Picea-abies/v1.1/gtf/Pabies1.1-gene.gtf
in=$(realpath ../data/Cufflinks)
out=$(realpath ../data/Cuffmerge)

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## run
#for f in $(find $in -name "transcripts.gtf");
#do
#    fnam=$(basename ${f})
    sbatch -e $out/cuffmerge.err -o $out/cuffmerge.txt \
    --mail-user $mail ../UPSCb-common/pipeline/runCuffmerge.sh \
    $in $out $fasta $gtf 
#done

