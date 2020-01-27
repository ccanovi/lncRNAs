#!/bin/bash -l

## global vars
mail="camilla.canovi@umu.se"
fasta=/mnt/picea/storage/reference/Picea-abies/v1.1/fasta/Pabies01-genome-collapsed-for-STAR.fa
gtf=/mnt/picea/storage/reference/Picea-abies/v1.1/gtf/Pabies1.1-gene.gtf


## create the out dir
in=~/Git/lncRNAs/data/STAR
out=~/Git/lncRNAs/data/Cufflinks

if [ ! -d $out ]; then
    mkdir -p $out
fi

## run
for f in `find $in -name "*_STAR.bam"`;
do
    fnam=$(basename ${f/.bam/})
    sbatch --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out \
    ../UPSCb-common/pipeline/runCufflinks.sh -n \
    -i 70000 $out $fasta $f $gtf 
done

