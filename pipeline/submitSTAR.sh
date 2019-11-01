#!/bin/bash -l

## global vars
mail="camilla.canovi@umu.se"
genome=/mnt/picea/storage/reference/Picea-abies/v1.1/indices/STAR/v2.5.2b/Pabies01-genome
fasta=/mnt/picea/storage/reference/Picea-abies/v1.1/fasta/Pabies01-genome-collapsed-for-STAR.fa
# with TEs to check
gff3=/mnt/picea/storage/reference/Picea-abies/v1.1/gff3/Pabies1.1-gene.gff3
# without TE
# gff3=/mnt/picea/storage/reference/Picea-abies/v1.1/gff3/Pabies1.1-gene-wo-intron.gff3

## create the out dir
in=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/trimmomatic
out=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/lncRNAs/STAR

if [ ! -d $out ]; then
    mkdir -p $out
fi

## run
for f in `find $in -name "*_trimmomatic_1.fq.gz"`;
do
    fnam=$(basename ${f/_1.fq.gz/})
    sbatch --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out --mem=180GB -n 30 \
    -J STAR-$fnam ../UPSCb-common/pipeline/runSTAR.sh -l 100000000000 \
    -f gff3 -g $gff3 -m 69000 $out $genome $fasta $f $in/${fnam}_2.fq.gz 
done

