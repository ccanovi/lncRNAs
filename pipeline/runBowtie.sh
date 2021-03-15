#/bin/bash -l
# SBATCH -p core
# SBATCH -n 8
# SBATCH -o bowtie.out
# SBATCH -e bowtie.err
# SBATCH -J bowtie

set -eux
CPU=8

module load bioinfo-tools samtools bowtie

in=$(realpath ../sRNA/ShortStack_genome/Pabies_SE_miRNA.mature.fa)
inx=$(realpath ../data/indices/bowtie)/lincRNA
out=$(realpath ../data/bowtie)/Pabies_SE_miRNA.mature

[[ ! -d $(dirname $out) ]] && mkdir -p $(dirname $out)

sed "s:U:T:g" $in > $out.fa

bowtie -f -v 0 -p $CPU -S -a -m 100 $inx $out.fa | samtools view -bS - | samtools sort -o ${out}.noMM.bam -

bowtie -f -v 1 -p $CPU -S -a -m 100 $inx $out.fa | samtools view -bS - | samtools sort -o ${out}.oneMM.bam -

