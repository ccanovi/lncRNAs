#/bin/bash -l
# SBATCH -p core
# SBATCH -n 8
# SBATCH -o bowtie-index.out
# SBATCH -e bowtie-index.err
# SBATCH -J bowtie-index

module load bioinfo-tools bowtie

in=$(realpath ../data/analysis/DE/linc_network.fasta)
out=$(realpath ../data/indices)/bowtie/lincRNA

[[ ! -d $out ]] && mkdir -p $out

bowtie-build $in $out