#!/bin/bash -l

email=camilla.canovi@umu.se
lncRNA=~/Git/lncRNAs/data/analysis/DE/lincRNAs_backbone2.fasta
out=~/Git/lncRNAs/data/lnctar
mRNA=~/Git/lncRNAs/data/analysis/DE/genes_backbone2_N_changed_to_A.fasta
proj=u2015037

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/lnctar.err -o $out/lnctar.out --mem=16G \
--mail-user=$email -A $proj runlnctar.sh $lncRNA $mRNA $out
