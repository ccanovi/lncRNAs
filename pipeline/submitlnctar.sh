#!/bin/bash -l

email=camilla.canovi@umu.se
#lncRNA=~/Git/lncRNAs/data/analysis/DE/lincRNAs_backbone2.fasta
in=~/Git/lncRNAs/data/lnctar
out=$in
mRNA=~/Git/lncRNAs/data/analysis/DE/genes_backbone2_N_changed_to_A.fasta
proj=u2015037

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/lnctar.err -o $out/lnctar.out -n 1 --mem=1G -a 1-1000 -t 48:00:00 \
-p rbx --mail-user=$email -A $proj runLnctarAsArray.sh $in/linc $mRNA $out
