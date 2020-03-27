#!/bin/bash -l
email=camilla.canovi@umu.se
out=~/Git/lncRNAs/data/PLncPRO
proj=u2015037

# create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbatch -e $out/plncpro_step1.err -o $out/plncpro_step1.out --mem=16G \
--mail-user=$email -A $proj -t 2-00:00:00 runPLncPRO_step1.sh
