#!/bin/bash -l
#SBATCH -p core -n 20
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -w picea

email=camilla.canovi@umu.se
out=~/Git/lncRNAs/data/PLncPRO
proj=u2015037

## load modules
module load bioinfo-tools diamond

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job
sbacth -e $out/plncpro_step2.err -o $out/plncpro_step2.out --mem=16G \
--mail-user=$email -A $proj -t 2-00:00:00 ../src/bash/runPLncPRO_step2.sh
