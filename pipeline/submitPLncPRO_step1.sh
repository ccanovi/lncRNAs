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
sbatch singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/delhomme-upscb-lncrna.simg python /opt/plncpro/bin/blastparse_mt3.py /mnt/picea/home/ccanovi/Git/lncRNAs/data/trinity/Trinity.fasta_blastres \
-e $out/plncpro_step1.err -o $out/plncpro_step1.out --mem=16G \
--mail-user=$email -A $proj -t 2-00:00:00
