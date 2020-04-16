#!/bin/bash -l
#SBATCH -p core -n 10
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -w picea

# submit the job
singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/delhomme-upscb-lncrna.simg python \
/opt/plncpro/bin/rf/predict.py /mnt/picea/home/ccanovi/Git/lncRNAs/data/trinity/Trinity.fasta_all_features \
/mnt/picea/storage/reference/PLncPRO /dicot.model /mnt/picea/home/ccanovi/Git/lncRNAs/data/PLncPRO/Trinity.txt false

e