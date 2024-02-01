#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -w picea
#SBATCH --mem=100G

# submit the job
seidr view -n ~/Git/lncRNAs/data/seidr/backbone/nodelist ~/Git/lncRNAs/data/seidr/backbone/backbone-2-percent.sf -b -o ~/Git/lncRNAs/data/seidr/backbone/nodelist-only.sf