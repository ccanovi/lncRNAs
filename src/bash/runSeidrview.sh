#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -w picea
#SBATCH --mem=100G
# submit the job
seidr view --nodelist ~/Git/lncRNAs/data/seidr/backbone/reduce_list_clean_noexcel.tsv --in-file ~/Git/lncRNAs/data/seidr/backbone/backbone-2-percent.sf -N -d $'\t' | cut -f 1,2,26 > ~/Git/lncRNAs/data/seidr/backbone/reducedIndexEdgeList.txt
