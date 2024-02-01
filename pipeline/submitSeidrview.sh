#!/bin/bash -l

email=camilla.canovi@umu.se
proj=u2015037

## load modules
module load bioinfo-tools clp seidr-devel


# submit the job
sbatch -e ~/Git/lncRNAs/data/seidr/backbone/seidrview_nodelist.err -o ~/Git/lncRNAs/data/seidr/backbone/seidrview_nodelist.out \
--mail-user=$email -A $proj -t 2-00:00:00 ../src/bash/runSeidrview_nodelist.sh
