#!/bin/bash -l

email=camilla.canovi@umu.se
proj=u2015037
go=~/Git/lncRNAs/functional_prediction/go.obo
in=~/Git/lncRNAs/functional_prediction/results
# submit the job
sbatch -e ~/Git/lncRNAs/functional_prediction/runnewGOA.err -o ~/Git/lncRNAs/functional_prediction/runnewGOA.out --mem=32G \
--mail-user=$email -A $proj \
runnewGOA.sh $in $go
