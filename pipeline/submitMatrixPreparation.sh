#!/bin/bash -l

email=camilla.canovi@umu.se
out=~/Git/lncRNAs/functional_prediction/results
proj=u2015037
go=~/Git/lncRNAs/functional_prediction/go.obo
annotation=~/Git/lncRNAs/functional_prediction/annotation.txt
network=~/Git/lncRNAs/functional_prediction/edgelist-weighted.txt

module load Armadillo/9.600.6
module load boost/1.72

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

sbatch -e ~/Git/lncRNAs/functional_prediction/matrixPreparation.err -o ~/Git/lncRNAs/functional_prediction/matrixPreparation.out --mem=32G \
--mail-user=$email -A $proj \
runMatrixPreparation.sh $go $annotation $network $out