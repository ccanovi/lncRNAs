#!/bin/bash -l

## be verbose and print
set -eux

## load the modules
module load bioinfo-tools diamond

## variables
proj=u2015037
mail=camilla.canovi@umu.se
query=~/Git/lncRNAs/data/Transdecoder/Trinity.fasta.transdecoder.pep
out=~/Git/lncRNAs/data/DIAMOND
ref=/mnt/picea/storage/reference/UniRef90/201908/indices/diamond/uniref90.dmnd
mapping=/mnt/picea/storage/reference/Taxonomy/20190825/prot.accession2taxid.gz
name=/mnt/picea/storage/reference/Taxonomy/20190825/names.dmp
node=/mnt/picea/storage/reference/Taxonomy/20190825/nodes.dmp

if [ ! -d $out ]; then
    mkdir -p $out
fi

## submit the job
#fnam=$(basename $3)_$(basename ${2//.f*a*/.blt})

sbatch --mail-user $mail -o $out/transdecoder.out -e $out/transdecoder.err \
-A $proj ../UPSCb-common/pipeline/runDiamond.sh -a -m $mapping -n $name -o $node blastp $query $ref $out
