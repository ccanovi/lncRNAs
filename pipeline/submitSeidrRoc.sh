#!/bin/bash -l

account=u2019016
email=camilla.canovi@umu.se

# process the argument
aggregate=$(realpath ../data/seidr/aggregated/aggregated.sf)
backbone=$(realpath ../data/seidr/backbone)
out=$(realpath ../data/seidr/roc)
gs_pos=$(realpath ../data/seidr_new/Picea-abies_KEGG-based-positive-gold-standard.tsv)
gs_neg=$(realpath ../data/seidr_new/Picea-abies_KEGG-based-negative-gold-standard.tsv)

if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
for j in {1..10}; do
# for backbone files
  sbatch -A $account -o $out/${j}-percent.out \
  -e $out/${j}-percent.err -J roc-${j}  \
  ../UPSCb-common/pipeline/runSeidrRoc.sh $backbone/backbone-${j}-percent.sf \
  $gs_pos $gs_neg $out/${j}-percent.roc
  done

# for aggregated files
sbatch -A $account --mail-user $email -o $out/aggregated.out \
  -e $out/aggregated.err -J roc-aggr  --mem=100G \
  ../UPSCb-common/pipeline/runSeidrRoc.sh $aggregate $gs_pos \
 $gs_neg $out/aggregated.roc
