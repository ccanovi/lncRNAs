#!/bin/bash -l
in=$(realpath ../data/seidr/results)
sf=$(realpath ../data/seidr/sf)
out=$(realpath ../data/seidr/aggregated)
mail=camilla.canovi@umu.se
account=u2019016

if [ ! -d $out ]; then
  mkdir -p $out
fi

if [ ! -d $sf ]; then
  mkdir -p $sf
  find $in -name "*.sf" -exec ln -s "{}" $sf \;
fi


sbatch --mem=128GB -o $out/aggregated.out -e $out/aggregated.err --mail-user=$mail -A $account \
$(realpath ../UPSCb-common/pipeline/runSeidrAggregate.sh) $out $sf/*.sf

