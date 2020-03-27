#!/bin/bash -l

email=camilla.canovi@umu.se
right=../data/trinity/insilico_read_normalization_altogether/right.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_pctSD10000.fq.gz
left=../data/trinity/insilico_read_normalization_altogether/left.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_pctSD10000.fq.gz
assembly=../data/trinity/Trinity.fasta
proj=u2015037
out=../data/Transrate
threads=4

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job2
sbatch -e $out/transrate.err -o $out/transrate.out --mem=16G \
--mail-user=$email -A $proj ../src/bash/runTransrate.sh $assembly $left $right $threads
