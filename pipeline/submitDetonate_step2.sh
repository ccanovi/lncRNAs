#!/bin/bash -l

email=camilla.canovi@umu.se
right=../data/trinity/insilico_read_normalization_altogether/right.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_pctSD10000.fq.gz
left=../data/trinity/insilico_read_normalization_altogether/left.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_pctSD10000.fq.gz
fasta=../data/trinity/Trinity.fasta
proj=u2015037
out=../data/detonate
length=238
param=../data/detonate/detonate.params
## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

# submit the job2
sbatch -e $out/detonate_step2.err -o $out/detonate_step2.out --mem=16G \
--mail-user=$email -A $proj ../src/bash/runDetonate_2step.sh $left $right $fasta $out $length $param
