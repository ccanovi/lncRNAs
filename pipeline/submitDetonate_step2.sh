#!/bin/bash -l

email=camilla.canovi@umu.se
right=~/Git/lncRNAs/data/trinity/insilico_read_normalization_altogether/right.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_pctSD10000.fq.gz
left=~/Git/lncRNAs/data/trinity/insilico_read_normalization_altogether/left.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_pctSD10000.fq.gz
fasta=~/Git/lncRNAs/data/trinity/Trinity.fasta
proj=u2015037

# submit the job
#sbatch -e $out/detonate_step2.err -o $out/detonate_step2.out --mem=16G \
#--mail-user=$email -A $proj $fasta
rsem-eval-calculate-score --bowtie2 -p 4 \
  --transcript-length-parameters detonate.param --paired-end \
  $right \
  $left \
  $fasta 139_ORP 500
