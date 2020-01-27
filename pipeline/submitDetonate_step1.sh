#!/bin/bash -l

email=camilla.canovi@umu.se
fasta=~/Git/lncRNAs/data/trinity/Trinity.fasta
proj=u2015037
+
# submit the job
#sbatch -e $out/detonate_step1.err -o $out/detonate_step1.out --mem=16G \
#--mail-user=$email -A $proj $fasta
rsem-eval-estimate-transcript-length-distribution \
  $fasta \
  detonate.param
