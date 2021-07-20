thing <- genes
library(Biostrings)
seq <- readDNAStringSet("/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all-phase.gff3.CDSandLTR-TE.fa")
names(seq) <- sub(" .*","",names(seq))
names(seq) <- sub("\\.1","",names(seq))
IDs <- thing$gene
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/genes_backbone2.fasta")

#for 1st degree neighbours
#stra <- as_tibble(bla_FDG3)
#genes3 <- stra[grep("MA_*",stra$value),]
IDs <- bla_FDG3[grep("MA",bla_FDG3)]
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/gene_FDN3.fasta")


thing <- miRNAs
library(Biostrings)
seq <- readDNAStringSet("sRNA/ShortStack_genome/Pabies_SE_miRNA.mature.fa")
names(seq) <- sub(" .*","",names(seq))
IDs <- thing$gene
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/miRNAs_backbone2.fasta")
seq <- readDNAStringSet("data/trinity/Trinity.fasta")

seq <- readDNAStringSet("data/trinity/Trinity.fasta")
names(seq) <- sub(" .*","",names(seq))
IDs <- bla_FDG3[grep("TRINITY",bla_FDG3)]
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/myguy3.fasta")
IDs <- "TRINITY_DN58094_c0_g2_i1"


mRNA <- readDNAStringSet("~/Git/lncRNAs/data/analysis/DE/genes_backbone2.fasta")
colSums(alphabetFrequency(mRNA))
sel <- alphabetFrequency(mRNA)[,"N"] != 0
mRNA[sel] <- DNAStringSet(gsub("N","A",mRNA[sel]))
writeXStringSet(mRNA,"~/Git/lncRNAs/data/analysis/DE/genes_backbone2_N_changed_to_A.fasta")
colSums(alphabetFrequency(mRNA))


  
