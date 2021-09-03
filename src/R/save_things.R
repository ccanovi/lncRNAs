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
IDs <- bla_FDG12[grep("MA",bla_FDG12)]
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/lncTar_12/gene_FDN12.fasta")


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
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/candidates.fasta")
IDs <- "TRINITY_DN37788_c0_g2_i2"
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/lncTar_12/myguy12.fasta")
IDs <- c("TRINITY_DN17439_c0_g1_i4","TRINITY_DN22829_c0_g1_i1","TRINITY_DN22966_c0_g2_i1","TRINITY_DN4106_c0_g1_i19",
         "TRINITY_DN4327_c0_g1_i8","TRINITY_DN1474_c0_g1_i2","TRINITY_DN3337_c0_g1_i6","TRINITY_DN39226_c0_g1_i1",
         "TRINITY_DN4293_c0_g1_i1","TRINITY_DN7954_c0_g2_i1","TRINITY_DN6930_c1_g1_i1","TRINITY_DN5946_c0_g1_i2")

mRNA <- readDNAStringSet("~/Git/lncRNAs/data/analysis/DE/genes_backbone2.fasta")
colSums(alphabetFrequency(seq))
sel <- alphabetFrequency(seq)[,"N"] != 0
seq[sel] <- DNAStringSet(gsub("N","A",seq[sel]))
writeXStringSet(mRNA,"~/Git/lncRNAs/data/analysis/DE/genes_backbone2_N_changed_to_A.fasta")
colSums(alphabetFrequency(mRNA))



  
