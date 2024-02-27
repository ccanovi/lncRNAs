thing <- genes
library(Biostrings)
seq <- readDNAStringSet("/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all-phase.gff3.CDSandLTR-TE.fa")
names(seq) <- sub(" .*","",names(seq))
names(seq) <- sub("\\.1","",names(seq))
IDs <- "MA_125713g0010"
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/target_gene.fasta")

thing <- miRNAs
library(Biostrings)
seq <- readDNAStringSet("sRNA/ShortStack_genome/Pabies_SE_miRNA.mature.fa")
names(seq) <- sub(" .*","",names(seq))
IDs <- mi
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/miRNAs_network.fasta")
seq <- readDNAStringSet("data/trinity/Trinity.fasta")

load(here("data/analysis/seidr/network.rda"))
seq <- readDNAStringSet("data/trinity/Trinity.fasta")
names(seq) <- sub(" .*","",names(seq))
IDs <- rownames(dat)[grep("TRINITY",rownames(dat))]
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/linc_network.fasta")

seq <- readDNAStringSet("data/trinity/Trinity.fasta.gz")
names(seq) <- sub(" .*","",names(seq))
IDs <- names(seq)[names(seq) !="TRINITY_DN25410_c1_g1_i3"]
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/transcripts_all.fasta")

mRNA <- readDNAStringSet("~/Git/lncRNAs/data/analysis/DE/genes_backbone2.fasta")
colSums(alphabetFrequency(seq))
sel <- alphabetFrequency(seq)[,"N"] != 0
seq[sel] <- DNAStringSet(gsub("N","A",seq[sel]))
writeXStringSet(mRNA,"~/Git/lncRNAs/data/analysis/DE/genes_backbone2_N_changed_to_A.fasta")
colSums(alphabetFrequency(mRNA))

