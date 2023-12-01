thing <- genes
library(Biostrings)
seq <- readDNAStringSet("/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all-phase.gff3.CDSandLTR-TE.fa")
names(seq) <- sub(" .*","",names(seq))
names(seq) <- sub("\\.1","",names(seq))
IDs <- "MA_125713g0010"
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/target_gene.fasta")

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
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/candidates_PR.fasta")
IDs <- "TRINITY_DN37788_c0_g2_i2"
writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/lncTar_12/myguy12.fasta")
IDs <- c("TRINITY_DN23326_c0_g1_i1","TRINITY_DN16429_c0_g1_i2","TRINITY_DN58094_c0_g2_i1","TRINITY_DN9869_c0_g1_i5",
         "TRINITY_DN19049_c0_g1_i1","TRINITY_DN52747_c1_g1_i1","TRINITY_DN24450_c0_g1_i1","TRINITY_DN58806_c0_g1_i1",
         "TRINITY_DN135093_c0_g1_i1","TRINITY_DN18510_c0_g1_i13")

mRNA <- readDNAStringSet("~/Git/lncRNAs/data/analysis/DE/genes_backbone2.fasta")
colSums(alphabetFrequency(seq))
sel <- alphabetFrequency(seq)[,"N"] != 0
seq[sel] <- DNAStringSet(gsub("N","A",seq[sel]))
writeXStringSet(mRNA,"~/Git/lncRNAs/data/analysis/DE/genes_backbone2_N_changed_to_A.fasta")
colSums(alphabetFrequency(mRNA))



#getting the TPM values for PlantGenie

#sapply(list.files(here("data/salmon"),recursive = TRUE,
 #                 pattern = "quant.sf",full.names = TRUE))

files <- list.files(here("data/salmon"), 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)

# filter the files
files <- files[match(samples$ScilifeID,sub("_sortm.*","",basename(dirname(files))))]


library(parallel)
something <-do.call(cbind,mclapply(files,function(f){
  
  read_tsv(f,show_col_types = FALSE) %>%
#    filter(grepl("^MA*",Name)) %>% 
    select(TPM) 
                         
                         #write_tsv(here(tpm,"data/analysis/salmon/tpm_values.tsv"))
                         
                         },mc.cores=16L))

samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)
something <- do.call(
  cbind,
  lapply(split.data.frame(t(something),
                          samples$ID),
         colSums))
colnames(something) <- c("Stage1","Stage1","Stage1","Stage2","Stage2","Stage2","Stage3","Stage3","Stage3","Stage4","Stage4","Stage4",
                         "Stage5","Stage5","Stage5","Stage6","Stage6","Stage6","Stage7","Stage7","Stage7","Stage8","Stage8","Stage8")
rownames(something) <- read_tsv(files[[length(files)]],show_col_types = FALSE) %>%
  select(Name) %>% unlist()

rownames(something) <- gsub("\\.\\d+$", "", rownames(something))
write.csv(something,file=here("data/analysis/salmon/tpm_values.csv"),row.names = TRUE)
laST_BJ <- read.csv(here("data/analysis/salmon/tpm_values.csv"),col.names = 
                      c("genes","Stage1","Stage1","Stage1","Stage2","Stage2","Stage2","Stage3","Stage3","Stage3","Stage4","Stage4","Stage4",
                        "Stage5","Stage5","Stage5","Stage6","Stage6","Stage6","Stage7","Stage7","Stage7","Stage8","Stage8","Stage8"))
