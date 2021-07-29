#' ---
#' title: "non-coding RNAs metadata gathering"
#' author: "Camilla Canovi & Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Libraries
suppressPackageStartupMessages({
  library(Biostrings)
  library(here)
  library(magrittr)
  library(parallel)
  library(tidyverse)
  library(tximport)
})
  
#' Helper functions
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))

#' # Data
#' ## Salmon
salmon <- list.files(here("data/Salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

counts <- suppressMessages(round(tximport(files = salmon, 
                                          type = "salmon",
                                          txOut=TRUE)$counts))

colnames(counts) <- sub(".*/","",sub("_L00[1,2]_sort.*","",
                                     sub("[1,2]_[1].*_.*DXX_","", 
                                         sub("_dual.*","",salmon))))

#' Combine the technical replicates
counts <- sapply(split.data.frame(t(as.data.frame(counts)),colnames(counts)),colSums)

#' ## stage specificity
load(here("data/analysis/DE/vst-aware.rda"))
samples_m <- read.csv(here("doc/samples.csv"))
time_expression <- expressionSpecificity(
  exp.mat = vsta,
  tissues = as.character(samples_m$Stages[match(colnames(vsta),samples_m$SampleID)]),
  output = "complete")

#' Visualisation
#' 
#' * expression is either very stage specific (peak at 1) or spreading between 0 and 0.6 (rather ubiquitous) 
plot(density(time_expression[!is.na(time_expression[,"score"]),"score"]),
     main="time specificity distribution")

#' * the expression for stage specific gene is on average low
boxplot(time_expression[time_expression[,"score"]==1,"maxn"],
        main="expression of the stage specific genes",ylab="vst expression")

#' * the expression of gene with a stage specificity >= 0.9 is as expected mostly within 
#' a single stage but can spread to as much as all (8) stages
boxplot(time_expression[time_expression[,"score"]>=0.9,"n"])

#' here is the example
time_expression[time_expression[,"score"]>=0.9 & time_expression[,"n"]==8,]

#' Some summary statistics 
means <- rowMeans2(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]),na.rm = TRUE)
median <- rowMedians(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]))
sd <- rowSds(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]))
mad <- rowMads(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]))
boxplot(list(mean=means,sd=sd,median=median,mad=mad))

#' Start a tibble to store the metadata
metadata <- left_join(as.data.frame(counts) %>% rownames_to_column("TRINITY_ID"),
                      as.data.frame(time_expression) %>% rownames_to_column("TRINITY_ID"),
                      by="TRINITY_ID") %>% as_tibble()

#' ## GMAP
#' Process all three file types, concatenate the loci if needed
#' and report per trinity ID, including a type column
metadata %<>% left_join(reduce(mclapply(list.files(here("data/GMAP"),
                             pattern="*_gene",
                             recursive=FALSE,
                             full.names=TRUE),
                  function(f){
                    read_table2(f,
                                col_names = c("Gene","X2","X3",
                                              "start","end","X6","strand","X8","TRINITY_ID"),
                                col_types = cols_only("Gene" = col_character(),
                                                      "start" = col_double(),
                                                      "end" = col_double(),
                                                      "strand" = col_character(),
                                                      "TRINITY_ID" = col_character())) %>%
                      unite(Gene,start,end,strand,col="gmap_loc",sep="_") %>% 
                      mutate(TRINITY_ID=sub(".*=","",TRINITY_ID)) %>% 
                      group_by(TRINITY_ID) %>% 
                      summarise(gmap_loc=paste(gmap_loc,collapse="|")) %>% 
                      mutate(gmap_type=sub("_gene","",basename(f)))
                  },mc.cores=3L),bind_rows),by="TRINITY_ID")


#' ## BedToolsIntersect
metadata %<>% left_join(read_table2(here("data/GMAP/BedToolsIntersect2/GMAP_all-Eugene-gene-only.tsv"),
                          col_names = c("scaffold","X2","X3",
                                        "X4","X5","X6","X7","X8","TRINITY_ID",
                                        "X10","X11","X12","start","end","X15",
                                        "strand","X17","GENE_ID","intersect_length"),
                          col_types = cols_only(scaffold=col_character(),
                                                "TRINITY_ID" = col_character(),
                                                "start" = col_double(),
                                                "end" = col_double(),
                                                "strand" = col_character(),
                                                "GENE_ID" = col_character(),
                                                "intersect_length" = col_double())) %>% 
  mutate(gene_intersect_percent=round(intersect_length/(end-start+1)*100,digits=2)) %>% 
  unite(scaffold,start,end,strand,col="intersect_loc",sep="_") %>% 
  mutate(TRINITY_ID=sub(".*=","",TRINITY_ID),
         GENE_ID=sub(".*=","",GENE_ID))%>% 
  group_by(TRINITY_ID) %>% 
  summarise_all(paste,collapse="|"),by="TRINITY_ID")


#' ## BedToolsSubtract
metadata %<>% left_join(read_table2(here("data/GMAP/BedToolsSubtract/GMAP_all-Eugene-gene-only_no-self_no-full-overlap.gff3"),
                        col_names = c("scaffold","X2","X3",
                                      "t_start","t_end","X6","X7","X8","TRINITY_ID",
                                      "X10","start","end",
                                      "X13","X14",
                                      "strand","subtract_length"),
                        col_types = cols_only("scaffold" = col_character(),
                                              "t_start" = col_double(),
                                              "t_end" = col_double(),
                                              "TRINITY_ID" = col_character(),
                                              "start" = col_double(),
                                              "end" = col_double(),
                                              "strand" = col_character(),
                                              "subtract_length" = col_double())) %>% 
  mutate(trinity_subtract_percent=round(subtract_length/(t_end-t_start+1)*100,digits=2)) %>% 
  unite(scaffold,start,end,strand,col="subtract_loc",sep="_") %>% 
  mutate(TRINITY_ID=sub(".*=","",TRINITY_ID))%>% 
  group_by(TRINITY_ID) %>% 
  summarise_all(paste,collapse="|"),by="TRINITY_ID")
  
#' There are about 200 thousands without any overlap, 
#' while there are trinity IDs that intersect a gene than
#' can be subtracted (i.e. some are wholly contained).
table(subtract=!is.na(metadata$trinity_subtract_percent),
      intersect=!is.na(metadata$gene_intersect_percent))

#' ## trmap
metadata %<>% left_join(reduce(mclapply(list.files(here("data/trmap/trinity"),
                           full.names=TRUE,
                           pattern="trmapIDs.txt",
                           recursive=TRUE),
                function(f){
                  tibble(TRINITY_ID=scan(f,what="character"),
                         trmap_from=basename(dirname(f)))
                },mc.cores=3L),bind_rows) %>% distinct(),by="TRINITY_ID")

#' All the trmap are also contained in the bedtools intersect,
#' however there are about 50,000 more that are not identified by trmap,
#' because the intersect is probably too small
table(trmap=!is.na(metadata$trmap_from),
      intersect=!is.na(metadata$gene_intersect_percent))

plot(density(as.numeric(
  unlist(strsplit(metadata$gene_intersect_percent[
    is.na(metadata$trmap_from) & 
      !is.na(metadata$gene_intersect_percent)],"\\|")))),
  col=3,main="distribution of the gene intersection length")
lines(density(as.numeric(
  unlist(strsplit(metadata$gene_intersect_percent[
    !is.na(metadata$trmap_from) & 
      !is.na(metadata$gene_intersect_percent)],"\\|")))),col=4)
legend("topright",lty=1,col=c(3,4),c("trmap only","both"))

#' ## Diamond
metadata %<>% left_join(read_table2(here("data/DIAMOND/uniref90.dmnd_Trinity.blt.gz"),
                       col_names = c("TRINITY_ID", "reference_ID", 
                                     "identity", "alignment_length", "mismatch",
                                     "gap_open", "start_trinity", "end_trinity", 
                                     "start_ref", "end_ref",
                                     "evalue", "bitscore", "trinity_length", 
                                     "ref_length", "taxonomy")) %>% 
  mutate(trinity_cov=alignment_length/trinity_length,
         ref_cov=alignment_length/ref_length,
         TRINITY_ID=sub("\\.p[0-9]+","",TRINITY_ID)) %>% 
  group_by(TRINITY_ID) %>% 
  arrange(desc(identity),desc(trinity_cov),desc(ref_cov)) %>% 
  slice(1) %>% ungroup(),by="TRINITY_ID")

#' ## Taxonomy
#' We extract the UniRef ID
uniref_id <- unique(metadata$reference_ID[!is.na(metadata$reference_ID)])

#' Create a function
f<-function(lines,pos){
    return(lines[lines$reference_ID %in% uniref_id,])
}

#' To parse the input data on the fly
#' (if the txt file is compressed, decompress it prior to 
#' running that line, as oytherwise R will do the decompression in mem, 
#' loosing the advantage of chunking)
tax_map <- read_delim_chunked(col_names=c("reference_ID","taxon","taxID"),
                            col_types=c("ccc"),
                            here("uniref/annotation/uniref90_id-table.txt"),
                            chunk_size=1e6,delim=" ",callback=DataFrameCallback$new(f))

#' A quick look at the top species
mar <- par("mar")
par(mar=c(10.1,4.1,0.1,0.1))
barplot(sort(table(tax_map$taxon),decreasing=TRUE)[1:20],las=2)
par(mar=mar)

#' Extend the metadata 
metadata %<>% left_join(tax_map,by="reference_ID")

#' ## Sequences
ref <- readDNAStringSet(here("data/trinity/Trinity.fasta"))
metadata %<>% left_join(tibble(TRINITY_ID=sub(" .*","",names(ref)),
                               GC=rowSums(alphabetFrequency(ref)[,c("C","G")]) / width(ref),
                               length=width(ref)),by="TRINITY_ID")

#' ## CNCI
metadata %<>% left_join(read_tsv(file = here("data/CNCI/CNCI.index"),
                                 col_names = c("TRINITY_ID", "CNCI_index", "CNCI_score", "CNCI_start", "CNCI_end", "CNCI_length"),
                                 col_types=cols(.default=col_character())) %>%
                          mutate(TRINITY_ID=gsub(">| .*","",TRINITY_ID)),by="TRINITY_ID")

#' ## PLEK
metadata %<>% left_join(read_tsv(file=here("data/PLEK/Trinity.txt"),
                                 col_names = c("PLEK_type","PLEK_score","TRINITY_ID"),
                                 col_types = c("cdc")) %>% 
                          mutate(TRINITY_ID=gsub(">| .*","",TRINITY_ID)),by="TRINITY_ID")

#' ## CPC2
metadata %<>% left_join(read_tsv(file = here("data/CPC2/results.txt"),
                                 comment="#",
                                 col_names = c("TRINITY_ID",
                                               "CPC2_transcript_length",
                                               "CPC2_peptide_length",
                                               "CPC2_Fickett_score",
                                               "CPC2_pI",
                                               "CPC2_ORF_integrity",
                                               "CPC2_coding_probability" ,
                                               "CPC2_label"),
                                 col_types=c("cddddddc")),by="TRINITY_ID")

#' ## PLncPRO
metadata %<>% left_join(read_table2(file = here("data/PLncPRO/Trinity.txt"),
                                    col_names = c("TRINITY_ID","PLncPRO_coding_potential","X3","X4"),
                                    col_types = cols_only("TRINITY_ID"=col_character(),
                                                          "PLncPRO_coding_potential"=col_number())),by="TRINITY_ID")

#' ## Transdecoder
metadata %<>% left_join(read_table2(here("data/Transdecoder/Trinity.fasta.transdecoder.IDs.txt"),
                                    col_names = c("X1","X2","X3",
                                                  "Transdecoder_type","Transdecoder_AA_length","Transdecoder_score","Coord"),
                                    col_types = cols_only("Transdecoder_type"=col_character(),
                                                          "Transdecoder_AA_length"=col_number(),
                                                          "Transdecoder_score"=col_character(),
                                                          "Coord"=col_character())) %>% 
                          separate(Coord,into=c("TRINITY_ID","Transdecoder_start","Transdecoder_end"),extra="drop",sep=":|=|-|\\(") %>% 
                          mutate(Transdecoder_type=factor(sub("type:","",Transdecoder_type)),
                                 Transdecoder_strand=gsub("\\(|\\).*","",Transdecoder_score),
                                 Transdecoder_score=parse_double(sub(".*=","",Transdecoder_score))),by="TRINITY_ID")

#' # Flagging
#' * Coding
metadata %<>% mutate(coding=length >= 200 & 
                       CNCI_index == "coding" & 
                       PLEK_type == "Coding" & 
                       CPC2_label == "coding" & 
                       PLncPRO_coding_potential == 1 & 
                       !is.na(Transdecoder_type))

#' * Non-coding
metadata %<>% mutate(non_coding=length >= 200 & 
                       CNCI_index == "noncoding" & 
                       PLEK_type == "Non-coding" & 
                       CPC2_label == "noncoding" & 
                       PLncPRO_coding_potential == 0)

#' * Network
metadata %<>% mutate(seidr=TRINITY_ID %in% scan(here("data/seidr/genes.tsv"),sep="\t",what="character"))

#' # Export
write_tsv(metadata,file=here("data/metadata.tsv"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

