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
  library(ggplot2)
  library(ggforce)
  library(purrr)
})
  
#' Helper functions
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
source(here("UPSCb-common/src/R/blastUtilities.R"))

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
load(here("data/analysis/DE/vst-aware_linc.rda"))
samples_m <- read.csv(here("doc/samples_final.csv"))
time_expression <- expressionSpecificity(
  exp.mat = vsta,
  tissues = as.character(samples_m$Stages[match(colnames(vsta),samples_m$SampleID)]),
  output = "complete")

colnames(time_expression) <- gsub("aij.","",colnames(time_expression))

vsta <- sapply(split.data.frame(t(as.data.frame(vsta)),colnames(vsta)),colSums)

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
#means <- rowMeans2(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]),na.rm = TRUE)
#median <- rowMedians(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]))
#sd <- rowSds(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]))
#mad <- rowMads(as.matrix(time_expression[,grep("aij",colnames(time_expression),value=TRUE)]))
#boxplot(list(mean=means,sd=sd,median=median,mad=mad))

#' Start a tibble to store the metadata
colnames(counts) <- paste0("raw_",colnames(counts))
colnames(vsta) <- paste0("vst_",colnames(vsta))
metadata <- left_join(left_join(as.data.frame(counts) %>% rownames_to_column("TRINITY_ID"),
                                as.data.frame(vsta) %>% rownames_to_column("TRINITY_ID"),by="TRINITY_ID"),
                      as.data.frame(time_expression) %>% rownames_to_column("TRINITY_ID"),
                      by="TRINITY_ID") %>% as_tibble()

#' ## GMAP
#' Process all three file types, concatenate the loci if needed
#' and report per trinity ID, including a type column
metadata %<>% left_join(purrr::reduce(mclapply(list.files(here("data/GMAP"),
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

#' Process the original files to extract the number of exons per alignments
metadata %<>% left_join(
  purrr::reduce(mclapply(list.files(here("data/GMAP"),
                             pattern="Pabies1\\.0-Trinity\\.2\\..*\\.gz",
                             recursive=FALSE,
                             full.names=TRUE),
                  function(f){read_tsv(f,
                                       col_names = c("X1","X2","Type",
                                                     "X4","X5","X6","X7","X8","TRINITY_ID"),
                                       col_types = cols_only("Type" = col_factor(),
                                                             "TRINITY_ID" = col_character()),comment="#") %>% 
                      filter(Type=="exon") %>% 
                      separate(TRINITY_ID,into=paste0("V",1:11),sep=";|=| ") %>% 
                      select(V6,TRINITY_ID=V8) %>% group_by(V6) %>% add_count(name="Exon") %>% 
                      group_by(TRINITY_ID) %>% select(TRINITY_ID,Exon) %>% summarise_all(paste,collapse="|") 
                  },mc.cores=3L),bind_rows),by="TRINITY_ID")

#' ##BedtoolsClosest
metadata %<>% left_join(read_table2(here("data/GMAP/BedToolsClosest/genes_all_sorted-ref_genes_all.tsv.gz"),
            col_names = c("scaffold","X2","X3",
                          "t_start","t_end","X6","strand","X8","ID",
                          "X10","X11","X12","start","end",
                          "X15","X16","X17",
                          "GENE_ID","closest_length"),
            col_types = cols_only("scaffold" = col_character(),
                                  "t_start" = col_double(),
                                  "t_end" = col_double(),
                                  "strand" = col_character(),
                                  "ID" = col_character(),
                                  "start" = col_double(),
                                  "end" = col_double(),
                                  "GENE_ID" = col_character(),
                                  "closest_length" = col_double())) %>% 
  mutate(trinity_closest_percent=round(closest_length/(t_end-t_start+1)*100,digits=2)) %>% 
  separate(ID,into=c("transcript","TRINITY_ID"), sep=";", convert = TRUE) %>% 
    mutate(TRINITY_ID=sub(".*=","",TRINITY_ID),
           GENE_ID=sub(".*=","",GENE_ID))%>% 
  select(-transcript) %>% 
  group_by(TRINITY_ID) %>% 
  summarise_all(paste,collapse="|"),by="TRINITY_ID")

#only_nc_BC <- metadata %>% filter(non_coding == TRUE &
  #                      (!is.na(GENE_ID)) &
  #                     (closest_length != 0) &
  #                     (closest_length != -1) &
  #                     (!is.na(score))) 
#only_nc_BC_500 <- only_nc_BC %>% filter(closest_length > 500)
#only_nc_BC_1000 <- only_nc_BC %>% filter(closest_length > 1000)

#ggplot(only_nc_BC, aes(x = closest_length)) +
  #geom_histogram(fill = "orange") + 
  #scale_x_continuous(trans='log10') +
  #labs(x = "distance_log10", y = "Frequency") +
  #theme_bw()

#ggplot(only_nc_BC, aes(x = closest_length)) +
  #geom_density(fill = "orange") +
  #scale_x_continuous(trans='log10') +
  #labs(x = "distance_log10", y = "Frequency") +
  #theme_classic() 

#' ## Diamond
metadata %<>% left_join(read_table(here("data/DIAMOND/uniref90.dmnd_Trinity.blt.gz"),
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
ref <- readDNAStringSet(here("data/trinity/Trinity.fasta.gz"))
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
metadata %<>% left_join(read_table(here("data/Transdecoder/Trinity.fasta.transdecoder.IDs.txt"),
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

#" ## Cd-hist-est
metadata %<>% left_join(read_tsv(here("data/analysis/cdhit/cd-hit-est_id80_cluster-membership.tsv"),
                                 show_col_types = FALSE) %>% rename(cdhit_cluster=cluster),by=c("TRINITY_ID"="ID"))

#' ## PLAZA
#' We may have to re-run the Blast at some point as not all IDs are well formatted
metadata %<>% left_join(readBlast(here("data/blastn_plaza/plaza_linc_network.blt"),
                                  plot=FALSE,verbose=FALSE,
                 format=BM8ext)$df %>% separate(col=subject.id,
                                 into=c("Plaza_scaffold","Plaza_species"),
                                 sep="\\|",fill="right",extra="drop") %>% 
                   rename_with(.fn=function(x){paste0("Plaza_",x)}) %>% 
  dplyr::rename(TRINITY_ID=Plaza_query.id) %>% group_by(TRINITY_ID) %>%
  summarise_all(paste,collapse="|"),by="TRINITY_ID")

#' ## miRNA
metadata %<>% left_join(readBlast(here("precursors/linc_network.fasta_Pabies_SE_miRNA.precursor.blt"),
          format=BM8ext,verbose=FALSE,plot=FALSE)$df %>% 
  rename_with(.fn=function(x){paste0("miRNA_",x)}) %>% 
  dplyr::rename(TRINITY_ID=miRNA_subject.id) %>% group_by(TRINITY_ID) %>%
  summarise_all(paste,collapse="|") ,by="TRINITY_ID")

#' ## Infomap
metadata %<>% left_join(read_tsv(here("data/seidr/backbone/infomapClusters.tsv"),
                                 show_col_types = FALSE) %>% 
  dplyr::rename(TRINITY_ID=gene,Infomap_cluster=cluster),by="TRINITY_ID")

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
metadata %<>% mutate(seidr=TRINITY_ID %in% scan(here("data/analysis/seidr/genes.tsv"),sep="\t",what="character"))


#' * lincRNAs

metadata %<>% mutate(lincRNAs= non_coding == TRUE &
                     !is.na(score) &
                     as.integer(sapply(strsplit(metadata$closest_length,"\\|"),min)) > 1000) 

lincRNAs <- metadata %>% filter(non_coding == TRUE &
                                !is.na(score) &
                                as.integer(sapply(strsplit(metadata$closest_length,"\\|"),min)) > 1000) 

#write_tsv(lincRNAs,file=here("doc/lincRNAs.tsv"))

#' # Export
write_tsv(metadata,file=here("data/metadata.tsv.gz"))
saveRDS(metadata,file=here("data/metadata.rds"))

#exclude batch 1, not used in the analyses
meta <- read_tsv(here("data/metadata.tsv.gz"))
meta <- select(me, -c(2,3,4,5,6,7,8,9,10,11,vst_P464_202,vst_P464_203,vst_P464_204,vst_P464_205,vst_P464_206,vst_P464_207B,vst_P464_208,vst_P464_209B))
write_tsv(meta,file=here("data/metadata_new.tsv.gz"))
saveRDS(meta,file=here("data/metadata_new.rds"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
