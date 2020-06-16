suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
  library(here)
  library(magrittr)
  })
dir(here("data/"))
ref <- readDNAStringSet("data/trinity/Trinity.fasta")
bla <- sub(" .*","",names(ref))
bla_new <- alphabetFrequency(ref)
bla_GC <- rowSums(alphabetFrequency(ref)[,c("C","G")]) / width(ref)
length_ref <- width(ref)
#' create the big tibble
CNCI_read <- read_delim(file = "data/CNCI/CNCI.index",
                        delim="\t",
                        skip = 1,
                        col_names = c("Transcript.ID", "index", "score", "start", "end", "length")) %>%
    mutate(Transcript.ID,Transcript.ID=gsub(">| .*","",Transcript.ID))
CNCI_tib <- CNCI_read$index

PLEK_read <- read_delim(file="data/PLEK/Trinity.txt",
                        delim="\t",
                        skip = 1,
                        col_names = c("Type","Score","Transcript.ID")) %>% 
    mutate(Transcript.ID,Transcript.ID=gsub(">| .*","",Transcript.ID))
PLEK_index <- left_join(CNCI_read, PLEK_read, by = NULL, copy=FALSE)
PLEK_tib <- PLEK_index$index

CPC2_read <- read_delim(file = "data/CPC2/results.txt",
                        delim="\t")
CPC2_tib <- CPC2_read$label

PLncPRO_read <- read_table2(file = "data/PLncPRO/Trinity.txt",
                            col_names = c("Transcript.ID","coding_potential","X3","X4"),
                            col_types = cols_only("Transcript.ID"=col_character(),
                                                  "coding_potential"=col_number())) %>%
  mutate(Transcript.ID,Transcript.ID=gsub(">| .*","",Transcript.ID))

PLncPRO_index <- left_join(CNCI_read, PLncPRO_read, by = NULL, copy=FALSE)

PLncPRO_tib <- PLncPRO_index$coding_potential

transdecoder <- read_table2("data/Transdecoder/Trinity.fasta.transdecoder.IDs.txt",
                            col_names = c("Transdecoder.ID","X2","X3",
                                          "Type","AALength","Score","Coord"),
                            col_types = cols_only("Transdecoder.ID"=col_character(),
                                                  "Type"=col_character(),
                                                  "AALength"=col_number(),
                                                  "Score"=col_character(),
                                                  "Coord"=col_character())) %>% 
  mutate(Transdecoder.ID,Transdecoder.ID=gsub(">| .*","",Transdecoder.ID)) %>%  
  mutate(Transdecoder.ID,Transcript.ID=gsub("\\.p\\d+","",Transdecoder.ID)) %>% 
  mutate(Type,Type=factor(sub("type:","",Type))) %>% 
  mutate(Score,Score=parse_double(sub(".*=","",Score),locale=locale(decimal_mark = "."))) %>% 
  mutate(Coord,Strand=factor(gsub(".*\\(|\\)","",Coord))) %>% 
  mutate(Coord,Start=parse_integer(gsub(".*:|-.*","",Coord))) %>% 
  mutate(Coord,End=parse_integer(gsub(".*[0-9]+-|\\(.*","",Coord)))
transdecoder_index <- left_join(CNCI_read, transdecoder, by = NULL, copy=FALSE)
transdecoder_tib <- transdecoder_index$Type

final_tibble <- tibble(Transcript.ID = bla, GC_content = bla_GC, length = length_ref, CNCI = CNCI_tib, PLEK = PLEK_tib, CPC2 = CPC2_tib, PLncPRO = PLncPRO_tib, Transdecoder = transdecoder_tib)
#remove the duplicated transcript detected by Salmon before to continue
final_tibble_bla <-  final_tibble[-c(164348),]

#adding trmap

trmap_uniq <- read_delim("data/trmap/trinity/uniq/trmapIDs.txt",
                         delim = "/",
                         col_names = "Transcript.ID")

trmap_mult <- read_delim("data/trmap/trinity/mult/trmapIDs.txt",
                         delim = "/",
                         col_names = "Transcript.ID")

trmap_transloc <- read_delim("data/trmap/trinity/transloc/trmapIDs.txt",
                             delim = "/",
                             col_names = "Transcript.ID")


blabla <- anti_join(final_tibble_bla, trmap_mult, by = NULL, copy=FALSE)
blabla2 <- anti_join(blabla, trmap_uniq, by = NULL, copy=FALSE)
final_tibble <- anti_join(blabla2, trmap_transloc, by = NULL, copy=FALSE)


# trying to add time_expression
load(here("data/analysis/DE/vst-aware.rda"))
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
samples_m <- read.csv("doc/samples.csv")
time_expression <- expressionSpecificity(exp.mat = vsta[,samples_m$SampleID],
                                         tissues = as.character(samples_m$Stages),
                                         output = "complete")

stuff <- as_tibble(time_expression)
stuff_new <- stuff %>% add_column(Transcript.ID = final_tibble_bla$Transcript.ID)
final_tibble_s <- left_join(final_tibble, stuff_new, by = NULL, copy=FALSE)

#filter only coding
coding <- final_tibble_s %>% filter(length >= 200, CNCI == "coding", PLEK == "coding", CPC2 == "coding", PLncPRO == "1", Transdecoder != "NA") %>% 
  rename_if(grepl("aij",colnames(.)),funs(str_replace(.,"aij\\.","")))

#filter only non-coding

non_coding <- final_tibble_s %>% filter(length >= 200, CNCI == "noncoding", PLEK == "noncoding", CPC2 == "noncoding", PLncPRO == "0") %>% 
  rename_if(grepl("aij",colnames(.)),funs(str_replace(.,"aij\\.","")))


#library(Biostrings)
#seq <- readDNAStringSet("data/trinity/Trinity.fasta")
#names(seq) <- sub(" .*","",names(seq))
#IDs <- non_coding$TRINITY_ID
#writeXStringSet(seq[IDs],file="/mnt/picea/projects/spruce/nstreet/spruce-lncRNA-network/fasta/lncRNA.fasta")



#time expression for only non-coding

time_expression_nc <- non_coding %>% select(Transcript.ID, score, S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, maxn, n, peak)
time_expression_nc_filtered <-  time_expression_nc %>% filter(score != "NA")

noNaN_batch1 <- time_expression_nc_filtered %>% filter(S0 != "NaN",S9 != "NaN")
no_batch1 <- anti_join(time_expression_nc_filtered,noNaN_batch1, by = NULL, copy=FALSE)
time_expression_nc_last <- no_batch1 %>% select(-S0,-S9)


write_delim(x = time_expression_nc_last, path = "doc/time_expression_nc_last.tsv", delim = " ")


par(bg="orange")
plot(density(time_expression_nc_last$score), 
     xlab="stage_specificity",
     ylab="density",
     main="NON_CODING",
     lwd=3,
     font.main=2,
     font.lab=2,
     cex.main=3,
     cex.lab=2.5)

barplot(table(time_expression_nc_last$peak))
barplot(table(time_expression_nc_last %>% filter(score > 0.9) %>% select(peak)))
time_expression_nc_last %<>% left_join(tibble(Transcript.ID=rownames(vsta),avgexp=rowMaxs(vsta)))
plot(density(as.matrix(time_expression_nc_last %>% filter(score > 0.9) %>% select(avgexp))))
plot(density(as.matrix(time_expression_nc_last %>% filter(score <= 0.9) %>% select(avgexp))))


#time expression for only coding
time_expression_c <- coding %>% select(Transcript.ID, score, S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, maxn, n, peak)
time_expression_c_filtered <-  time_expression_c %>% filter(score != "NA")

noNaN_batch1 <- time_expression_c_filtered %>% filter(S0 != "NaN",S9 != "NaN")
no_batch1 <- anti_join(time_expression_c_filtered,noNaN_batch1, by = NULL, copy=FALSE)
time_expression_c_last <- no_batch1 %>% select(-S0,-S9)

write_delim(x = time_expression_c_last, path = "doc/time_expression_c_last.tsv", delim = " ")


par(bg="plum")
plot(density(time_expression_c_filtered$score), 
     xlab="stage_specificity",
     ylab="density",
     main="CODING",
     lwd=3,
     font.main=2,
     font.lab=2,
     cex.main=3,
     cex.lab=2.5)
barplot(table(time_expression_c_last$peak))
barplot(table(time_expression_c_last %>% filter(score > 0.9) %>% select(peak)))

time_expression_c_last %<>% left_join(tibble(Transcript.ID=rownames(vsta),avgexp=rowMaxs(vsta)))
plot(density(as.matrix(time_expression_c_last %>% filter(score > 0.9) %>% select(avgexp))))
plot(density(as.matrix(time_expression_c_last %>% filter(score <= 0.9 ) %>% select(avgexp))))

par(bg="white")
boxplot(list(non_coding_specific=as.matrix(time_expression_nc_last %>% filter(score > 0.9) %>% select(avgexp)),
  non_coding_aspecific=as.matrix(time_expression_nc_last %>% filter(score <= 0.9) %>% select(avgexp)),
  coding_specific=as.matrix(time_expression_c_last %>% filter(score > 0.9) %>% select(avgexp)),
  coding_aspecific=as.matrix(time_expression_c_last %>% filter(score <= 0.9 ) %>% select(avgexp))),log="y",
  ylab="expression",
  col=c(rep("orange",1,2),rep("plum",3,4)),
  font.lab=2,
  cex.lab=1.5)



#checking GC_content

GC_50 <- non_coding %>% filter(GC_content >= 0.50)
GC_40 <- non_coding %>% filter(GC_content >= 0.40)
GC_less_30 <- non_coding %>% filter(GC_content < 0.30)
GC_3040 <- non_coding %>% filter(GC_content >= 0.30, GC_content < 0.40)

GC_50c <- coding %>% filter(GC_content >= 0.50)
GC_40c <- coding %>% filter(GC_content >= 0.40)
GC_less_30c <- coding %>% filter(GC_content < 0.30)
GC_3040c <- coding %>% filter(GC_content >= 0.30, GC_content < 0.40)

