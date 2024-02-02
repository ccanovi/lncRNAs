
library(here)
library(tidyverse)

source(here("UPSCb-common/src/R/gopher.R"))
pal12 <- brewer.pal(12,"Paired")
pal_grey <- brewer.pal(9,"Greys")
# First we declare where the file is
statsFile <- here("data/seidr/backbone/statistics.tsv")

# We read the file as a table
stats <- read.table(statsFile, header = T, stringsAsFactors = F, sep='\t')

# Create two data frames, one with all source node information and one with target node information
# We need to do it like this because most of the genes will appear as source AND target but a 
# few genes might appear only as source OR as target.
statsDF1 <- data.frame(gene = stats$Source,
                       PageRank = stats$PageRank_Source,
                       Betweenness = stats$Betweenness_Source,
                       Strength = stats$Strength_Source,
                       Eigenvector = stats$Eigenvector_Source,
                       Katz = stats$Katz_Source)
statsDF2 <- data.frame(gene = stats$Target,
                       PageRank = stats$PageRank_Target,
                       Betweenness = stats$Betweenness_Target,
                       Strength = stats$Strength_Target,
                       Eigenvector = stats$Eigenvector_Target,
                       Katz = stats$Katz_Target)

# Bind both data frames
statsDF <- rbind(statsDF1, statsDF2)

# Remove duplicated rows
stats.res <- statsDF[!duplicated(statsDF),]


# Using tidyverse and dividing by RNA classes
total_tibble <- as_tibble(stats.res)

lincRNAs <- total_tibble %>% filter(grepl("TRINITY",gene)) 
coding <- total_tibble %>% filter(grepl("MA_",gene)) 
genes <- coding$gene
miRNAs <- total_tibble %>% filter(grepl("miRNA",gene))

lincRNAs <- lincRNAs %>% add_column(type = "lincRNAs")
miRNAs <- miRNAs %>% add_column(type = "miRNAs")

TF <- read.csv(file = "doc/TF_list.txt",sep = "")
TF_list <- TF$TF_ID

library(VennDiagram)
library(ggplot2)
oveer <- calculate.overlap(list(coding=coding$gene,TF=TF_list,filename=NULL))
TF_last <- oveer$a2
TF_last <- as_tibble(TF_last)
TF_last <- TF_last %>% add_column(type = "TF")
TF_last <- rename(TF_last, gene = value)
TFs <- semi_join(coding,TF_last,by = NULL, copy=FALSE)
TFs <- TFs %>% add_column(type = "TFs")
coding <- anti_join(coding,TF_last,by = NULL, copy=FALSE)
coding <- coding %>% add_column(type = "coding")

test <- rbind(lincRNAs,miRNAs)
test2 <- rbind(test,TFs)
test_last <-  rbind(test2,coding)

#Violin plot based on network stats
library(ggpubr)
my_comparisons <- list(c("coding","lincRNAs"),c("lincRNAs","miRNAs"),c("miRNAs","TFs"),c("TFs","coding"))
p <-  ggplot(test_last, aes(x=type,y=PageRank,colour=type)) +  
  geom_violin(size=1.5) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=13)) +
  theme(legend.position = "None") +
  scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999")) +
  stat_compare_means(comparison = my_comparisons) +
  stat_compare_means(label.y.npc = "top",label.x.npc = "centre")
p

b <-  ggplot(test_last, aes(x=type,y=Betweenness,colour=type)) +
  geom_violin(size=1.5) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=13)) +
  theme(legend.position = "None") +
  scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999")) +
  stat_compare_means(comparison = my_comparisons) +
  stat_compare_means(method = "kruskal.test",label.y.npc = "top",label.x.npc = "centre")
b

s <-  ggplot(test_last, aes(x=type,y=Strength,colour=type)) +
  geom_violin(size=1.5) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=13)) +
  theme(legend.position = "None") +
  scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999")) +
  stat_compare_means(comparison = my_comparisons) +
  stat_compare_means(method = "kruskal.test",label.y.npc = "top",label.x.npc = "centre")

s

k <-  ggplot(test_last, aes(x=type,y=Katz,colour=type)) +
  geom_violin(size=1.5) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=13)) +
  theme(legend.position = "None") +
  scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999")) +
  stat_compare_means(comparison = my_comparisons) +
  stat_compare_means(method = "kruskal.test",label.y.npc = "top",label.x.npc = "centre")
k

e <-  ggplot(test_last, aes(x=type,y=Eigenvector)) +
  geom_violin() +
  theme_classic()
e

library(cowplot)
all <- plot_grid(p,b,s,labels=c("A","B","C"),vjust = 1.3,ncol=3)#hjust=-20,)
plot(all)
ggsave(filename=here("data/analysis/figures//Network_characteristics_by_col_col.png"),device ="png",dpi = 600)

source(here("UPSCb-common/Rtoolbox/R/plotEigengene.R"))
# Get first degree neighbours of the first 5 PageRank lincRNAs. Remember to run plotEigengene.R before.
load(here("data/analysis/seidr/network.rda")) #or to be seen
samples <- read_csv(here("doc/samples_B2.csv")) #or B2
edgeList <- read.table(here("data/seidr/backbone/edgeList.txt"))

getGeneFDN <- function(edgeList, gene, source.col=1, target.col=2) {
  # TODO check for data type
  # TODO check that all genes are in the edgelist
  
  s2t <- edgeList[edgeList[source.col] == gene,][,target.col]
  t2s <- edgeList[edgeList[target.col] == gene,][,source.col]
  res <- union(s2t,t2s)
  
  return(res)
}

#run plotEigengene.R in src/R to plot profiles
FDG_c1 <- getGeneFDN(edgeList,"TRINITY_DN96496_c0_g1_i1")
lincRNA_c1 <-  FDG_c1[grep("TRINITY",FDG_c1)]
FDG_c1[grep("MA",FDG_c1)]
c1 <- plotEigengene(t(dat),"TRINITY_DN96496_c0_g1_i1",rep("stages",nrow(t(dat))),samples$Stages)
plot(c1)
ggsave(filename=here("data/analysis/figures//candidate1_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(FDG_c1, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes <- unlist(SDN_list) %>% unique()

FDG_c1_enr <- gopher(FDG_c1, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")

FDG_c2 <- getGeneFDN(edgeList,"TRINITY_DN12829_c0_g1_i5")
lincRNA_c2 <-  FDG_c2[grep("TRINITY",FDG_c2)]
FDG_c2[grep("MA",FDG_c2)]
c2 <- plotEigengene(t(dat), "TRINITY_DN12829_c0_g1_i5",rep("stages", nrow(t(dat))),samples$Stages)
plot(c2)
ggsave(filename=here("data/analysis/figures//candidate2_PR.png"),device ="png",dpi = 600)

SDN_list2 <- lapply(FDG_c2, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes2 <- unlist(SDN_list2) %>% unique()

FDG_c2_enr <- gopher(FDG_c2, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")

FDG_c3 <- getGeneFDN(edgeList,"TRINITY_DN52747_c1_g1_i1")
lincRNA_c3 <-  FDG_c3[grep("TRINITY",FDG_c3)]
FDG_c3[grep("MA",FDG_c3)]
FDG_c3[grep("miRNA",FDG_c3)]
c3 <- plotEigengene(t(dat), "TRINITY_DN52747_c1_g1_i1",rep("stages", nrow(t(dat))),samples$Stages)
plot(c3)
ggsave(filename=here("data/analysis/figures//candidate3_PR.png"),device ="png",dpi = 600)

SDN_list3 <- lapply(FDG_c3, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes3 <- unlist(SDN_list3) %>% unique()

FDG_c3_enr <- gopher(FDG_c3, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")

FDG_c4 <- getGeneFDN(edgeList,"TRINITY_DN6985_c1_g1_i1")
lincRNA_c4 <-  FDG_c4[grep("TRINITY",FDG_c4)]
FDG_c4[grep("MA",FDG_c4)]
c4 <- plotEigengene(t(dat), "TRINITY_DN6985_c1_g1_i1",rep("stages", nrow(t(dat))),samples$Stages)
plot(c4)
ggsave(filename=here("data/analysis/figures//candidate4_PR.png"),device ="png",dpi = 600)

SDN_list4 <- lapply(FDG_c4, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes4 <- unlist(SDN_list) %>% unique()

FDG_c4_enr <- gopher(FDG_c4, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")

FDG_c5 <- getGeneFDN(edgeList,"TRINITY_DN20571_c0_g1_i1")
lincRNA_c5 <-  FDG_c5[grep("TRINITY",FDG_c5)]
FDG_c5[grep("MA",FDG_c5)]
FDG_c5[grep("miRNA",FDG_c5)]
c5 <- plotEigengene(t(dat), "TRINITY_DN20571_c0_g1_i1",rep("stages", nrow(t(dat))),samples$Stages)
plot(c5)
ggsave(filename=here("data/analysis/figures//candidate5_PR.png"),device ="png",dpi = 600)

SDN_list5 <- lapply(FDG_c5, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes5 <- unlist(SDN_list) %>% unique()

FDG_c5_enr <- gopher(FDG_c5, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")

# enrichment
FDG_enr5 <- gopher(genes=FDG_c5, alpha = 0.05, task=list("go", "mapman","kegg","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
FDG_enr5_2nd <- gopher(genes=SDN_genes5, alpha = 0.05, task=list("go", "mapman","kegg","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
plotEnrichedTreemap(x = FDG_c5, enrichment = "mapman") #,namespace = "BP")

#adding newGOA

source(here("~/Git/aspleaf_candidates/lncRNA_pipeline_smk/UPSCb-common/src/R/percentile.R"))

NewGOA_linc <- read_delim("functional_prediction/results/lincAnnotation.txt", delim = ",",
                          col_names = c("gene","GOterm","OriginalAnnotation","PredictionScore"))

NewGOA_genes <- read_delim("functional_prediction/results/geneAnnotation.txt", delim = ",",
                           col_names = c("gene","GOterm","OriginalAnnotation","PredictionScore"))

par(mgp=c(2.3,1,0))
pdf("data/analysis/figures/NewGOA_last.pdf")
png("data/analysis/figures/NewGOA_last.png",width = 600,height = 350)
plot(density(-log2(NewGOA_genes$PredictionScore)), main= NA,ylim=c(0,0.07),col="#F8766D",
     xlab="-log2(PredictionScore)",ylab="density",cex.lab=1.0,lwd=3,bty="l")
title(main="NewGOA selection",adj=0,cex.main= 1.5)
lines(density(-log2(NewGOA_linc$PredictionScore)),col="#00BFC4",lwd=3,bty="l")
#cor.test(percentile(NewGOA_genes$PredictionScore),percentile(NewGOA_linc$PredictionScore))
abline(v=quantile(-log2(NewGOA_genes$PredictionScore),probs=seq(0,.1,0.03)),lty=2,col="#F8766D",lwd=1.0)
abline(v=quantile(-log2(NewGOA_linc$PredictionScore),probs=seq(0,.1,0.03)),lty=2, col="#00BFC4",lwd=1.0)
legend("right", bty = "n",
       legend=c("coding","lincRNA"), border = c("#F8766D","#00BFC4"),
       col= c("#F8766D","#00BFC4"),pch=22,pt.lwd=1.5,
       title = "type",title.adj= 0.15,cex = 1.0)
dev.off()

NewGOA_genes <- NewGOA_genes %>% add_column(type = "coding")
NewGOA_linc <- NewGOA_linc %>% add_column(type = "lincRNAs")
newGOA_total <- rbind(NewGOA_genes,NewGOA_linc)

selected_linc <- subset(NewGOA_linc,PredictionScore >= 2^-13)
selected_genes <- subset(NewGOA_genes,PredictionScore >= 2^-13)
selected_miRNAs <- subset(NewGOA_miRNAs,PredictionScore >= 2^-13)
length(unique(selected_linc$gene))
only_selected_miRNAs <- unique(selected_miRNAs$gene)
only_selected_linc <- unique(selected_linc$gene)
only_selected_genes <- unique(selected_genes$gene)
selected_linc_ord <- selected_linc[order(selected_linc$PredictionScore,decreasing = TRUE), ]
NewGOA_cluster9_genes <-  subset(only_selected_genes, only_selected_genes %in% gene_9$gene)
NewGOA_cluster9_linc <-  subset(only_selected_linc, only_selected_linc %in% lincRNA_9$gene)
NewGOA_cluster7_miRNAs <-  subset(only_selected_miRNAs, only_selected_miRNAs %in% miRNA_7$gene)

#checking how many New_GOA annotation there are in each clusters
NewGOA_miRNAs <- read_delim("functional_prediction/results/miRNA_Annotation.txt", delim = ",",
                            col_names = c("gene","GOterm","OriginalAnnotation","PredictionScore"))
