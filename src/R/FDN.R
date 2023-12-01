
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
tibble <- read_tsv(here("data/metadata.tsv.gz"))
removing <- false_linc_network$TRINITY_ID
save(removing,file=here("doc/extra_filtering.rda"))
stats.res_new <- stats.res[! stats.res$gene %in% removing,]

#trying to do a violin plot based on PageRank
stats.res$gene <- as.factor(stats.res$gene)
stats.res$PageRank <- as.factor(stats.res$PageRank)
geneColors <- rep("#000000", length(unique(test_last$type)))
names(geneColors) <- unique(test_last$type)
geneColors[coding] <-  color="#60B266"
geneColors[lincRNAs] <- "#FEF65B"
geneColors[miRNAs] <-  "#009AA9"

#Then in ggplot
scale_color_manual(values=geneColors) 
p <-  ggplot(test_last, aes(x=type,y=PageRank,colour=type)) +  
  geom_violin(size=1.5) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=13)) +
  theme(legend.position = "None") +
  scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999"))
p

#p2 +scale_fill_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999"))
p + scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999"))


b <-  ggplot(test_last, aes(x=type,y=Betweenness,colour=type)) +
  geom_violin(size=1.5) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=13)) +
  theme(legend.position = "None") +
  scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999"))
b

s <-  ggplot(test_last, aes(x=type,y=Strength)) +
  geom_violin() +
  scale_y_log10() +
  theme_classic()
s

k <-  ggplot(test_last, aes(x=type,y=Katz,colour=type)) +
  geom_violin(size=1.5) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=13)) +
  theme(legend.position = "None") +
  scale_color_manual(values=c("coding" = "#F8766D","lincRNAs" ="#619CFF","miRNAs" ="#F0E442","TFs" ="#999999"))
k

e <-  ggplot(test_last, aes(x=type,y=Eigenvector)) +
  geom_violin() +
  theme_classic()
e

library(cowplot)
all <- plot_grid(p,b,k,labels=c("A","B","C"),vjust = 1.3,ncol=3)#hjust=-20,)
plot(all)
ggsave(filename=here("data/analysis/figures_new//Network_characteristics_by_col_col.png"),device ="png",dpi = 600)
linc_ord <- lincRNAs[order(lincRNAs$PageRank,decreasing = TRUE), ]
newGOA_ord <- as_tibble(newGOA_ord)
linc2col <- linc_ord[,1:2]
linc_10 <- linc2col[1:10,]

# Using tidyverse and dividing by RNA classes
bla <- as_tibble(stats.res_new)
bla2 <- column_to_rownames(bla, var = "gene")

lincRNAs2 <- bla2 %>% filter(grepl("TRINITY",gene)) #%>% 
  #mutate(gene,gene=gsub("TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+_i[0-9]+","linc",lincRNAs$gene))

lincRNAs2 <- column_to_rownames(lincRNAs, var = "gene")
lalala <- lincRNAs2[order(decreasing = TRUE,lincRNAs2$PageRank), ]

#sort(lincRNAs2$PageRank)
coding <- bla %>% filter(grepl("MA_",gene)) #%>% 
  #mutate(gene,gene=gsub("MA_[0-9]+g[0-9]+","gene",genes$gene))
miRNAs <- bla %>% filter(grepl("miRNA",gene)) #%>% 
  #mutate(gene,gene=gsub("miRNA_[0-9]+-[0-9]p","miRNA",miRNAs$gene))



lincRNAs <- lincRNAs %>% add_column(type = "lincRNAs")
miRNAs <- miRNAs %>% add_column(type = "miRNAs")
test <- rbind(lincRNAs,miRNAs)
test2 <- rbind(test,TFs)
test_last <-  rbind(test2,coding)

TF <- read.csv(file = "doc/TF_list.txt",sep = "")
TF_list <- TF$TF_ID

library(VennDiagram)
library(ggplot2)
oveer <- calculate.overlap(list(coding=coding$gene,TF=TF_list,filename=NULL))

xw <-  list(A=NewGOA_genes1$gene,
            B=in_both$gene)
venn_object <- venn.diagram(xw, filename = NULL)
grid.draw(venn_object)
TF_last <- oveer$a2
TF_last <- as_tibble(TF_last)
TF_last <- TF_last %>% add_column(type = "TF")
TF_last <- rename(TF_last, gene = value)
# if it is not working, then colnames(TF_last)[colnames(TF_last) == "value"] <- "gene"
TFs <- semi_join(coding,TF_last,by = NULL, copy=FALSE)
TFs <- TFs %>% add_column(type = "TFs")
coding <- anti_join(coding,TF_last,by = NULL, copy=FALSE)
coding <- coding %>% add_column(type = "coding")


common <- plot.new()
grid.draw(venn.diagram(list(
  genes=genes$gene,
  TF=TF_list),
  filename=NULL))

load("data/analysis/seidr/network.rda")
network <- read_tsv("doc/network_matrix.tsv")
S1_m1 <- median(c(-3,531791,-3.460242,-3.480264))
S2_m1 <- median(c(-3.621773,-3.404639,-3.060226))
S3_m1 <- median(c(3.973096,-3.681682,-3.809513))
S4_m1 <- median(c(4.260013,4.332053,4.551479))
S5_m1 <- median(c(6.407743,5.645703,5.968173))
S6_m1 <- median(c(5.566387,4.651838,3.784862))
S7_m1 <- median(c(3.887946,0,4.794728))
S8_m1 <- median(c(-4.015017,-4.168086,-4.115969))

median(c(1.542634,1.711371,1.725953))
median(c(1.686414,1.25093,1.632859))
median(c(2.124223,2.247199,2.115126))
median(c(1.927139,2.095128,1.785452))
median(c(1.774646,1.800945,1.760106))
median(c(1.467236,1.612047,1.918805))
median(c(2.860756,3.043042,2.784779))
median(c(2.417958,2.365149,2.318025))

cor.test(x=c(-3.230121,-3.404639,-3.681682,4.332053,5.968173,4.651838,3.887946,-4.115969),
         y=c(1.711371,1.632859,2.124223,1.927139,1.774646,1.612047,2.860756,2.365149),
         method = "pearson")


Venn_2vs1_med_3vs2_med <- plot.new()
grid.draw(venn.diagram(list(
  stage2vs1_med=rownames(res_2vs1_med[res_2vs1_med$padj <= 0.01 & abs(res_2vs1_med$log2FoldChange) >= 0.5 & !is.na(res_2vs1_med$padj),]),
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & abs(res_3vs2_med$log2FoldChange) >= 0.5 & !is.na(res_3vs2_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))
# Write the results table 
#write.table(stats.res, file="<myNodeStats.tsv>", row.names = F, sep='\t', quote = F)

# Get first degree neighbours
getGeneFDN <- function(edgeList, gene, source.col=1, target.col=2) {
  # TODO check for data type
  # TODO check that all genes are in the edgelist
  
  s2t <- edgeList[edgeList[source.col] == gene,][,target.col]
  t2s <- edgeList[edgeList[target.col] == gene,][,source.col]
  res <- union(s2t,t2s)
  
  return(res)
}
edgeList <- read.table(here("data/seidr/backbone/edgelist2.txt"))
bla_FDG <- getGeneFDN(edgeList,"TRINITY_DN23326_c0_g1_i1")
bla_FDG <- bla_FDG[! bla_FDG %in% removing]
lincRNA_bla1 <-  bla_FDG[grep("TRINITY",bla_FDG)]
c1 <- plotEigengene(dat, "TRINITY_DN23326_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c1)
ggsave(filename=here("data/analysis/figures_new//candidate1_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes <- unlist(SDN_list) %>% unique()

bla_FDG2 <- getGeneFDN(edgeList,"TRINITY_DN16429_c0_g1_i2")
bla_FDG2 <- bla_FDG2[! bla_FDG2 %in% removing]
lincRNA_bla2 <-  bla_FDG2[grep("TRINITY",bla_FDG2)]
c2 <- plotEigengene(dat, "TRINITY_DN16429_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plot(c2)
ggsave(filename=here("data/analysis/figures_new//candidate2_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG2, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes2 <- unlist(SDN_list) %>% unique()
bla_FDG3 <- getGeneFDN(edgeList,"TRINITY_DN58094_c0_g2_i1")
bla_FDG3 <- bla_FDG3[! bla_FDG3 %in% removing]
lincRNA_bla3 <-  bla_FDG3[grep("TRINITY",bla_FDG3)]
c3 <- plotEigengene(dat, "TRINITY_DN58094_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c3)
ggsave(filename=here("data/analysis/figures_new//candidate3_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG3, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes3 <- unlist(SDN_list) %>% unique()

bla_FDG4 <- getGeneFDN(edgeList,"TRINITY_DN9869_c0_g1_i5")
c4 <- plotEigengene(dat, "TRINITY_DN9869_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plot(c4)
ggsave(filename=here("data/analysis/figures_new//candidate4_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG4, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes4 <- unlist(SDN_list) %>% unique()
bla_FDG5 <- getGeneFDN(edgeList,"TRINITY_DN19049_c0_g1_i1")
bla_FDG5 <- bla_FDG5[! bla_FDG5 %in% removing]
lincRNA_bla5 <-  bla_FDG5[grep("TRINITY",bla_FDG5)]
c5 <- plotEigengene(dat, "TRINITY_DN19049_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c5)
ggsave(filename=here("data/analysis/figures_new//candidate5_PR.png"),device ="png",dpi = 600)
write_excel_csv2(jk,file = "doc/FDNs/genes11.csv",col_names = F,quote = "none")
SDN_list <- lapply(bla_FDG5, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes5 <- unlist(SDN_list) %>% unique()
bla_FDG6 <- getGeneFDN(edgeList,"TRINITY_DN52747_c1_g1_i1")
bla_FDG6 <- bla_FDG6[! bla_FDG6 %in% removing]
lincRNA_bla6 <-  bla_FDG6[grep("TRINITY",bla_FDG6)]
c6 <- plotEigengene(dat, "TRINITY_DN52747_c1_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c6)
ggsave(filename=here("data/analysis/figures_new//candidate6_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG6, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes6 <- unlist(SDN_list) %>% unique()
bla_FDG7 <- getGeneFDN(edgeList,"TRINITY_DN24450_c0_g1_i1")
bla_FDG7 <- bla_FDG7[! bla_FDG7 %in% removing]
lincRNA_bla7 <-  bla_FDG7[grep("TRINITY",bla_FDG7)]
c7 <- plotEigengene(dat, "TRINITY_DN24450_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c7)
ggsave(filename=here("data/analysis/figures_new//candidate7_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG7, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes7 <- unlist(SDN_list) %>% unique()
bla_FDG8 <- getGeneFDN(edgeList,"TRINITY_DN58806_c0_g1_i1")
bla_FDG8 <- bla_FDG8[! bla_FDG8 %in% removing]
lincRNA_bla8 <-  bla_FDG8[grep("TRINITY",bla_FDG8)]
c8 <- plotEigengene(dat, "TRINITY_DN58806_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c8)
ggsave(filename=here("data/analysis/figures_new//candidate8_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG8, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes8 <- unlist(SDN_list) %>% unique()
bla_FDG9 <- getGeneFDN(edgeList,"TRINITY_DN135093_c0_g1_i1")
bla_FDG9 <- bla_FDG9[! bla_FDG9 %in% removing]
lincRNA_bla9 <-  bla_FDG9[grep("TRINITY",bla_FDG9)]
c9 <- plotEigengene(dat, "TRINITY_DN135093_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c9)
ggsave(filename=here("data/analysis/figures_new//candidate9_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG9, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes9 <- unlist(SDN_list) %>% unique()
bla_FDG10 <- getGeneFDN(edgeList,"TRINITY_DN18510_c0_g1_i13")
c10 <- plotEigengene(dat, "TRINITY_DN18510_c0_g1_i13",rep("bla", nrow(dat)),samples$Stages)
plot(c10)
ggsave(filename=here("data/analysis/figures_new//candidate10_PR.png"),device ="png",dpi = 600)

SDN_list <- lapply(bla_FDG10, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes10 <- unlist(SDN_list) %>% unique()
bla_FDG11 <- getGeneFDN(edgeList,"TRINITY_DN15918_c0_g1_i1")
bla_FDG11 <- bla_FDG11[! bla_FDG11 %in% removing]
lincRNA_bla11 <-  bla_FDG11[grep("TRINITY",bla_FDG11)]
c11 <- plotEigengene(dat, "TRINITY_DN15918_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c11)
ggsave(filename=here("data/analysis/figures_new//candidate9_PR.png"),device ="png",dpi = 600)
SDN_list <- lapply(bla_FDG11, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes11 <- unlist(SDN_list) %>% unique()
bla_FDG12 <- getGeneFDN(edgeList,"TRINITY_DN3664_c0_g2_i18")
bla_FDG12 <- bla_FDG12[! bla_FDG12 %in% removing]
lincRNA_bla12 <-  bla_FDG12[grep("TRINITY",bla_FDG12)]
c12 <- plotEigengene(dat, "TRINITY_DN3664_c0_g2_i18",rep("bla", nrow(dat)),samples$Stages)
plot(c12)
ggsave(filename=here("data/analysis/figures_new//candidate10_PR.png"),device ="png",dpi = 600)
SDN_list <- lapply(bla_FDG12, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes12 <- unlist(SDN_list) %>% unique()

#Target genes
plotEigengene(dat, "TRINITY_DN23326_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages,title = "TRINITY_DN23326_c0_g1_i1")
plotEigengene(dat, "MA_10223157g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_10223157g0010")
plotEigengene(dat, "MA_55488g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_55488g0010")
plotEigengene(dat, "MA_7940170g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_7940170g0010")

plotEigengene(dat, "TRINITY_DN16429_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN16429_c0_g1_i2")
plotEigengene(dat, "MA_10257300g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_10257300g0010")

plotEigengene(dat, "TRINITY_DN58094_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN58094_c0_g2_i1")
plotEigengene(dat, "MA_10276052g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_10276052g0010")
plotEigengene(dat, "MA_11119g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_11119g0010")
plotEigengene(dat, "MA_2868651g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_2868651g0010")

plotEigengene(dat, "TRINITY_DN9869_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN9869_c0_g1_i5")
plotEigengene(dat, "MA_19047g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_19047g0010")

plotEigengene(dat, "TRINITY_DN19049_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN19049_c0_g1_i1")
plotEigengene(dat, "MA_32409g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_32409g0010")
plotEigengene(dat, "MA_4763496g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_4763496g0010")

plotEigengene(dat, "TRINITY_DN52747_c1_g1_i1",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN52747_c1_g1_i1")
plotEigengene(dat, "MA_118973g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_118973g0010")
plotEigengene(dat, "MA_8167263g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_8167263g0010")
plotEigengene(dat, "MA_338744g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_338744g0010")

#7 no targets
plotEigengene(dat, "TRINITY_DN24450_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN24450_c0_g1_i1")

plotEigengene(dat, "TRINITY_DN58806_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN58806_c0_g1_i1")
plotEigengene(dat, "MA_2868651g0010",rep("bla", nrow(dat)),samples$Stages,title = "MA_2868651g0010")
plotEigengene(dat, "MA_361914g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_361914g0010")
plotEigengene(dat, "MA_6122032g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_6122032g0010")
plotEigengene(dat, "MA_136565g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_136565g0010")

plotEigengene(dat, "TRINITY_DN135093_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN135093_c0_g1_i1")
plotEigengene(dat, "MA_45153g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_45153g0010")
plotEigengene(dat, "MA_1119638g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_1119638g0010")
plotEigengene(dat, "MA_13740g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_13740g0010")
plotEigengene(dat, "MA_8894812g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_8894812g0010")

plotEigengene(dat, "TRINITY_DN18510_c0_g1_i13",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN18510_c0_g1_i13")
plotEigengene(dat, "MA_4541173g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_4541173g0010")

#11 no targets

plotEigengene(dat, "TRINITY_DN37788_c0_g2_i2",rep("bla", nrow(dat)),samples$Stages, title = "TRINITY_DN37788_c0_g2_i2")
plotEigengene(dat, "MA_118973g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_118973g0010")
plotEigengene(dat, "MA_16299g0010",rep("bla", nrow(dat)),samples$Stages, title = "MA_16299g0010")



# enrichment
bla_FDG_enr9 <- gopher(genes=bla_FDG9, alpha = 0.05, task=list("go", "mapman","kegg","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
bla_FDG_enr_2nd <- gopher(genes=SDN_genes12, alpha = 0.05, task=list("go", "mapman","kegg","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
plotEnrichedTreemap(x = bla_FDG_enr9, enrichment = "mapman") #,namespace = "BP")


load(here("data/analysis/seidr/network.rda"))
samples <- read_csv(here("doc/samples_B2.csv"))
gigetto <- read_tsv("data/analysis/DE/vst-aware_genes+TEs.tsv")
gigi <- read_tsv(here("doc/network_matrix.tsv")) 
plotEigengene(dat, "TRINITY_DN10144_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN101668_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN10359_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN1136_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN11557_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN11723_c0_g1_i14",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN11855_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN11942_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN12056_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN12264_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN1230_c0_g1_i21",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN12781_c0_g5_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN13496_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN13815_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN13931_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN14121_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN1443_c0_g2_i12",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN14642_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN14682_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN1474_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN14880_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15083_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN15131_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15417_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN15464_c0_g3_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN15464_c0_g3_i4",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN15486_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN1559_c0_g2_i5",rep("bla", nrow(dat)),samples$Stages)
###plotEigengene(dat, "TRINITY_DN15730_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN15792_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN15949_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16146_c1_g2_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16253_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16354_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16354_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16363_c5_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16497_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN164_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16538_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16538_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16731_c0_g1_i15",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN16933_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN17064_c0_g3_i1",rep("bla", nrow(dat)),samples$Stages)
#####plotEigengene(dat, "TRINITY_DN17065_c0_g2_i7",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN17186_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN17204_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17439_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN17445_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN17784_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN18097_c0_g1_i18",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN18384_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN18828_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN19087_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN19293_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN1948_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN19939_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN19939_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN20003_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2005_c1_g1_i4",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2080_c0_g2_i5",rep("bla", nrow(dat)),samples$Stages)
#####plotEigengene(dat, "TRINITY_DN20999_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN21012_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN21012_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN21200_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN21484_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN21509_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN21597_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN22622_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22742_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN22795_c0_g2_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22829_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22966_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN23658_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN24162_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2483_c1_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN25022_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2554_c0_g1_i14",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN25721_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN2613_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN26435_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN27258_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2752_c0_g1_i10",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2763_c0_g1_i16",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN276_c1_g2_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN28108_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2980_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN2980_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN29951_c2_g1_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN30140_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN3060_c1_g1_i3",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN30747_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN30822_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN3092_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN30963_c0_g2_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3110_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN31319_c1_g1_i7",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN32172_c0_g1_i12",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN32364_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN3337_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN34138_c0_g2_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN34386_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN34539_c0_g2_i5",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN3517_c0_g1_i12",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN3800_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN38395_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN39226_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN3985_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN4063_c1_g1_i8",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4106_c0_g1_i19",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4216_c0_g1_i10",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4293_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4297_c0_g1_i13",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4327_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN4543_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4994_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN5287_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN5317_c1_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN5339_c3_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN56194_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN56988_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN5831_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN6068_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN6070_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN6119_c0_g1_i12",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN6342_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN82_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN8704_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN8845_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN9613_c0_g1_i14",rep("bla", nrow(dat)),samples$Stages)
##plotEigengene(dat, "TRINITY_DN9996_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)





plotEigengene(dat, bla_FDG12, time = samples$Stages,rep("bla", nrow(dat)))#, multiline = F)


test12 <- as.tibble(dat)
test12<- test12 %>% mutate(samples=rownames(dat))%>% select(samples,bla_FDG12)
test12<-test12 %>% inner_join(samples, by=c("samples"="ID")) %>% select(Stages, bla_FDG12)
test.df12 <- t(as.data.frame(scale(test12%>% select(bla_FDG12))))
colnames(test.df12) = test12$Stages

test <- as.tibble(dat)
test<- test %>% mutate(samples=rownames(dat))%>% select(samples,bla_FDG)
test<-test %>% inner_join(samples, by=c("samples"="ID")) %>% select(Stages, bla_FDG)
test.df <- t(as.data.frame(scale(test%>% select(bla_FDG))))
colnames(test.df) = test$Stages

hm <- heatmap.2(test.df12,
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = colnames(test.df12),
                col=hpal)
#' guy1 trying LncTAr
blabla2 <- as_tibble(bla_FDG)
xx2 <- blabla2 %>% filter(grepl("MA_*", value))
zz2 <- blabla2 %>% filter(grepl("TRINITY", value))

#adding newGOA
newGOA_12 <- read_csv("functional_prediction/results/SDN12_NewGOA.txt",col_names = TRUE)
newGOA_ord <- newGOA_12[order(newGOA$PredictionScore,decreasing = TRUE), ]
newGOA_ord <- as_tibble(newGOA_ord)
newGOA2col <- newGOA_ord[,1:2]
newGOA_30 <- newGOA2col[1:30,]

#comparing newGOA and gopher
load("functional_prediction/semantic_similarity.RData")
GOINFO <- computeIC()
newGOA12c <- c("GO:0006810","GO:0008150","GO:0031410","GO:0015095","GO:1902749","GO:0009132","GO:0043934",
               "GO:0009750","GO:0048830","GO:0071704","GO:0051938","GO:0016049","GO:0003951","GO:1901570",
               "GO:0008484","GO:0048766","GO:0003871","GO:0007129","GO:0019751","GO:0006617","GO:0050789",
               "GO:0046653","GO:0045490","GO:0080058","GO:0070825","GO:0015101","GO:0008865","GO:0045935",
               "GO:0010558","GO:0019405")
gopher11c <- "GO:0045735"
gopher2nd_12 <- c("GO:0032993", "GO:0044815")
cluster3 <- c ("GO:0005506", "GO:0046906", "GO:0048037", "GO:0016020", "GO:0045229", "GO:0016491", "GO:0055114", "GO:0009579",
               "GO:0003824", "GO:0004497", "GO:0009812", "GO:0016705", "GO:0009225", "GO:0008233", "GO:0016825", "GO:0016759",
               "GO:0009536", "GO:0019748", "GO:0016757", "GO:0071554", "GO:1901564", "GO:0010876", "GO:0070001", "GO:0031984",
               "GO:0043169", "GO:0006508", "GO:0009657", "GO:0006063", "GO:0015979", "GO:0043167", "GO:0031224", "GO:0051273",
               "GO:0016758", "GO:0043227", "GO:0006011", "GO:0042214", "GO:0030246", "GO:0004175", "GO:1901361", "GO:0005215",
               "GO:0005976", "GO:0006720", "GO:0016301", "GO:0045552", "GO:0010287", "GO:0006793", "GO:0005975", "GO:0018298",
               "GO:0019538", "GO:0042349", "GO:0004089", "GO:0005737", "GO:0071944", "GO:0042440", "GO:0006979", "GO:0033729",
               "GO:0008171")

t<-mSemSim(newGOA12c,cluster2,GOINFO)
t
BMA(t)


library(ggplot2)
t.melt <- reshape2::melt(t)
p<-ggplot(t.melt, aes(Var1, Var2, fill= value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  ylab("gopher") +
  xlab("newGOA")+
  theme_classic()
p
pdf("functional_prediction/results/guy12_cluster2.pdf", width = 30, height = 30)
p 
graphics.off()
termVector <- c("GO:00001662","GO:0003674","GO:0003824","GO:0005488","GO:0005575","GO:0006807",
                "GO:0008150","GO:0008152","GO:0009987","GO:0016020","GO:0036094","GO:0043167",
                "GO:0043170","GO:0043226","GO:0043227","GO:0043229","GO:0043231","GO:0044237",
                "GO:0044238","GO:0044249","GO:0044260","GO:0050789","GO:0050794","GO:0065007",
                "GO:0071704","GO:0097159","GO:0110165","GO:1901363","GO:1901564","GO:1901576",
                "GO:0009317", "GO:0045735", "GO:0016717", "GO:0006629", "GO:0016885",
                "GO:0004645", "GO:0044281", "GO:0004075", "GO:0009653")
remove_ancestors(termVector,GOINFO)

NewGOA <- read_delim("functional_prediction/results/predictedAnnotation.txt",delim = ",")
range(NewGOA$PredictionScore)
quantile(NewGOA$PredictionScore,probs=seq(0,1,0.01))
source(here("UPSCb-common/src/R/percentile.R"))
percentile(NewGOA$PredictionScore)
plot(density(-log2(NewGOA$PredictionScore)))
abline(v=quantile(-log2(NewGOA$PredictionScore),probs=seq(0,.1,0.01)),lty=2)
#NewGOA <- as_tibble(NewGOA)
#bla_guys <- NewGOA$PredictionScore >= quantile(-log2(NewGOA$PredictionScore),probs=0.02) 
#stupido <- NewGOA %>% filter (quantile(-log2(PredictionScore),probs=0.02))
#scemo <- NewGOA %>% filter (PredictionScore == 1)
#bla_guys <- subset(NewGOA,PredictionScore >= quantile(-log2(NewGOA$PredictionScore),probs=0.02 & PredictionScore <= 1))
NewGOA_linc <- read_delim("functional_prediction/results/lincAnnotation.txt", delim = ",")

NewGOA_linc <- NewGOA_linc[! NewGOA_linc$gene %in% removing,]



range(NewGOA_linc$PredictionScore)
percentile(NewGOA_linc$PredictionScore)
plot(density(-log2(NewGOA_linc$PredictionScore)))
abline(v=quantile(-log2(NewGOA_linc$PredictionScore),probs=seq(0,.1,0.01)),lty=2)
NewGOA_genes <- read_delim("functional_prediction/results/geneAnnotation_new.txt", delim = ",",
                           col_names = c("gene","GOterm","OriginalAnnotation","PredictionScore"))
range(NewGOA_genes$PredictionScore)
percentile(NewGOA_genes$PredictionScore)
plot(density(-log2(NewGOA_genes$PredictionScore)))
abline(v=quantile(-log2(NewGOA_genes$PredictionScore),probs=seq(0,.1,0.01)),lty=2)
par(mgp=c(2.3,1,0))
plot(density(-log2(NewGOA_genes$PredictionScore)), main= NA,ylim=c(0,0.07),col="#F8766D",
     xlab="-log2(PredictionScore)",ylab="density",cex.lab=1.0,lwd=3,bty="l")
title(main="NewGOA selection",adj=0,cex.main= 1.5)
New_GOA_sel = recordPlot()
plot.new()
New_GOA_sel
lines(density(-log2(NewGOA_linc$PredictionScore)),col="#00BFC4",lwd=3,bty="l")
New_GOA2 = New_GOA_sel = recordPlot()
plot.new()
New_GOA2
cor.test(percentile(NewGOA_genes$PredictionScore),percentile(NewGOA_linc$PredictionScore))
abline(v=quantile(-log2(NewGOA_genes$PredictionScore),probs=seq(0,.1,0.03)),lty=2,col="#F8766D",lwd=1.0)
abline(v=quantile(-log2(NewGOA_linc$PredictionScore),probs=seq(0,.1,0.03)),lty=2, col="#00BFC4",lwd=1.0)
legend("right", bty = "n",
       legend=c("coding","lincRNA"), border = c("#F8766D","#00BFC4"),
       col= c("#F8766D","#00BFC4"),pch=22,pt.lwd=1.5,
       title = "type",title.adj= 0.15,cex = 1.0)
NewGOA_last.pdf = New_GOA2 = New_GOA_sel = recordPlot()
plot.new()
NewGOA_last.pdf
dev.off()
pdf("data/analysis/figures_new/NewGOA_last2.pdf")
png("data/analysis/figures_new/NewGOA_last.png",width = 600,height = 350)
NewGOA_genes <- NewGOA_genes %>% add_column(type = "coding")
NewGOA_linc <- NewGOA_linc %>% add_column(type = "lincRNAs")
newGOA_total <- rbind(NewGOA_genes,NewGOA_linc)



library(scales)
show_col(hue_pal()(2))


library(ggplot2)
library(ggridges)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
micio <- ggplot(newGOA_total,aes(x=-log2(PredictionScore),colour=type)) +
  geom_density(size=1.5) +
  labs(colour = NULL) +
  labs(title = "NewGOA selection") +
  theme(text=element_text(size=15)) +
  theme(plot.title = element_text(face = "bold",size=20)) +
  theme_classic()  
  #geom_vline(xintercept=-log2(newGOA_total$PredictionScore),seq(0,.1,by=0.03),col=type)
plot(micio)
ggsave(filename=here("data/analysis/figures_new//NewGOA_selection.png"),device ="png",dpi = 600)
dev.off()

mi2 <- micio + geom_vline(xintercept=-log2(NewGOA_linc$PredictionScore),seq(0,.1,by=0.03)) #,colour=newGOA_total$type)
plot(mi2)
plot(rank(log2(NewGOA_linc$PredictionScore)), 
     rank(log2(NewGOA_genes$PredictionScore)), 
     cex = 0.1, xlab = "linc", ylab = "genes")

#So selecting everything between the value of quantile(-log2(bla$PredictionScore),probs=0.02)  and 1
ex <- NewGOA_linc %>% filter(PredictionScore >= quantile(-log2(PredictionScore),probs=0.01))
ex2 <- subset(NewGOA_linc,-log2(PredictionScore) >= quantile(-log2(PredictionScore),probs=0.01))
length(unique(ex2$gene))
quantile(-log2(NewGOA_linc$PredictionScore),probs=0.01)
2^-13
selected_linc <- subset(NewGOA_linc,PredictionScore >= 2^-13)
selected_genes <- subset(NewGOA_genes,PredictionScore >= 2^-13)
selected_miRNAs <- subset(NewGOA_miRNAs,PredictionScore >= 2^-13)
length(unique(selected_miRNAs$gene))
only_selected_miRNAs <- unique(selected_miRNAs$gene)
selected_linc_ord <- selected_linc[order(selected_linc$PredictionScore,decreasing = TRUE), ]
NewGOA_cluster3_miRNAs <-  subset(only_selected_miRNAs, only_selected_miRNAs %in% miRNA_3$gene)

cl1 <- matrix(c(1037,796,1019,526,23,13),ncol=2,byrow = TRUE)
colnames(cl1) <- c("numbers","numbers_NewGOA")
rownames(cl1) <- c("genes","lincRNAs","miRNAs")
cl1 <- as_tibble(cl1)
#checking how many New_GOA annotation there are in each clusters
NewGOA_miRNAs <- read_delim("functional_prediction/results/miRNA_Annotation.txt", delim = ",")
