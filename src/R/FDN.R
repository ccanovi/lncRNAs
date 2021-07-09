
library(here)
library(tidyverse)

source(here("UPSCb-common/src/R/gopher.R"))

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
bla <- as_tibble(stats.res)
bla2 <- column_to_rownames(bla, var = "gene")

lincRNAs <- bla %>% filter(grepl("TRINITY",gene))
lincRNAs2 <- column_to_rownames(lincRNAs, var = "gene")
lalala <- lincRNAs2[order(decreasing = TRUE,lincRNAs2$PageRank), ]
#sort(lincRNAs2$PageRank)
genes <- bla %>% filter(grepl("MA_",gene))
miRNAs <- bla %>% filter(grepl("miRNA",gene))

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
SDN_list <- lapply(bla_FDG, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes <- unlist(SDN_list) %>% unique()

bla_FDG2 <- getGeneFDN(edgeList,"TRINITY_DN16429_c0_g1_i2")
SDN_list <- lapply(bla_FDG2, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes2 <- unlist(SDN_list) %>% unique()
bla_FDG3 <- getGeneFDN(edgeList,"TRINITY_DN58094_c0_g2_i1")
SDN_list <- lapply(bla_FDG3, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes3 <- unlist(SDN_list) %>% unique()

bla_FDG4 <- getGeneFDN(edgeList,"TRINITY_DN9869_c0_g1_i5")
SDN_list <- lapply(bla_FDG4, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes4 <- unlist(SDN_list) %>% unique()
bla_FDG5 <- getGeneFDN(edgeList,"TRINITY_DN19049_c0_g1_i1")
SDN_list <- lapply(bla_FDG5, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes5 <- unlist(SDN_list) %>% unique()
bla_FDG6 <- getGeneFDN(edgeList,"TRINITY_DN52747_c1_g1_i1")
SDN_list <- lapply(bla_FDG6, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes6 <- unlist(SDN_list) %>% unique()
bla_FDG7 <- getGeneFDN(edgeList,"TRINITY_DN24450_c0_g1_i1")
SDN_list <- lapply(bla_FDG7, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes7 <- unlist(SDN_list) %>% unique()
bla_FDG8 <- getGeneFDN(edgeList,"TRINITY_DN58806_c0_g1_i1")
SDN_list <- lapply(bla_FDG8, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes8 <- unlist(SDN_list) %>% unique()
bla_FDG9 <- getGeneFDN(edgeList,"TRINITY_DN135093_c0_g1_i1")
SDN_list <- lapply(bla_FDG9, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes9 <- unlist(SDN_list) %>% unique()
bla_FDG10 <- getGeneFDN(edgeList,"TRINITY_DN18510_c0_g1_i13")
SDN_list <- lapply(bla_FDG10, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes10 <- unlist(SDN_list) %>% unique()
bla_FDG11 <- getGeneFDN(edgeList,"TRINITY_DN15918_c0_g1_i1")
SDN_list <- lapply(bla_FDG11, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes11 <- unlist(SDN_list) %>% unique()
bla_FDG12 <- getGeneFDN(edgeList,"TRINITY_DN37788_c0_g2_i2")
SDN_list <- lapply(bla_FDG12, function(gene){
  SDN <- getGeneFDN(edgeList, gene)
})
SDN_genes12 <- unlist(SDN_list) %>% unique()

# enrichment
bla_FDG_enr <- gopher(genes=bla_FDG11, alpha = 0.05, task=list("go", "mapman","kegg","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
bla_FDG_enr_2nd <- gopher(genes=SDN_genes12, alpha = 0.05, task=list("go", "mapman","kegg","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
plotEnrichedTreemap(x = bla_FDG_enr, enrichment = "go",namespace = "BP")


plotEigengene(dat, "TRINITY_DN10144_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN101668_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16354_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN10359_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN1136_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN11557_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN11723_c0_g1_i14",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN11855_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN11942_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN12056_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN12264_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN1230_c0_g1_i21",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN12781_c0_g5_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN13496_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN13815_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN13931_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN14121_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN1443_c0_g2_i12",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN14642_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN14682_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN1474_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN14880_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15083_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15131_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15417_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15464_c0_g3_i2",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN15464_c0_g3_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15486_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN1559_c0_g2_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15730_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15792_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN15949_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16146_c1_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16253_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16354_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16354_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16363_c5_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16497_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN164_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16538_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16538_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN16731_c0_g1_i15",rep("bla", nrow(dat)),samples$Stages)
#plotEigengene(dat, "TRINITY_DN16933_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17064_c0_g3_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17065_c0_g2_i7",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17186_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17204_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17439_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17445_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN17784_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN18097_c0_g1_i18",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN18384_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN18828_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN19087_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN19293_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN1948_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN19939_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN19939_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN20003_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2005_c1_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2080_c0_g2_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN20999_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN21012_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN21012_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN21200_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN21484_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN21509_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN21597_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22622_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22742_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22795_c0_g2_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22829_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN22966_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN23658_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN24162_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2483_c1_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN25022_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2554_c0_g1_i14",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN25721_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2613_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN26435_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN27258_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2752_c0_g1_i10",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2763_c0_g1_i16",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN276_c1_g2_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN28108_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2980_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN2980_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN29951_c2_g1_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN30140_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3060_c1_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN30747_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN30822_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3092_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN30963_c0_g2_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3110_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN31319_c1_g1_i7",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN32172_c0_g1_i12",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN32364_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3337_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN34138_c0_g2_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN34386_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN34539_c0_g2_i5",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3517_c0_g1_i12",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3800_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN38395_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN39226_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN3985_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4063_c1_g1_i8",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4106_c0_g1_i19",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4216_c0_g1_i10",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4293_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4297_c0_g1_i13",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4327_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4543_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN4994_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN5287_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN5317_c1_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN5339_c3_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN56194_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN56988_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN5831_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN6068_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN6070_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN6119_c0_g1_i12",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN6342_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN82_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN8704_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN8845_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN9613_c0_g1_i14",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN9996_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)





plotEigengene(dat, "MA_11119g0010",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "MA_2868651g0010",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, bla_FDG3, time = samples$Stages,rep("bla", nrow(dat)))#, multiline = F)


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

blah <- read_delim("functional_prediction/results/predictedAnnotation.txt",delim = ",")
range(blah$PredictionScore)
quantile(blah$PredictionScore,probs=seq(0,1,0.01))
source(here("UPSCb-common/src/R/percentile.R"))
percentile(blah$PredictionScore)
plot(density(-log1p(blah$PredictionScore)))
plot(density(blah$PredictionScore))
abline(v=quantile(-log2(blah$PredictionScore),probs=seq(0,.1,0.01)),lty=2)
blah <- as_tibble(blah)
bla_guys <- blah$PredictionScore >= quantile(-log2(blah$PredictionScore),probs=0.02) 
stupido <- blah %>% filter (quantile(-log2(PredictionScore),probs=0.02))
scemo <- blah %>% filter (PredictionScore == 1)
bla_guys <- subset(blah,PredictionScore >= quantile(-log2(blah$PredictionScore),probs=0.02 & PredictionScore <= 1))
blah_linc <- read_delim("functional_prediction/results/lincAnnotation.txt", delim = ",")
range(blah_linc$PredictionScore)
percentile(blah_linc$PredictionScore)
plot(density(blah_linc$PredictionScore))
plot(density(-log2(blah_linc$PredictionScore)))
abline(v=quantile(-log2(blah_linc$PredictionScore),probs=seq(0,.1,0.01)),lty=2)
blah_genes <- read_delim("functional_prediction/results/geneAnnotation.txt", delim = ",")
range(blah_genes$PredictionScore)
percentile(blah_genes$PredictionScore)
plot(density(blah_genes$PredictionScore))
plot(density(-log2(blah_genes$PredictionScore)))
abline(v=quantile(-log2(blah_genes$PredictionScore),probs=seq(0,.1,0.01)),lty=2)
plot(density(-log2(blah_genes$PredictionScore)))
lines(density(-log2(blah_linc$PredictionScore),col=2))
cor.test(percentile(blah_genes$PredictionScore),percentile(blah_linc$PredictionScore))
abline(v=quantile(-log2(blah_genes$PredictionScore),probs=seq(0,.1,0.03)),lty=2)
abline(v=quantile(-log2(blah_linc$PredictionScore),probs=seq(0,.1,0.03)),lty=2, col=2)
plot(rank(log2(blah_linc$PredictionScore)), 
     rank(log2(blah_genes$PredictionScore)), 
     cex = 0.1, xlab = "linc", ylab = "genes")

#So selecting everything between the value of quantile(-log2(bla$PredictionScore),probs=0.02)  and 1
ex <- blah_linc %>% filter(PredictionScore >= quantile(-log2(PredictionScore),probs=0.01))
ex2 <- subset(blah_linc,-log2(PredictionScore) >= quantile(-log2(PredictionScore),probs=0.01))
length(unique(ex2$gene))
quantile(-log2(blah_linc$PredictionScore),probs=0.01)
2^-13
ex3 <- subset(blah_linc,PredictionScore >= 2^-13)
length(unique(ex3$gene))
ex3_ord <- ex3[order(ex3$PredictionScore,decreasing = TRUE), ]

load(here("data/analysis/DE/vst-aware_B2.rda"))
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
samples_m <- read.csv("doc/samples_final.csv")
time_expression <- expressionSpecificity(exp.mat = vsta[,samples_m$ID],
                                         tissues = as.character(samples_m$Stages),
                                         output = "complete")

stuff <- as_tibble(time_expression)
stuff_new <- stuff %>% add_column(gene = ex3$gene)
final_tibble_nc <- left_join(non_coding, stuff_new, by = NULL, copy=FALSE) %>% 
  rename_if(grepl("aij",colnames(.)),funs(str_replace(.,"aij\\.",""))) %>% 
  add_column(gene = ex3$gene)