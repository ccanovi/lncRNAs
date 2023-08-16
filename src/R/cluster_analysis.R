suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(treemap)
  library(jsonlite)
})
source(here("UPSCb-common/src/R/gopher.R"))

#Correlation between genes 
network <- read.delim("data/seidr/backbone/edgelist2.txt", header=TRUE, stringsAsFactors=FALSE)

#Genes in Cluster 
InfomapClusters <- read.delim("data/seidr/backbone/infomapClusters.tsv", stringsAsFactors = FALSE)
InfomapClusters <- InfomapClusters[! InfomapClusters$gene %in% removing,]
#interesting clusters in the nettwork
cluster1 <- subset(InfomapClusters, cluster == "Cluster1")
cluster2 <- subset(InfomapClusters, cluster == "Cluster2")
cluster3 <- subset(InfomapClusters, cluster == "Cluster3")
cluster4 <- subset(InfomapClusters, cluster == "Cluster4")
cluster5 <- subset(InfomapClusters, cluster == "Cluster5")
cluster6 <- subset(InfomapClusters, cluster == "Cluster6")
cluster7 <- subset(InfomapClusters, cluster == "Cluster7")
cluster8 <- subset(InfomapClusters, cluster == "Cluster8")
cluster9 <- subset(InfomapClusters, cluster == "Cluster9")
#Gopher

cl1_enr <- gopher(genes=cluster1$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl2_enr <- gopher(genes=cluster2$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl3_enr <- gopher(genes=cluster3$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl4_enr <- gopher(genes=cluster4$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl5_enr <- gopher(genes=cluster5$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl6_enr <- gopher(genes=cluster6$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl7_enr <- gopher(genes=cluster7$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl8_enr <- gopher(genes=cluster8$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")

#treemap

plotEnrichedTreemap(x = cl1_enr, enrichment = "go") #, namespace = "CC")


clusterTreemapColors <- rep(c("#9E1F63","#662D91","#2B3990","#1B75BC","#27AAE1",
                              "#2AB592","#035E31","#009444","#65BC46","#A5CE42",
                              "#F9ED32","#FBB040","#F15A29","#EF4136","#BE1E2D"),15)
clusterTreemapText <- rep(c("white","white","white","white","white",
                            "white","white","white","white","white",
                            "black","black","white","white","white"),15)
plotEnrichedTreemap <- function(x, enrichment = c('go','mapman', 'kegg', 'pfam', 'ko_pathway', 'ko', 'kog','cog'), 
                                namespace = c('none', 'BP', 'MF', 'CC'), 
                                title = "", 
                                de = c("none", "up", "down"), 
                                clusterColor = "#9E1F63", clusterText='black',
                                nameCol="name",
                                namespaceCol="namespace",
                                sizeCol="padj",
                                colorCol="padj", 
                                convertSize=TRUE,legend=TRUE) {
  require(treemap)
  require(RColorBrewer)
  enrichment <- match.arg(enrichment)
  namespace <- match.arg(namespace)
  de <- match.arg(de)
  enrData <- if (is.data.frame(x)){
    x 
  } else if (is.list(x)) {
    x[[enrichment]]
  } else {
    stop("x object must be either a list or a data.frame")
  }
  # calculate the size based on padj
  enrData$size <- if(convertSize) {
    abs(log10(enrData[[colorCol]]))
  } else {
    enrData[[colorCol]]
  }
  #default treemap
  index = nameCol
  fontcolor.labels=clusterText
  fontface.labels=c(2)
  fontsize.labels=c(12)
  inflate.labels = TRUE
  align.labels=list(c("center", "center"))
  position.legend <- if (!legend) {"none"} else {"bottom"}
  #vColor = 'name'
  border.col='black'
  border.lwds=c(7,2)
  palette <- colorRampPalette(c("white", clusterColor))
  palette <- palette(10)
  vColor = c("size")
  if (sizeCol=="padj") {
    vSize = vColor
  } else {
    vSize = sizeCol
  }
  type = "value"
  title.legend="abs(log10(pAdj))"
  bg.labels= 0
  if(enrichment =='go' ){
    if(namespace=='none') {
      index = c('namespace', 'name')
      palette = "Set1"
      fontcolor.labels=c("#FF000000", "black")
      fontface.labels=c(1, 2)
      fontsize.labels=c(1, 20)
      inflate.labels=TRUE
      align.labels=list(c("left", "top"), c("center", "center") )
      type = "index"
      border.col=c('black', 'white')
      title.legend="GO namespace"
    } else {
      enrData <- enrData[enrData[[namespaceCol]]==namespace,]
    }
  }
  # mapman name fix for better visualization
  if(enrichment =='mapman') {
    enrData[[nameCol]] <- gsub("\\."," ",enrData[[nameCol]])
  } 
  # Paint it properly if it is up or down regulated
  if(de !='none') {
    if(de=='up') {
      palette = "OrRd"
    }
    if(de=='down') {
      palette = "GnBu"
    }
  } 
  # generate treemap
  treemap(enrData, 
          index = index,
          vSize = vSize, 
          palette = palette,
          type = type, 
          vColor = vColor,
          title=title, 
          fontcolor.labels=fontcolor.labels, 
          fontface.labels=fontface.labels,     
          #fontsize.labels=fontsize.labels,
          bg.labels=bg.labels, 
          inflate.labels = inflate.labels ,
          lowerbound.cex.labels = 0, 
          position.legend = position.legend,
          border.col=border.col,         
          border.lwds=border.lwds,
          align.labels=align.labels,
          title.legend=title.legend,
          overlap.labels = 1
  )
}
#try to do Eigen plot

library("devtools")
devtools::install_github("loalon/Rtoolbox")
library(rBiobox)
plotEigengene(toyData$expressionData, toyData$geneCluster, toyData$time, toyData$conditions)

?plotEigengene()
load(here("data/analysis/seidr/network.rda"))
samples <- read_csv(here("doc/samples_B2.csv"))
plotEigengene(dat, cluster7$gene, samples$Stages, rep("bla", nrow(dat)))
plotEigengene(dat, try, samples$Stages, rep("bla", nrow(dat)))

#plot for each cluster and each different class of RNA
plotEigengene(dat, cluster1$gene, samples$Stages, rep("bla", nrow(dat)))
gene8 <- cluster8$gene[grep("MA",cluster8$gene)]
plotEigengene(dat, gene8, samples$Stages, rep("bla", nrow(dat)))
miRNA5 <- cluster5$gene[grep("miRNA",cluster5$gene)] 
plotEigengene(dat, miRNA5, samples$Stages, rep("bla", nrow(dat)))
lincRNA8 <- cluster8$gene[grep("TRINITY",cluster8$gene)]
plotEigengene(dat, lincRNA8, samples$Stages, rep("bla", nrow(dat)))

plotEigengene(dat, cluster1$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster1",colors = c("#1B9E77"))
ggsave(filename=here("data/analysis/figures_new//cluster1.png"),device ="png",dpi = 600)
plotEigengene(dat, cluster2$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster2",colors = c("#D95F02"))
ggsave(filename=here("data/analysis/figures_new//cluster2.png"),device ="png",dpi = 600)
plotEigengene(dat, cluster3$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster3",colors = c("#A6761D"))
ggsave(filename=here("data/analysis/figures_new//cluster3.png"),device ="png",dpi = 600)
plotEigengene(dat, cluster4$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster4",colors = c("#E7298A"))
ggsave(filename=here("data/analysis/figures_new//cluster4.png"),device ="png",dpi = 600)
plotEigengene(dat, cluster5$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster5",colors = c("#66A61E"))
ggsave(filename=here("data/analysis/figures_new//cluster5.png"),device ="png",dpi = 600)
plotEigengene(dat, cluster6$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster6",colors = c("#E6AB02"))
ggsave(filename=here("data/analysis/figures_new//cluster6.png"),device ="png",dpi = 600)
plotEigengene(dat, cluster7$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster7",colors = c("#7570B3"))
ggsave(filename=here("data/analysis/figures_new//cluster7.png"),device ="png",dpi = 600)
plotEigengene(dat, cluster8$gene, rep("bla", nrow(dat)),samples$Stages,title = "Cluster8",colors = c("#666666"))
ggsave(filename=here("data/analysis/figures_new//cluster8.png"),device ="png",dpi = 600)



#doing bar_plots
xar <- read_tsv(here("doc/tissue_specificity_network.tsv"))

gene_1 <-  cluster1[grep("MA",cluster1$gene),]
gene_1 <- gene_1 %>% add_column(type = "gene")
miRNA_1 <-  cluster1[grep("miRNA",cluster1$gene),]
miRNA_1 <- miRNA_1 %>% add_column(type = "miRNA")
lincRNA_1 <-  cluster1[grep("TRINITY",cluster1$gene),]
lincRNA_1 <- lincRNA_1 %>% add_column(type = "linc")
colnames(gene_1)[1] <- "network$ID"
bla <- gene_1 %>% add_row(lincRNA_1)
Cluster1 <- bla %>% add_row(miRNA_1)

gene_2 <-  cluster2[grep("MA",cluster2$gene),]
gene_2 <- gene_2 %>% add_column(type = "gene")
miRNA_2 <-  cluster2[grep("miRNA",cluster2$gene),]
miRNA_2 <- miRNA_2 %>% add_column(type = "miRNA")
lincRNA_2 <-  cluster2[grep("TRINITY",cluster2$gene),]
lincRNA_2 <- lincRNA_2 %>% add_column(type = "linc")
colnames(lincRNA_2)[1] <- "network$ID"
bla2 <- gene_2 %>% add_row(lincRNA_2)
Cluster2 <- bla2 %>% add_row(miRNA_2)

gene_3 <-  cluster3[grep("MA",cluster3$gene),]
gene_3 <- gene_3 %>% add_column(type = "gene")
miRNA_3 <-  cluster3[grep("miRNA",cluster3$gene),]
miRNA_3 <- miRNA_3 %>% add_column(type = "miRNA")
lincRNA_3 <-  cluster3[grep("TRINITY",cluster3$gene),]
lincRNA_3 <- lincRNA_3 %>% add_column(type = "linc")
colnames(lincRNA_3)[1] <- "network$ID"
bla3 <- gene_3 %>% add_row(lincRNA_3)
Cluster3 <- bla3 %>% add_row(miRNA_3)

gene_4 <-  cluster4[grep("MA",cluster4$gene),]
gene_4 <- gene_4 %>% add_column(type = "gene")
miRNA_4 <-  cluster4[grep("miRNA",cluster4$gene),]
miRNA_4 <- miRNA_4 %>% add_column(type = "miRNA")
lincRNA_4 <-  cluster4[grep("TRINITY",cluster4$gene),]
lincRNA_4 <- lincRNA_4 %>% add_column(type = "linc")
colnames(lincRNA_4)[1] <- "network$ID"
bla4 <- gene_4 %>% add_row(lincRNA_4)
Cluster4 <- bla4 %>% add_row(miRNA_4)

gene_5 <-  cluster5[grep("MA",cluster5$gene),]
gene_5 <- gene_5 %>% add_column(type = "gene")
miRNA_5 <-  cluster5[grep("miRNA",cluster5$gene),]
miRNA_5 <- miRNA_5 %>% add_column(type = "miRNA")
lincRNA_5 <-  cluster5[grep("TRINITY",cluster5$gene),]
lincRNA_5 <- lincRNA_5 %>% add_column(type = "linc")
colnames(gene_5)[1] <- "network$ID"
bla5 <- gene_5 %>% add_row(lincRNA_5)
Cluster5 <- bla5 %>% add_row(miRNA_5)

gene_6 <-  cluster6[grep("MA",cluster6$gene),]
gene_6 <- gene_6 %>% add_column(type = "gene")
miRNA_6 <-  cluster6[grep("miRNA",cluster6$gene),]
lincRNA_6 <-  cluster6[grep("TRINITY",cluster6$gene),]
lincRNA_6 <- lincRNA_6 %>% add_column(type = "linc")
colnames(lincRNA_6)[1] <- "network$ID"
Cluster6 <- gene_6 %>% add_row(lincRNA_6)


gene_7 <-  cluster7[grep("MA",cluster7$gene),]
gene_7 <- gene_7 %>% add_column(type = "gene")
miRNA_7 <-  cluster7[grep("miRNA",cluster7$gene),]
lincRNA_7 <-  cluster7[grep("TRINITY",cluster7$gene),]
lincRNA_7 <- lincRNA_7 %>% add_column(type = "linc")
colnames(gene_7)[1] <- "network$ID"
Cluster7 <- gene_7 %>% add_row(lincRNA_7)

gene_8 <-  cluster8[grep("MA",cluster8$gene),]
gene_8 <- gene_8 %>% add_column(type = "gene")
miRNA_8 <-  cluster8[grep("miRNA",cluster8$gene),]
lincRNA_8 <-  cluster8[grep("TRINITY",cluster8$gene),]
lincRNA_8 <- lincRNA_8 %>% add_column(type = "linc")
colnames(lincRNA_8)[1] <- "network$ID"
Cluster8 <- gene_8 %>% add_row(lincRNA_8)

bla_2 <- Cluster1 %>% add_row(Cluster2)
bla_3 <- bla_2 %>% add_row(Cluster3)
bla_4 <- bla_3 %>% add_row(Cluster4)
bla_5 <- bla_4 %>% add_row(Cluster5)
bla_6 <- bla_5 %>% add_row(Cluster6)
bla_7 <- bla_6 %>% add_row(Cluster7)
clusters <- bla_7 %>% add_row(Cluster8)

xar_b2 <- left_join(xar, clusters, by = "network$ID",copy=FALSE)
xar_b2 <- xar_b2  %>% filter(cluster != "NA")
xar_1 <- left_join(xar,Cluster1,by = "network$ID",copy=FALSE)
xar_1 <- xar_1  %>% filter(cluster != "NA")
xar_2 <- left_join(xar,Cluster2,by = "network$ID",copy=FALSE)
xar_2 <- xar_2  %>% filter(cluster != "NA")
xar_3 <- left_join(xar,Cluster3,by = "network$ID",copy=FALSE)
xar_3 <- xar_3  %>% filter(cluster != "NA")
xar_4 <- left_join(xar,Cluster4,by = "network$ID",copy=FALSE)
xar_4 <- xar_4  %>% filter(cluster != "NA")
xar_5 <- left_join(xar,Cluster5,by = "network$ID",copy=FALSE)
xar_5 <- xar_5  %>% filter(cluster != "NA")
xar_6 <- left_join(xar,Cluster6,by = "network$ID",copy=FALSE)
xar_6 <- xar_6  %>% filter(cluster != "NA")
xar_7 <- left_join(xar,Cluster7,by = "network$ID",copy=FALSE)
xar_7 <- xar_7  %>% filter(cluster != "NA")
xar_8 <- left_join(xar,Cluster8,by = "network$ID",copy=FALSE)
xar_8 <- xar_8  %>% filter(cluster != "NA")


library(ggplot2)

f <- ggplot(xar_b2) +
  geom_bar(aes(x=peak,fill=type))
f 
f + facet_wrap(~type)

g <- ggplot(xar_b2) +
  geom_bar(aes(x=peak,fill=cluster))
g
g + facet_wrap(~cluster)

h <- ggplot(xar_3) +
  geom_bar(aes(x=peak,fill=type))
h
h + facet_wrap(~type)
