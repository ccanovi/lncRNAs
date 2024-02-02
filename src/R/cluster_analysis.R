suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(treemap)
  library(jsonlite)
})
source(here("UPSCb-common/src/R/gopher.R"))

#Correlation between genes 
network <- read.delim("data/seidr/backbone/edgeList.txt", header=TRUE, stringsAsFactors=FALSE)

#Genes in Cluster 
InfomapClusters <- read.delim("data/seidr/backbone/infomapClusters.tsv")
cluster1 <- subset(InfomapClusters, cluster == "Cluster1")
cluster2 <- subset(InfomapClusters, cluster == "Cluster2")
cluster3 <- subset(InfomapClusters, cluster == "Cluster3")
cluster4 <- subset(InfomapClusters, cluster == "Cluster4")
cluster5 <- subset(InfomapClusters, cluster == "Cluster5")
cluster6 <- subset(InfomapClusters, cluster == "Cluster6")
cluster7 <- subset(InfomapClusters, cluster == "Cluster7")
cluster8 <- subset(InfomapClusters, cluster == "Cluster8")
cluster9 <- subset(InfomapClusters, cluster == "Cluster9")
cluster10 <- subset(InfomapClusters, cluster == "Cluster10")
cluster11 <- subset(InfomapClusters, cluster == "Cluster11")
cluster12 <- subset(InfomapClusters, cluster == "Cluster12")
cluster13 <- subset(InfomapClusters, cluster == "Cluster13")
cluster14 <- subset(InfomapClusters, cluster == "Cluster14")
cluster15 <- subset(InfomapClusters, cluster == "Cluster15")

#Gopher
cl1_enr <- gopher(genes=cluster1$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl2_enr <- gopher(genes=cluster2$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl3_enr <- gopher(genes=cluster3$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl4_enr <- gopher(genes=cluster4$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl5_enr <- gopher(genes=cluster5$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl6_enr <- gopher(genes=cluster6$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl7_enr <- gopher(genes=cluster7$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl8_enr <- gopher(genes=cluster8$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl9_enr <- gopher(genes=cluster9$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl10_enr <- gopher(genes=cluster10$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl11_enr <- gopher(genes=cluster11$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl12_enr <- gopher(genes=cluster12$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl13_enr <- gopher(genes=cluster13$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl14_enr <- gopher(genes=cluster14$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")
cl15_enr <- gopher(genes=cluster15$gene, alpha = 0.05, task=list("go", "mapman","pfam"), background = InfomapClusters$gene, url="pabies", endpoint = "enrichment")

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

#run the function plotEigengene.R in src/R
load(here("data/analysis/seidr/network.rda"))
rownames(dat) <- sub("\\.1","",(rownames(dat)))
samples <- read_csv(here("doc/samples_B2.csv"))

#plot for each cluster and each different class of RNA
plotEigengene(t(dat), cluster1$gene,samples$Stages,rep("stages", nrow(t(dat))))
gene14 <- cluster14$gene[grep("MA",cluster14$gene)]
plotEigengene(dat, gene8, samples$Stages, rep("stages", nrow(t(dat))))
miRNA14 <- cluster14$gene[grep("miRNA",cluster14$gene)] 
plotEigengene(dat, miRNA5, samples$Stages, rep("bla", nrow(dat)))
lincRNA14 <- cluster14$gene[grep("TRINITY",cluster14$gene)]

#plotting clusters
plotEigengene(t(dat), cluster1$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster1",colors = c("#1B9E77"))
ggsave(filename=here("data/analysis/figures//cluster1.png"),device ="png",dpi = 600)
pdf("data/analysis/figures//cluster2.pdf")
plotEigengene(t(dat), cluster2$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster2",colors = c("#A6761D"))
dev.off()
ggsave(filename=here("data/analysis/figures//cluster2.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster3$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster3",colors = c("#89CE00"))
ggsave(filename=here("data/analysis/figures//cluster3.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster4$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster4",colors = c("#E7298A"))
ggsave(filename=here("data/analysis/figures//cluster4.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster5$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster5",colors = c("#66A61E"))
ggsave(filename=here("data/analysis/figures//cluster5.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster6$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster6",colors = c("#E6AB02"))
ggsave(filename=here("data/analysis/figures//cluster6.png"),device ="png",dpi = 600)
pdf("data/analysis/figures//cluster7.pdf")
plotEigengene(t(dat), cluster7$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster7",colors = c("#7570B3"))
dev.off()
ggsave(filename=here("data/analysis/figures//cluster7.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster8$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster8",colors = c("#666666"))
ggsave(filename=here("data/analysis/figures//cluster8.png"),device ="png",dpi = 600)
pdf("data/analysis/figures//cluster9.pdf")
plotEigengene(t(dat), cluster9$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster9",colors = c("#56B4E9"))
dev.off()
ggsave(filename=here("data/analysis/figures//cluster9.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster10$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster10",colors = c("#661100"))
ggsave(filename=here("data/analysis/figures//cluster10.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster11$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster11",colors = c("#CC79A7"))
ggsave(filename=here("data/analysis/figures//cluster11.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster12$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster12",colors = c("#DDCC77"))
ggsave(filename=here("data/analysis/figures//cluster12.png"),device ="png",dpi = 600)
pdf("data/analysis/figures//cluster13.pdf")
plotEigengene(t(dat), cluster13$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster13",colors = c("#F8B8D0"))
dev.off()
ggsave(filename=here("data/analysis/figures//cluster13.png"),device ="png",dpi = 600)
plotEigengene(t(dat), cluster14$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster14",colors = c("#E7298A"))
ggsave(filename=here("data/analysis/figures//cluster14.png"),device ="png",dpi = 600)
pdf("data/analysis/figures//cluster15.pdf")
plotEigengene(t(dat), cluster15$gene, rep("stages", nrow(t(dat))),samples$Stages,title = "Cluster15",colors = c("#332288"))
dev.off()
ggsave(filename=here("data/analysis/figures//cluster15.png"),device ="png",dpi = 600)

#counting elements in each cluster
gene_1 <-  cluster1[grep("MA",cluster1$gene),]
miRNA_1 <-  cluster1[grep("miRNA",cluster1$gene),]
lincRNA_1 <-  cluster1[grep("TRINITY",cluster1$gene),]

gene_2 <-  cluster2[grep("MA",cluster2$gene),]
miRNA_2 <-  cluster2[grep("miRNA",cluster2$gene),]
lincRNA_2 <-  cluster2[grep("TRINITY",cluster2$gene),]

gene_3 <-  cluster3[grep("MA",cluster3$gene),]
miRNA_3 <-  cluster3[grep("miRNA",cluster3$gene),]
lincRNA_3 <-  cluster3[grep("TRINITY",cluster3$gene),]

gene_4 <-  cluster4[grep("MA",cluster4$gene),]
miRNA_4 <-  cluster4[grep("miRNA",cluster4$gene),]
lincRNA_4 <-  cluster4[grep("TRINITY",cluster4$gene),]

gene_5 <-  cluster5[grep("MA",cluster5$gene),]
miRNA_5 <-  cluster5[grep("miRNA",cluster5$gene),]
lincRNA_5 <-  cluster5[grep("TRINITY",cluster5$gene),]

gene_6 <-  cluster6[grep("MA",cluster6$gene),]
miRNA_6 <-  cluster6[grep("miRNA",cluster6$gene),]
lincRNA_6 <-  cluster6[grep("TRINITY",cluster6$gene),]

gene_7 <-  cluster7[grep("MA",cluster7$gene),]
miRNA_7 <-  cluster7[grep("miRNA",cluster7$gene),]
lincRNA_7 <-  cluster7[grep("TRINITY",cluster7$gene),]

gene_8 <-  cluster8[grep("MA",cluster8$gene),]
miRNA_8 <-  cluster8[grep("miRNA",cluster8$gene),]
lincRNA_8 <-  cluster8[grep("TRINITY",cluster8$gene),]

gene_9 <-  cluster9[grep("MA",cluster9$gene),]
miRNA_9 <-  cluster9[grep("miRNA",cluster9$gene),]
lincRNA_9 <-  cluster9[grep("TRINITY",cluster9$gene),]

gene_10 <-  cluster10[grep("MA",cluster10$gene),]
miRNA_10 <-  cluster10[grep("miRNA",cluster10$gene),]
lincRNA_10 <-  cluster10[grep("TRINITY",cluster10$gene),]

gene_11 <-  cluster11[grep("MA",cluster11$gene),]
miRNA_11 <-  cluster11[grep("miRNA",cluster11$gene),]
lincRNA_11 <-  cluster11[grep("TRINITY",cluster11$gene),]

gene_12 <-  cluster12[grep("MA",cluster12$gene),]
miRNA_12 <-  cluster12[grep("miRNA",cluster12$gene),]
lincRNA_12 <-  cluster12[grep("TRINITY",cluster12$gene),]

gene_13 <-  cluster13[grep("MA",cluster13$gene),]
miRNA_13 <-  cluster13[grep("miRNA",cluster13$gene),]
lincRNA_13 <-  cluster13[grep("TRINITY",cluster13$gene),]

gene_14 <-  cluster14[grep("MA",cluster14$gene),]
miRNA_14 <-  cluster14[grep("miRNA",cluster14$gene),]
lincRNA_14 <-  cluster14[grep("TRINITY",cluster14$gene),]

gene_15 <-  cluster15[grep("MA",cluster15$gene),]
miRNA_15 <-  cluster15[grep("miRNA",cluster15$gene),]
lincRNA_15 <-  cluster15[grep("TRINITY",cluster15$gene),]

