#' ---
#' title: "Mfuzz clustering"
#' author: "Camilla Canovi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(matrixStats)
  library(Mfuzz)
  library(readr)
})

suppressMessages(source(here("UPSCb-common/src/R/featureSelection.R")))
suppressMessages(source(here("UPSCb-common/src/R/gopher.R")))

#' * data
load(here("data/analysis/DE/vst-aware_linc.rda"))

#' sample info
samples <- read_csv(here("doc/samples_B2.csv"))

#' # Fuzzy clustering
#' ##Prep
#' Create the eset
conds <- apply(samples[match(colnames(vsta),samples$ID),
                       "Stages"],
               1,paste,collapse="-")
eset <- ExpressionSet(sapply(split.data.frame(t(vsta),conds),colMeans))

#' Remove genes with too little a variation
#' 
#' First we look at the SD distribution
plot(density(rowSds(vsta)))
plot(density(rowSds(vsta)),xlim=c(-0.2,0.5))

#' A cutoff at 0.2 seems adequate
eset <- filter.std(eset,min.std=0.2)

#' Standardise the values
eset <- standardise(eset) 

#' ## Clustering
#' parameter estimation
m1 <- mestimate(eset)

#' check how many clusters we want
Dmin(eset, m=m1, crange=seq(2,25,1), repeats=3, visu=TRUE)

#' cluster
cl <- mfuzz(eset,m1,c=16)

#' ##Plot
colnames(eset)
dir.create(here("data/analysis/Mfuzz"),showWarnings = FALSE)
pdf(file=here("data/analysis/Mfuzz/clusters16.pdf"),width = 16,height=24)
mfuzz.plot2(eset,cl,x11 = FALSE,mfrow=c(3,4),time.labels = colnames(eset),
            centre = TRUE,las=2)
dev.off()

#' * data
load(here("data/analysis/DE/vst-aware_genes+TEs.rda"))
TEs <- vsta[grepl("^MA_", rownames(vsta)) == FALSE, ]
genes <- vsta[grepl("^MA_", rownames(vsta)) == TRUE, ]
#' sample info
samples <- read_csv(here("doc/samples_B2.csv"))

#' # Fuzzy clustering
#' ##Prep
#' Create the eset
conds_genes <- apply(samples[match(colnames(genes),samples$ID),
                       "Stages"],
               1,paste,collapse="-")
eset_genes <- ExpressionSet(sapply(split.data.frame(t(genes),conds_genes),colMeans))

#' Remove genes with too little a variation
#' 
#' First we look at the SD distribution
plot(density(rowSds(genes)))
plot(density(rowSds(genes)),xlim=c(-0.2,1.5))

#' A cutoff at 0.2 seems adequate
eset_genes <- filter.std(eset_genes,min.std=0.3)

#' Standardise the values
eset_genes <- standardise(eset_genes) 

#' ## Clustering
#' parameter estimation
m1 <- mestimate(eset_genes)

#' check how many clusters we want
Dmin(eset_genes, m=m1, crange=seq(2,25,1), repeats=3, visu=TRUE)

#' cluster
cl <- mfuzz(eset_genes,m1,c=16)

#' ##Plot
colnames(eset_genes)
dir.create(here("data/analysis/Mfuzz"),showWarnings = FALSE)
pdf(file=here("data/analysis/Mfuzz/clusters16_genes.pdf"),width = 16,height=24)
mfuzz.plot2(eset_genes,cl,x11 = FALSE,mfrow=c(3,4),time.labels = colnames(eset),
            centre = TRUE,las=2)
dev.off()


#' # Fuzzy clustering
#' ##Prep
#' Create the eset
conds_TEs <- apply(samples[match(colnames(TEs),samples$ID),
                             "Stages"],
                     1,paste,collapse="-")
eset_TEs <- ExpressionSet(sapply(split.data.frame(t(TEs),conds_TEs),colMeans))

#' Remove genes with too little a variation
#' 
#' First we look at the SD distribution
plot(density(rowSds(TEs)))
plot(density(rowSds(TEs)),xlim=c(-0.2,0.4))

#' A cutoff at 0.2 seems adequate
eset_TEs <- filter.std(eset_TEs,min.std=0.2)

#' Standardise the values
eset_TEs <- standardise(eset_TEs) 

#' ## Clustering
#' parameter estimation
m1 <- mestimate(eset_TEs)

#' check how many clusters we want
Dmin(eset_TEs, m=m1, crange=seq(2,25,1), repeats=3, visu=TRUE)

#' cluster
cl <- mfuzz(eset_TEs,m1,c=2)

#' ##Plot
colnames(eset_TEs)
dir.create(here("data/analysis/Mfuzz"),showWarnings = FALSE)
pdf(file=here("data/analysis/Mfuzz/clusters2_TEs.pdf"),width = 16,height=24)
mfuzz.plot2(eset_genes,cl,x11 = FALSE,mfrow=c(3,4),time.labels = colnames(eset),
            centre = TRUE,las=2)
dev.off()


#' ##Membership
str(cl)
barplot(cl$size)

# genes for cluster 20
names(cl$cluster)[cl$cluster == 20]

# enrichment of all clusters
background <- rownames(vst)[featureSelect(vst,unlist(samples[match(colnames(vst),samples$ID),"CONDITION"]),exp=0.1)]

enr.list <- lapply(1:length(cl$size),function(i,clc){
  lapply(names(clc)[clc==i],gopher,background=background,task="go",url="athaliana")
},cl$cluster)

# TODO
# export enr.list as in DifferentialExpression.R to files

# cluster membership
# find the highest membership score per gene
max.membership <- sapply(1:nrow(cl$membership),function(i,m,p){
  m[i,p[i]]
},cl$membership,cl$cluster)

# plot it
plot(density(max.membership))

# find the genes in every cluster that is above a given membership value (0.5 in the example below)
cluster.genes <- lapply(1:length(cl$size),function(i){
  dat <- cl$membership[cl$cluster == i,i]
  plot(density(dat),main=paste("cluster",i))
  abline(v=0.75,lty=2,lwd=2,col="gray")
  names(cl$cluster)[cl$cluster == i][cl$membership[cl$cluster == i,i] >= 0.5]
})

# assess how many genes are still shared between clusters
sum(duplicated(unlist(cl.specific.enr.list)))

# run the enrichment
cl.specific.enr.list <- lapply(cluster.genes,gopher,background=background,task="go",url="athaliana")

# run all enrichment at once
# cl.specific.enr.list <- lapply(cluster.genes,gopher,background=background,task=c("go","kegg","mapman"),url="athaliana")


dev.null <- lapply(1:length(cl.specific.enr.list),function(i){
  r <- cl.specific.enr.list[[i]]
  write_tsv(r$go,path=file.path(file.path(here("data/analysis/Mfuzz",
                                                   paste0("cluster_",i,"_GO-enrichment.tsv")))))
  write_tsv(r$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/Mfuzz",
                                                                    paste0("cluster_",i,"_GO-enrichment_for-REVIGO.tsv")))))
  # write_tsv(r$kegg,path=file.path(file.path(here("data/analysis/Mfuzz",
  #                                              paste0("cluster_",i,"_KEGG-enrichment.tsv")))))
  # write_tsv(r$kegg[,c("id","padj")],path=file.path(file.path(here("data/analysis/Mfuzz",
  # ....                                                              
  
})

# Overlap analysis - check how much overlap there is between clusters and how similar clusters are
# useful to refine the number of clusters
O <- overlap(cl)
Ptmp <- overlap.plot(cl,over=O,thres=0.05)

overlap.plot(cl,over=O,thres=0.1)

# More IDEAS
# get enzyme code enrichment: lapply(cluster.genes,gopher,background=background,task="kegg",url="athaliana")
# retrieve pathway from KEGG (manually) or using KEGGREST (ask Nico)
# Check out the pathview R package
# 
# get mapman enrichment: lapply(cluster.genes,gopher,background=background,task="mapman",url="athaliana") 
# install mapman to create the map 

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```
