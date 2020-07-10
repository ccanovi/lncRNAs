#' ---
#' title: "non-coding RNAs B2 Biological QA"
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
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(pander)
  library(plotly)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Metadata
#' Sample information

samples <- read_csv(here("doc/samples_final.csv"))

#' # Raw data

filelist <- list.files(here("data/Salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Sanity check to ensure that the data is sorted according to the sample info

filelist <- filelist[match(samples$ScilifeID,sub("_sortmerna.*","",basename(dirname(filelist))))]

stopifnot(all(match(sub("_sortmerna.*","",basename(dirname(filelist))),
                    samples$ScilifeID) == 1:nrow(samples)))

#' name the file list vector
names(filelist) <- samples$ID

#' Read the expression at the gene level

counts <- suppressMessages(round(tximport(files = filelist, type = "salmon",txOut=TRUE)$counts))

#' combine technical replicates

samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)
counts <- do.call(
  cbind,
  lapply(split.data.frame(t(counts),
                          samples$ID),
         colSums))

#' write the combined samples 
csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts),csamples$ID),]
write_csv(csamples,path=here("doc/samples_B2.csv"))

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0

sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by Stages

dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(csamples)

ggplot(dat,aes(x,y,fill=csamples$Stages)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected, considering we are using a pool of non-coding RNAs
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by Stages. 

dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Stages=csamples$Stages[match(ind,csamples$ID)])

ggplot(dat,aes(x=values,group=ind,col=Stages)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data_B2.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 

csamples$Stages <- factor(csamples$Stages)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = csamples,
  design = ~Stages)

save(dds,file=here("data/analysis/salmon/dds_B2.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)
save(vst,file=here("data/analysis/DE/vst-blind_B2.rda"))

vsda <- varianceStabilizingTransformation(dds, blind=FALSE)
vsta <- assay(vsda)
vsta <- vsta - min(vsta)
save(vsta,file=here("data/analysis/DE/vst-aware_B2.rda"))

#' * Validation
#' 
#' The variance stabilisation could be worse, considering is a pool of non-coding RNAs
 
meanSdPlot(vsta[rowSums(vsta)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vsta))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=1

#' And the number of possible combinations

nlevel=nlevels(dds$Stages)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    PC4=pc$x[,4],
                    csamples)

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=dds$Stages,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

p <- ggplot(pc.dat,aes(x=PC2,y=PC3,col=dds$Stages,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC2 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[2],"%)",sep="")))

p <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=dds$Stages,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[2],"%)",sep="")))


p <- ggplot(pc.dat,aes(x=PC1,y=PC4,col=dds$Stages,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC4 (",percent[2],"%)",sep="")))

p <- ggplot(pc.dat,aes(x=PC2,y=PC4,col=dds$Stages,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC2 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC4 (",percent[2],"%)",sep="")))


#' ### Heatmap
#' 
#' Filter for noise
conds <- factor(csamples$Stages)
sels <- rangeFeatureSelect(counts=vsta,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 1

#' * Heatmap of "all" genes

par(mar=c(0.05,0.05,0.05,0.05)) 
mar <- par("mar")
hm <- heatmap.2(t(scale(t(vsta[sels[[vst.cutoff+1]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="",
     labels=conds)

#' ## Conclusion
# The Biological QA is fine, considering we are looking at non-coding RNAs. 
# The sequencing depth is good. Also looking at the PCAs we don't have any outliers. 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
