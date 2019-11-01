#' ---
#' title: "non-coding RNAs Biological QA"
#' author: "Camilla Canovi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(matrixStats)
  library(parallel)
  library(pander)
  library(plotly)
  library(RColorBrewer)
  library(scatterplot3d)
  library(tidyverse)
  library(tximport)
  library(vsn)
})
#' Helper functions
source(here("UPSCb-common/src/R/plot.multidensity.R"))
source(here("UPSCb-common/src/R/featureSelection.R"))

#' Graphics
#pal <- brewer.pal(8,"Dark2")
pal <- brewer.pal(10,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' Metadata
#' Sample information

samples <- read_csv(here("doc/samples.csv"))

samples$ID <- sub(".*/","",sub("_L00[1,2]","",sub("[1,2]_[1].*_.*DXX_","", sub("_dual.*","",samples$ScilifeID))))  

samples$Batch <- factor(sprintf("B%d",grepl("P7614",samples$ScilifeID)+1))
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

# For the PCA, Use the vst data instead and use colMedians or colMeans
 tst <- sapply(split.data.frame(t(counts),
                         samples$Stages),
         colSums)

#' ## Quality Control
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' Let us take a look at the sequencing depth, colouring by Stages

dat_s <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat_s,aes(x,y,fill=samples$Stages)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by Batch. 

dat_b <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Batch=samples$Batch[match(ind,samples$ID)]) %>% 
  mutate(Stages=samples$Stages[match(ind,samples$ID)])

ggplot(dat_b,aes(x=values,group=ind,col=Batch)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

ggplot(dat_b,aes(x=values,group=ind,col=Stages)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 

samples$Stages <- factor(samples$Stages)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Batch+Stages)

save(dds,file=here("data/analysis/salmon/dds.rda"))

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

dir.create(here("data/analysis/DE"))
save(vst,file=here("data/analysis/DE/vst-blind.rda"))

vsda <- varianceStabilizingTransformation(dds, blind=FALSE)
vsta <- assay(vsda)
vsta <- vsta - min(vsta)
save(vsta,file=here("data/analysis/DE/vst-aware.rda"))

#' Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=2

nlevel=nlevels(dds$Stages) * nlevels(dds$Batch)

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
  
#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)

#' The PCA shows that a large fraction of the variance is 
#' explained by both variables.
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(dds$Stages)],
              pch=c(17,19)[as.integer(dds$Batch)])
legend("topleft",
       fill=pal[1:nlevels(dds$Stages)],
       legend=levels(dds$Stages))

legend("topright",
       pch=c(17,19),
       legend=levels(dds$Batch))

par(mar=mar)

#' ### 2D
pc.dat <- cbind(pc$x,samples)

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=dds$Stages,shape=dds$Batch,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

p <- ggplot(pc.dat,aes(x=PC2,y=PC3,col=dds$Stages,shape=dds$Batch,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC2 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[2],"%)",sep="")))

p <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=dds$Stages,shape=dds$Batch,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[2],"%)",sep="")))


p <- ggplot(pc.dat,aes(x=PC1,y=PC4,col=dds$Stages,shape=dds$Batch,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC4 (",percent[2],"%)",sep="")))

p <- ggplot(pc.dat,aes(x=PC2,y=PC4,col=dds$Stages,shape=dds$Batch,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC2 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC4 (",percent[2],"%)",sep="")))



#' ### Heatmap
#' 
#' Filter for noise
#' 
conds <- factor(paste(samples$Stages,samples$Batch))
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 10

#' Heatmap of "all" genes
#' 
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="")


#' ## Conclusion
#' CHANGEME
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
