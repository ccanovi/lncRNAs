#finding candidates based on expression profile

library(here)
library(tidyverse)
source(here("UPSCb-common/src/R/gopher.R"))
load(here("data/analysis/seidr/network.rda"))
samples <- read_csv(here("doc/samples_B2.csv"))


plotEigengene(dat, "TRINITY_DN10144_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)   #NO #in transloc -> OK for scaffold
plotEigengene(dat, "TRINITY_DN11942_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold  #too low expression
plotEigengene(dat, "TRINITY_DN14682_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)   #NO #in transloc -> OK for scaffold, but 2nd path is borderline for the start
plotEigengene(dat, "TRINITY_DN15417_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)   #not found in gmap! #too low expression
plotEigengene(dat, "TRINITY_DN16354_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold #too low expression
plotEigengene(dat, "TRINITY_DN16731_c0_g1_i15",rep("bla", nrow(dat)),samples$Stages)  #NO #in mult -> OK for scaffold, but 1st path is borderline for start and 3rd for end
plotEigengene(dat, "TRINITY_DN17439_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)   #!!!  #in uniq -> OK for scaffold
plotEigengene(dat, "TRINITY_DN18384_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)   #NO #in transloc -> OK for scaffold
plotEigengene(dat, "TRINITY_DN22742_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold #too low expression
plotEigengene(dat, "TRINITY_DN22829_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)   #?   #in uniq -> OK for scaffold
plotEigengene(dat, "TRINITY_DN22966_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)   #YES  #not found in gmap!
plotEigengene(dat, "TRINITY_DN23658_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)   #NO #YES  #in transloc -> OK for scaffold
plotEigengene(dat, "TRINITY_DN30822_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)   #NO #in transloc -> OK for scaffold, not sure about 1st path
plotEigengene(dat, "TRINITY_DN30963_c0_g2_i2",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold #too low expression
plotEigengene(dat, "TRINITY_DN32364_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold #too low expression
plotEigengene(dat, "TRINITY_DN4106_c0_g1_i19",rep("bla", nrow(dat)),samples$Stages)   #YES  #in uniq -> OK for scaffold
plotEigengene(dat, "TRINITY_DN4327_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages)    #!!!  #in uniq -> OK for scaffold
plotEigengene(dat, "TRINITY_DN4994_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)    #NO #in mult -> OK for scaffold
plotEigengene(dat, "TRINITY_DN5317_c1_g1_i1",rep("bla", nrow(dat)),samples$Stages)    #NO #?   #in transloc -> OK for scaffold, check 1st path the end
plotEigengene(dat, "TRINITY_DN56194_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages)   #NO #in mult -> OK for scaffold, check 2nd and 3rd path for the end
plotEigengene(dat, "TRINITY_DN6342_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)    #in uniq -> OK for scaffold #too low expression


#plotEigengene(dat, "TRINITY_DN11855_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)  #NO #in mult -> OK for scaffold
#plotEigengene(dat, "TRINITY_DN12264_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)  #in uniq -> OK for scaffold, borderline the end #low expression
plotEigengene(dat, "TRINITY_DN1474_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)    #!!!  #not found in gmap!
plotEigengene(dat, "TRINITY_DN14880_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages)   #NO #?    #in mult -> OK for scaffold
#plotEigengene(dat, "TRINITY_DN15083_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)  #NO #in transloc, I think they are too close to the end of the scaffold
#plotEigengene(dat, "TRINITY_DN15486_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)  #in uniq -> OK for scaffold #low expression
#plotEigengene(dat, "TRINITY_DN17064_c0_g3_i1",rep("bla", nrow(dat)),samples$Stages)  #in uniq -> OK for scaffold #low expression
#plotEigengene(dat, "TRINITY_DN17445_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)  #in uniq, I think it is too close to the end of the scaffold #low expression
#plotEigengene(dat, "TRINITY_DN19939_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)  #in uniq -> OK for scaffold #low expression
#plotEigengene(dat, "TRINITY_DN19939_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)  #in uniq -> OK for scaffold #low expression
#plotEigengene(dat, "TRINITY_DN22795_c0_g2_i3",rep("bla", nrow(dat)),samples$Stages)  #in uniq -> OK for scaffold #low expression
##plotEigengene(dat, "TRINITY_DN30747_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages) #not found in gmap! #low expression
#plotEigengene(dat, "TRINITY_DN3092_c0_g2_i4",rep("bla", nrow(dat)),samples$Stages)   #NO #in mult -> OK for scaffold
plotEigengene(dat, "TRINITY_DN3337_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)    #!!!  #in uniq -> OK for scaffold
plotEigengene(dat, "TRINITY_DN39226_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)   #?    #in uniq -> OK for scaffold
#plotEigengene(dat, "TRINITY_DN3985_c0_g1_i9",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold #low expression
#plotEigengene(dat, "TRINITY_DN4216_c0_g1_i10",rep("bla", nrow(dat)),samples$Stages)  #not found in gmap! #low expression
plotEigengene(dat, "TRINITY_DN4293_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)    #?    #in uniq -> OK for scaffold
#plotEigengene(dat, "TRINITY_DN4297_c0_g1_i13",rep("bla", nrow(dat)),samples$Stages)  #NO #in transloc, 2nd path is at the and of the scaffold
#plotEigengene(dat, "TRINITY_DN4543_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold #low expression
#plotEigengene(dat, "TRINITY_DN6070_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)    #?    #in uniq, but it is at the and of the scaffold
#plotEigengene(dat, "TRINITY_DN8704_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold #low expression
#plotEigengene(dat, "TRINITY_DN8845_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)   #in uniq -> OK for scaffold, check the beginning #low expression


# checking the absolute expression

raw_data <- read_csv(here("data/analysis/salmon/raw-unormalised-gene-expression_data_linc.csv"))


# checking the candidates are not too close to the start or ending of a scaffold (remember the 2nd col is the one with the length of the scaffold)

all_scaffold <- read_tsv("/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome.fa.fai.gz",col_names = FALSE)

# get gmap files from other_tibble.R and check on them

start_distance <- c((5778 -1 + 1)/11990, (5390 -1 +1)/200048, (26454 -1 +1)/56540, (15415 -1 +1)/30584, (135 -1 +1)/12455, (17375 -1 +1)/23244, (589 -1 +1)/20040, (1990 -1 +1)/19580,
                    (15554 -1 +1)/16518, (928 -1 +1)/10785, (1063 -1 +1)/20953, (14233 -1 +1)/16358, (8573 -1 +1)/28295, (20058 -1 +1)/28295, (3119 -1 +1)/24032, (9109 -1 +1)/18228,
                    (7195 -1 +1)/11330, (5621 -1 +1)/16483, (34233 -1 +1)/35773, (29129 -1 +1)/47320, (5208 -1 +1)/32280, (9732 -1 +1)/29352, (7390 -1 +1)/25351, (7573 -1 +1)/12030,
                    (903 -1 +1)/22786, (8350 -1 +1)/14462, (7106 -1 +1)/18908, (19633 -1 +1)/28277, (42149 -1 +1)/63097, (49864 -1 +1)/51619, (8863 -1 +1)/10597, (6096 -1 +1)/12469,
                    (3882 -1 +1)/10032, (2849 -1 +1)/15320, (17098 -1 +1)/28838, (3250 -1 +1)/11564, (14078 -1 +1)/17630, (2815 -1 +1)/40558, (26100 -1 +1)/43506, (18634 -1 +1)/19003,
                    (18033 -1 +1)/19003, (3190 -1 +1)/9819, (7384 -1 +1)/22466, (15356 -1 +1)/18316, (10126 -1 +1)/19564, (10794 -1 +1)/19564, (21866 -1 +1)/57903, (10263 -1 +1)/32519,
                    (11700 -1 +1)/21344, (4394 -1 +1)/16616, (8414 -1 +1)/29691, (6059 -1 +1)/56477, (21001 -1 +1)/30292, (10636 -1 +1)/23599, (14631 -1 +1)/17025, (5252 -1 +1)/63173,
                    (11885 -1 +1)/14434, (5029 +1 -1)/11950, (1113 -1 +1)/10332)
end_distance <-  c((11990 -6876 +1)/11990, (200048 -7849 +1)/200048, (56540 - 28425 +1)/56540, (30584 -15882 +1)/30584, (12455 -358 +1)/12455, (23244 -17809 +1)/23244,
                   (20040 -1176 +1)/20040, (19580 -2632 +1)/19580, (16518 -16205 +1)/16518, (10785 -1640 +1)/10785, (20953 -1314 +1)/20953, (16358 -14558 +1)/16358,
                   (28295 -19656 +1)/28295, (28295 -21678 +1)/28295, (24032 -4003 +1)/24032, (18228 -13683 +1)/18228, (11330 -7514 +1)/11330, (16483 -5771 +1)/16483,
                   (35773 -34821 +1)/35773, (47320 -29488 +1)/47320, (32280 -6830 +1)/32280, (29352 -9812 +1)/29352, (25351 -10175 +1)/25351, (12030 -8630 +1)/12030,
                   (22786 -5278 +1)/22786, (14462 -8447 +1)/14462, (18908 -16034 +1)/18908, (28277 -19723 +1)/28277, (63097 -48699 +1)/63097, (51619 -50667 +1)/51619,
                   (10597 -9666 +1)/10597, (12469 -6347 +1)/12469, (10032 -4457 +1)/10032, (15320 -3755 +1)/15320, (28838 -24218 +1)/28838, (11564 -4413 +1)/11564,
                   (17630 -16412 +1)/17630, (40558 -9168 +1)/40558, (43506 -32467 +1)/43506, (19003 -18825 +1)/19003, (19003 -18755 +1)/19003, (9819 -5664 +1)/9819,
                   (22466 -7668 +1)/22466, (18316 -17903 +1)/18316, (19564 -12607 +1)/19564, (19564 -12607 +1)/19564, (57903 -24685 +1)/57903, (32519 -15955 +1)/32519,
                   (21344 -12331 +1)/21344, (16616 -1 +1)/16616, (29691 -9107 +1)/29691, (56477 -6815 +1)/56477, (30292 -21846 +1)/30292, (23599 -11181 +1)/23599,
                   (17025 -17021 +1)/17025, (63173 -7396 +1)/63173, (14434 -14434 +1)/14434, (11950 -5954 +1)/11950, (10332 -1559 +1)/10332)

plot(density(start_distance))
plot(density(end_distance))
source(here("UPSCb-common/src/R/percentile.R"))
percentile(end_distance)
abline(v=quantile(end_distance,probs=seq(.7,1,0.03),lty=2))
abline(v=quantile(start_distance,probs=seq(0,.3,0.03),lty=2))     

### adding a few guys
plotEigengene(dat, "TRINITY_DN7954_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages,colors = c("#E5C494")) #non male #in uniq
#plotEigengene(dat, "TRINITY_DN7987_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN6930_c1_g1_i1",rep("bla", nrow(dat)),samples$Stages) #boh?
plotEigengene(dat, "TRINITY_DN6930_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages) #??????
plotEigengene(dat, "TRINITY_DN6998_c2_g1_i1",rep("bla", nrow(dat)),samples$Stages) #mult!
#plotEigengene(dat, "TRINITY_DN59945_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN5946_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages) #non male #in uniq
plotEigengene(dat, "TRINITY_DN49482_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages) #mult
plotEigengene(dat, "TRINITY_DN49596_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages) #not found in gmap! #low expression
plotEigengene(dat, "TRINITY_DN4968_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages) #not found in gmap! #not bad expression
plotEigengene(dat, "TRINITY_DN4968_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages) #transloc
plotEigengene(dat, "TRINITY_DN3958_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages) #low expression #in uniq
plotEigengene(dat, "TRINITY_DN29128_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages) #low expression #in uniq
plotEigengene(dat, "TRINITY_DN29475_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages) #not found in gmap! #low expression
plotEigengene(dat, "TRINITY_DN2961_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages) #not found in gmap! #low expression
plotEigengene(dat, "TRINITY_DN19125_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages) #in mult
plotEigengene(dat, "TRINITY_DN19125_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages) #low expression #in uniq
plotEigengene(dat, "TRINITY_DN19125_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages) #in mult
plotEigengene(dat, "TRINITY_DN19201_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages) #in mult
plotEigengene(dat, "TRINITY_DN19262_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages) #not found in gmap! #so and so expression #avoid
#plotEigengene(dat, "TRINITY_DN19654_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plotEigengene(dat, "TRINITY_DN39015_c0_g2_i2",rep("bla", nrow(dat)),samples$Stages) #in transloc

plotEigengene(dat, "TRINITY_DN28175_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages) #in uniq #not bad expression
plotEigengene(dat, "TRINITY_DN3892_c2_g1_i1",rep("bla", nrow(dat)),samples$Stages) #in uniq #the expression can work
plotEigengene(dat, "TRINITY_DN4877_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages) #in uniq #never saw such a high expression. That's the guy!
plotEigengene(dat, "TRINITY_DN4877_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages) #in uniq #pretty high expression, but the other isoform is better
plotEigengene(dat, "TRINITY_DN4968_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages) #in uniq #not bad expression
plotEigengene(dat, "TRINITY_DN58421_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages) #in uniq #not too bad expression 
#plotEigengene(dat, "TRINITY_DN6816_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages) #not in uniq #the expression is good


plotEigengene(dat, "miRNA_5305-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_5305-5p")
ggsave(filename=here("data/analysis/figures_new//miRNA2.png"),device ="png",dpi = 600)


bla <- selected_linc[order(selected_linc$PredictionScore,decreasing = TRUE), ]
bla2col <- bla[,1:2]
bla <- as_tibble(bla)
bla2col <- bla[,1:2]
bla_30 <- bla2col[1:30,]
#final guys selected
TRINITY_DN17439_c0_g1_i4 #also in Infomap cluster2
TRINITY_DN22829_c0_g1_i1 #not in Infomap clusters
TRINITY_DN22966_c0_g2_i1 #not in Infomap clusters
TRINITY_DN4106_c0_g1_i19 #also in Infomap cluster2
TRINITY_DN4327_c0_g1_i8  #also in Infomap cluster, but not in the first 8
TRINITY_DN1474_c0_g1_i2  #also in Infomap cluster5
TRINITY_DN3337_c0_g1_i6  #also in Infomap cluster2
TRINITY_DN39226_c0_g1_i1 #also in Infomap cluster5
TRINITY_DN4293_c0_g1_i1  #also in Infomap cluster2
TRINITY_DN4877_c0_g1_i1  #also in Infomap cluster3
TRINITY_DN4968_c0_g1_i3  #also in Infomap cluster3
TRINITY_DN58421_c0_g1_i2 #also in Infomap cluster3
TRINITY_DN5946_c0_g1_i2  #also in Infomap cluster3

#enrichment of 1st degree neighbours of the 6 final candidates
getGeneFDN <- function(edgeList, gene, source.col=1, target.col=2) {
  # TODO check for data type
  # TODO check that all genes are in the edgelist
  
  s2t <- edgeList[edgeList[source.col] == gene,][,target.col]
  t2s <- edgeList[edgeList[target.col] == gene,][,source.col]
  res <- union(s2t,t2s)
  
  return(res)
}
edgeList <- read.table(here("data/seidr/backbone/edgelist2.txt"))
candidate1_FDN <- getGeneFDN(edgeList,"TRINITY_DN17439_c0_g1_i4")
c1 <- plotEigengene(dat, "TRINITY_DN17439_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages)
plot(c1)
ggsave(filename=here("data/analysis/figures_new//candidate1.png"),device ="png",dpi = 600)
MA_19918g0010 #in cluster4           
TRINITY_DN4994_c0_g1_i3 #in cluster2
MA_10426731g0010 #in cluster2       
MA_10428011g0010 #in cluster2      
MA_10432599g0010 #in cluster2        
MA_10433793g0010 #in cluster2      
MA_122121g0010 #in cluster2          
MA_18112g0010  #in cluster2         
MA_480751g0010 #in cluster2         
TRINITY_DN18384_c0_g1_i4 #in cluster1
TRINITY_DN21012_c0_g1_i1 #in cluster2

candidate2_FDN <- getGeneFDN(edgeList,"TRINITY_DN22966_c0_g2_i1")
c2 <- plotEigengene(dat, "TRINITY_DN22966_c0_g2_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c2)
ggsave(filename=here("data/analysis/figures_new//candidate2.png"),device ="png",dpi = 600)

candidate3_FDN <- getGeneFDN(edgeList,"TRINITY_DN3337_c0_g1_i6")
c3 <- plotEigengene(dat, "TRINITY_DN3337_c0_g1_i6",rep("bla", nrow(dat)),samples$Stages)
plot(c3)
ggsave(filename=here("data/analysis/figures_new//candidate3.png"),device ="png",dpi = 600)
TRINITY_DN28199_c0_g1_i1 #in cluster4
MA_120984g0010 #in cluster2
TRINITY_DN6988_c0_g1_i2 #in cluster2
TRINITY_DN7346_c0_g1_i1 #in cluster2
#miRNA_14300-3p #in cluster2

candidate4_FDN <- getGeneFDN(edgeList,"TRINITY_DN4293_c0_g1_i1")
c4 <- plotEigengene(dat, "TRINITY_DN4293_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c4)
ggsave(filename=here("data/analysis/figures_new//candidate4.png"),device ="png",dpi = 600)
TRINITY_DN11909_c1_g1_i5 #in cluster1
MA_18951g0030 #in cluster2
MA_8491653g0010 #in cluster2

candidate5_FDN <- getGeneFDN(edgeList,"TRINITY_DN4877_c0_g1_i1")
c5 <- plotEigengene(dat, "TRINITY_DN4877_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c5)
ggsave(filename=here("data/analysis/figures_new//candidate5.png"),device ="png",dpi = 600)
MA_10432090g0010 #in cluster3
MA_130040g0010 #in cluster3
MA_24589g0010 #in cluster3
MA_7182809g0010 #in cluster3
MA_9918531g0010 #in cluster3
MA_10425895g0010 #in cluster3
TRINITY_DN21921_c0_g1_i6 #in cluster3
TRINITY_DN49482_c0_g1_i1 #in cluster3

candidate6_FDN <- getGeneFDN(edgeList,"TRINITY_DN5946_c0_g1_i2")
c6 <- plotEigengene(dat, "TRINITY_DN5946_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plot(c6)
ggsave(filename=here("data/analysis/figures_new//candidate6.png"),device ="png",dpi = 600)
MA_102205g0010 #in cluster3
MA_5475932g0010 #in cluster3

 en6 <- gopher(genes = candidate6_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")
the_guy <- gopher(genes = "MA_125713g0010",task=c("mapman"),url="pabies",endpoint = "gene-to-term")
the_guy_FDN <- getGeneFDN(edgeList,"MA_125713g0010")
en_guy_FDN <- gopher(genes = the_guy_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")
library(RColorBrewer)
#trying to do a multiline plot
plotEigengene(dat, c("MA_10435426g0010","MA_10427079g0030","MA_70814g0010","MA_10434973g0010","MA_10428192g0010","MA_9486006g0010"), rep("bla", nrow(dat)),samples$Stages,multiline = TRUE,
              title="miRNA_44911-5p target genes")
ggsave(filename=here("data/analysis/figures_new//miRNA1_targets-2.png"),device ="png",dpi = 600)
plotEigengene(dat, "MA_10435426g0010",rep("bla", nrow(dat)),samples$Stages)

plotEigengene(dat, c("MA_10435426g0010","MA_10427079g0030","MA_10434973g0010","MA_9486006g0010"), rep("bla", nrow(dat)),samples$Stages,multiline = TRUE,
              title = "miRNA_5305-5p target genes")
ggsave(filename=here("data/analysis/figures_new//miRNA2_targets.png"),device ="png",dpi = 600)
plotEigengene(dat, c("MA_163067g0010","MA_10433283g0020","MA_24916g0010","MA_28575g0010",
                     "MA_10432894g0010","MA_10135850g0010","MA_941687g0010","MA_128563g0010",
                     "MA_10380652g0010","MA_724730g0010","MA_10429833g0010","MA_2249g0010",
                     "MA_10231644g0010","MA_99694g0010","MA_20215g0010","MA_38220g0010",
                     "MA_45741g0010","MA_20447g0010","MA_116129g0010"), rep("bla", nrow(dat)),samples$Stages,multiline = TRUE)


plotEigengene(dat, c("MA_10433283g0020","MA_24916g0010","MA_941687g0010","MA_724730g0010"), rep("bla", nrow(dat)),samples$Stages,multiline = TRUE)
              
MA_163067g0010
MA_10433283g0020
MA_24916g0010
MA_28575g0010
MA_10432894g0010
MA_10135850g0010
MA_941687g0010
MA_128563g0010
MA_10380652g0010
MA_724730g0010
MA_10429833g0010
MA_2249g0010
MA_10231644g0010
MA_99694g0010
MA_20215g0010
MA_38220g0010
MA_45741g0010
MA_20447g0010
MA_116129g0010

MA_10433283g0020
MA_24916g0010
MA_941687g0010
MA_724730g0010

the_guy2 <- gopher(genes = "MA_5475932g0010",task=c("mapman"),url="pabies",endpoint = "gene-to-term")
the_guy2_FDN <- getGeneFDN(edgeList,"MA_5475932g0010")
en_guy2_FDN <- gopher(genes = the_guy2_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")

#other candidates
candidate7_FDN <- getGeneFDN(edgeList,"TRINITY_DN22829_c0_g1_i1")
MA_10276622g0010 #as well not in Infomap clusters
MA_10434511g0010 #as well not in Infomap clusters
c7 <- plotEigengene(dat, "TRINITY_DN22829_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c7)
ggsave(filename=here("data/analysis/figures_new//candidate7.png"),device ="png",dpi = 600)

candidate8_FDN <- getGeneFDN(edgeList,"TRINITY_DN4106_c0_g1_i19")
MA_10427558g0010 #in cluster2
MA_10432863g0010 #in cluster2
MA_26708g0010    #in cluster2
TRINITY_DN37936_c0_g1_i3 #in cluster2
TRINITY_DN56194_c0_g1_i7 #in cluster2
en8 <- gopher(genes = candidate8_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")
c8 <- plotEigengene(dat, "TRINITY_DN4106_c0_g1_i19",rep("bla", nrow(dat)),samples$Stages)
plot(c8)
ggsave(filename=here("data/analysis/figures_new//candidate8.png"),device ="png",dpi = 600)

candidate9_FDN <- getGeneFDN(edgeList,"TRINITY_DN4327_c0_g1_i8")
MA_99816g0010 #not in the first 8 clusters
MA_893876g0010 #not in the first 8 clusters
TRINITY_DN30963_c0_g2_i2 #not in the first 8 clusters
en9 <- gopher(genes = candidate9_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")
c9 <- plotEigengene(dat, "TRINITY_DN4327_c0_g1_i8",rep("bla", nrow(dat)),samples$Stages)
plot(c9)
ggsave(filename=here("data/analysis/figures_new//candidate9.png"),device ="png",dpi = 600)

candidate10_FDN <- getGeneFDN(edgeList,"TRINITY_DN1474_c0_g1_i2")
MA_245940g0010 #in cluster5
MA_358441g0010 #in cluster5
MA_4911g0010 #in cluster5
MA_10432733g0020 #in cluster2
MA_10435887g0040 #in cluster5
MA_524371g0010 #in cluster5
MA_8933673g0010 #in cluster5
en10 <- gopher(genes = candidate10_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")
c10 <- plotEigengene(dat, "TRINITY_DN1474_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plot(c10)
ggsave(filename=here("data/analysis/figures_new//candidate10.png"),device ="png",dpi = 600)

candidate11_FDN <- getGeneFDN(edgeList,"TRINITY_DN39226_c0_g1_i1")
MA_10435757g0010 #in cluster5
MA_191579g0010 #in cluster5
MA_66681g0010 #in cluster4
TRINITY_DN28199_c0_g1_i1 #in cluster4
#miRNA_5003-3p #in cluster5
en11 <- gopher(genes = candidate11_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")
c11 <- plotEigengene(dat, "TRINITY_DN39226_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages)
plot(c11)
ggsave(filename=here("data/analysis/figures_new//candidate11.png"),device ="png",dpi = 600)

candidate12_FDN <- getGeneFDN(edgeList,"TRINITY_DN4968_c0_g1_i3")
MA_412851g0010 #in cluster3
c12 <- plotEigengene(dat, "TRINITY_DN4968_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages)
plot(c12)
ggsave(filename=here("data/analysis/figures_new//candidate12.png"),device ="png",dpi = 600)

candidate13_FDN <- getGeneFDN(edgeList,"TRINITY_DN58421_c0_g1_i2")
#we have 64 FDNs: 22 lincRNAs + 2 miRNAs
en13 <- gopher(genes = candidate13_FDN,background = InfomapClusters$gene,task = list("go","mapman","pfam","kegg"),alpha = 0.05,url="pabies",endpoint = "enrichment")
c13 <- plotEigengene(dat, "TRINITY_DN58421_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages)
plot(c13)
ggsave(filename=here("data/analysis/figures_new//candidate13.png"),device ="png",dpi = 600)




