# libraries
library(dplyr)
library(here)
library(readr)

# data
dat <- read_tsv(here("doc/network_matrix.tsv")) %>% 
  tibble::column_to_rownames("ID") %>% as.data.frame()

# sample
smpl <- read_csv(here("doc/samples_final.csv")) %>% 
  filter(!duplicated(ID)) %>% select(ID,Stages)

# sanity
stopifnot(all(colnames(dat) == smpl$ID))

# calc medians
mat <- sapply(split.data.frame(t(dat),smpl$Stages),matrixStats::colMedians)
rownames(mat) <- rownames(dat)

# goi
goi <- c("miRNA_12146-5p", "TRINITY_DN10084_c0_g5_i1", "TRINITY_DN64598_c0_g1_i3",
"MA_49382g0010")

# permutations
res <- t(apply(combn(1:length(goi),2),2,function(inx){
  sp <- cor.test(mat[goi,][inx[1],],mat[goi,][inx[2],],method="spearman")
  pe <- cor.test(mat[goi,][inx[1],],mat[goi,][inx[2],],method="pearson")
  
  c(from=goi[inx[1]],
    to=goi[inx[2]],
    spearman=sp$estimate,
    spearman_pval=sp$p.value,
    pearson=pe$estimate,
    pearson_pval=pe$p.value)
}))


