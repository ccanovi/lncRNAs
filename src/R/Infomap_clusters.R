## SETUP
# load necessary librarieslibrary(data.table)
# load necessary libraries
library(data.table)
library(reshape2)
library(dplyr)
library(tidyverse)

## PREPARE THE DATA
# TASK - Load the seidrResolve file into the infomapTable variable
infomapTable <- read.delim("data/seidr/backbone/clustersResolve.txt", 
                           header = FALSE, stringsAsFactors=FALSE)

# TASK - Calculate the number of levels
# numberLevels =  total columns - 3 (info columns) - 1 (path column)
# Expected output: an integer
numberLevels <- length(colnames(infomapTable)) - 4

# TASK- create a of vector names for the levels
# ex. Level1, Level2, etc
# Expected output: [1] "Level1" "Level2" "Level3" "Level4" ...
levels <- sapply(1:numberLevels, function(x){
  paste0("Level", x)
})

# Add the names to the data frame
colnames(infomapTable) <- c('Path', levels, 'Flow', 'Index', 'Gene')

# create a new results data frame
res <- data.frame(Gene = infomapTable[, ncol(infomapTable)],
                  Path = infomapTable[, 1],
                  stringsAsFactors = F)

# CHECKPOINT
head(res)
# If everthing is ok the previous command should look similar to this
# These geneIDs are aspen genes, for spruce they would start with "PICAB"
#              Gene    Path
# 1   Potra2n1c2185 1:1:1:1
# 2  Potra2n5c11770 1:1:1:2
# 3 Potra2n16c30312 1:1:1:3
# 4 Potra2n15c28283 1:1:1:4

# Infomap data always have level 1 and none of the data is NA,
# we can add it directly
res$Level1 <- as.character(infomapTable[, 2])

# loop though the rest of the levels and attach the name from the previous ones
for (level in 2:numberLevels) {
  
  currentLevel <- paste0("Level", level)
  prevLevel <- paste0("Level",(level-1))
  
  # join names
  res[[currentLevel]] <- paste0(res[[prevLevel]], ":", infomapTable[[currentLevel]])
  
  # if there is an NA inside the current name, that gene doesn't belong to a cluster in that level, it turns into NA
  res[[currentLevel]] <- ifelse(res[[currentLevel]] %like% "NA", NA, res[[currentLevel]])  
}

## CLUSTERING LEVEL SELECTION
## Check the decision tree in the lecture material

# TASK - we start at Level1 if it doesn't match our criteria come back to this
# line and change the level to Level2 and so on until the 
# criteria has been reached
level <- 'Level1'

# How many cluster does this level has?
print(paste(level,"clusters:", length(unique(res[[level]])))) 

# Using dplyr we obtain the number of genes per cluster
counts <- res %>% 
  group_by(!!sym(level)) %>%
  summarise(geneNumber=n()) %>%
  arrange(desc(geneNumber))

# TASK - get the amount of clusters in the current level
# Use the unique and length commands
# Expected output: integer
clusterNumber <- unique(res$Level1)

# TASK - get amount of genes in the top 20 clusters
# Expected output: integer
genesInTop20 <- sum(counts[1:20,2])

# TASK - calculate the total amount of genes in the network
# The total genes in the network matches the rows of the data
# Expected output - integer
genesTotal <- sum(counts[,2])

# Let's print the information so far
print(paste("Genes in the top 20 clusters:", genesInTop20))
print(paste("Genes in the network", genesTotal))

# CHECKPOINT
# The following line will tell you if the current level matches our
# selection criteria.
print(paste("Percentage of genes in the top 20 clusters:", round(genesInTop20/genesTotal, digits = 4)*100))

# If the previous percentage doesn't match our criteria, come back
# to the CLUSTER LEVEL SELECTION section and increase the current 
# level in one unit

# If the criteria was matched then continue
# 
## CONSTRUCT THE cluster LIST
# Extract the gene names for each cluster of interest
# 
# TASK - subset the res matrix and keep only the Gene column
# and the current level column
# 
infomapClusters <- res[,c("Gene", "Level1")]

# CHECKPOINT
# 
# Let's check the current cluster list
head(infomapClusters)
# If everything is ok the previous command should give a result similar to this
# Gene Level1
# 1   Potra2n1c2185      1
# 2  Potra2n5c11770      1
# 3 Potra2n16c30312      1

# TASK - add the the string "Cluster" as prefix to the cluster number
# 
infomapClusters$Level1 <- paste0("Cluster", infomapClusters$Level1)

# change the column names to gene and cluster before exporting  
names(infomapClusters) <- c("gene", "cluster")


# CHECKPOINT
# 
# Let's check the current edgeList
head(infomapClusters)
# If everything is ok the previous command should give a result similar to this
#              gene  cluster
# 1   Potra2n1c2185 Cluster1
# 2  Potra2n5c11770 Cluster1
# 3 Potra2n16c30312 Cluster1

# FILTERING
# Lets filter the clusters and only keep those that have at least 15 genes.
minSize <- 15

infomapClusters.filter <- infomapClusters %>% 
  group_by(cluster) %>% 
  filter(n() >= minSize)

# Let's check the result
infomapClusters.filter %>% tally()


## SAVE THE RESULTS

# write the results in a tsv file
write.table(infomapClusters.filter,
            file = "~/Git/lncRNAs/data/seidr/backbone/infomapClusters.tsv",
            row.names = FALSE, sep = '\t')

write.csv(infomapClusters.filter,
       file = "~/Git/lncRNAs/data/seidr/backbone/infomapClusters.csv",
       row.names = FALSE)
