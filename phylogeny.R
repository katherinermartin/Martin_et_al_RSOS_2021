# Script to create figures for Bayesian and Maximum Likelihood phylogenies
library(ggtree)
library(treeio)
library(ape)
library(ips)
library(ggplot2)
library(phytools)
library(svglite)
library(tidyr)
library(tidyverse)

# read in data

# Bayesian tree:
tree <-read.beast("classI_nucleotide_alm.nex.con.tre")
tree

data <- read.csv("allele_counts_by_species.csv") # individuals per allele, by site (Florida or Cape Verde) and species (C. mydas and C. caretta)

data$Cm[which(data$Cm > 0)] <- "Cm"
data$Cc[which(data$Cc > 0)] <- "Cc"
data$Cc_CV[which(data$Cc_CV > 0)] <- "Cc_CV"
data[data == 0] <- NA

data$Cm <- as.factor(data$Cm)
data$Cc <- as.factor(data$Cc)
data$Cc_CV <- as.factor(data$Cc_CV)


# Initialize tree
p <- ggtree(tree, ladderize= TRUE, right = TRUE)
p

# tree with colored labels based upon species and location
p2 <- p %<+% data +
  geom_tippoint(color = "#008000", fill = "#008000", aes(x=x+0.01, subset = Cm =="Cm"), size = 1) +
  geom_tippoint(color = "#ff8c00", fill = "#ff8c00", aes(x=x+0.022, subset = Cc =="Cc"), size = 1) +
  geom_tippoint(color = "#663399", fill ="#663399", aes(x=x+0.034, subset = Cc_CV =="Cc_CV"), size = 1) +
  geom_treescale(width = 0.2)

# add internal nodes colored by posterior probability
p3 <- p2 +
  geom_nodepoint(color = "black", fill = "black", aes(subset = prob >= 0.995), size = 1, shape = 21) + #posterior prob between 0.995- 1
  geom_nodepoint(color = "black", fill = "#808080", aes(subset = prob < 0.995 & prob >= 0.945), size = 1, shape = 21) + #posterior prob between 0.945-0.994
  geom_nodepoint(color = "black", fill = "#dcdcdc", aes(subset = prob < 0.945 & prob >= 0.895), size = 1, shape = 21) # posterior prob between 0.895-0.944
p3

# optional: above tree but with tip labels
p4 <- p3 + geom_tiplab(color="black", size = 1.5, hjust=-.6) # adds tip labels
p4 # tree with correctly colored nodes and labels too.



# Maximum likelihood tree:

tree <- read.raxml("RAxML_bipartitionsBranchLabels.classI_alm_raxml")

# label bootstrap values 90 and above
p <- ggtree(tree, ladderize= TRUE, right = TRUE) +
  geom_tiplab(color="black", size = 2, hjust = -0.2) +
  geom_nodepoint(color = "black", fill = "black", aes(subset = bootstrap >= 90), size = 1, shape = 21)
p # bootstrap labels above 90 labeled with solid black dot
