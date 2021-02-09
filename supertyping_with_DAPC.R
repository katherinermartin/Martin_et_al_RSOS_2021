# Supertyping using DAPC to describe alleles in functional categories
# this script includes k-means clustering, 2-step cross validation, and evaluation of k 3 through 11 at optimized values.

library(adegenet)
library(ggplot2)
library(reshape2)
library(scales)
library(unikn)

setwd("~/Documents/UCF/Research/MHC_Class_I/Analysis/03.05.19_analysis/supertyping/Jan21/")

# Data frame of amino acid residue values for 124 MHC alleles (includes Martin 2021 alleles and those alleles from Stiebens et al. 2013)

sea_turtle_full <- read.csv("MHC_amino_acid_matrix_values_for_DAPC.csv") # load in data matrix; electrochemical values for each amino acid in the full pepetide binding region, based on 
sea_turtle_full <- sea_turtle_full[,-1] # remove first column (allele identifiers) so that it can be read into find.clusters (below)

# the following primer from the Grunwald lab was extremely helpful: http://grunwaldlab.github.io/Population_Genetics_in_R/

### K MEANS CLUSTERING

# "K" is the number of clusters to group the data into and so choosing a k is very important. First, perform K-means clustering over a number of values of K and repeat 10 times for each value to explore each k.

maxK <- 25 # argument for max.n.clust, which is an integer indicating the maximum number of clusters to be tried; find.clusters will evaluate 1 through k clusters (here, 11)

myMat <- matrix(nrow=100, ncol=maxK) # create empty matrix to be filled by the for loop below
colnames(myMat) <- 1:ncol(myMat) # give column names to matrix


for(i in 1:nrow(myMat)){
  grp <- find.clusters(sea_turtle_full, n.pca = 50, choose.n.clust = FALSE,  max.n.clust = maxK) # retains 50 PCAs (n.pca), user doesn't choose cluster number (choose.n.clust = FALSE), and the maximum number of clusters to be tried (max.n.clust) is 11 (so it does 1-11)
  myMat[i,] <- grp$Kstat # fill matrix with Kstats from groups
}

# Visualizing k means clustering

my_df <- melt(myMat) # turn matrix into a df
colnames(my_df)[1:3] <- c("Group", "K", "BIC") # column headings for df
my_df$K <- as.factor(my_df$K) # make the K values a factor
head(my_df) # lists the groups per each K and the BIC

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1 # shows how BIC changes with the number of groups (K) that's chosen

## There is no clear elbow point of minimizes BIC values, but 3-11 seems to encompass the majority of the change, so evaluate k = 3:11

### CROSS-VALIDATION
# This is a two-step validation procedure, where the number of principle components to retain is evaluated at each k 3:11
# Prior to validation, run find.clusters to get the group membership (required for the validation steps)
# Then, run a broad validation step using xvalDAPC, where 1-50 principle components are retained for 30 iterations. This will yield a principle component number that minimizes the mean square error (and is therefore a good choice for that k value)
# Finally, run a narrower validation step using xvalDAPC, but this time instead of retaining 1-50 principle components, center the calculation on the result from above. For example, if PC = 10 was the value that minimized MSE in the first validation, then run the second validation on PCs 1 to 20 (where 10 is the center of that distribution), and run for 1000 iterations.
# The result will be the number of PCs to retain at a particular k, which can then be used to run DAPC().


# k = 3
# Set up group membership
grp_3 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 3) # evaluated to 25 clusters; retain 50 PCs. at k = 3 clusters

set.seed(1)
xval_3 <- xvalDapc(sea_turtle_full, grp = grp_3$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs asociated iwth the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.

xval_3[2:6] # 10

# 3_a: evaluate 10 as center of optimal distribution:

set.seed(1)
xval_3_a <- xvalDapc(sea_turtle_full, grp = grp_3$grp, n.pca = 1:20, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 20 PCAs, since 10 is the center.

xval_3_a[2:6] # 14 is optimal number of PCAs at 3 clusters for correctly predicting subsamples with the lowest error.

# Result: 3 clusters, 14 PCs

dapc_3 <- dapc(sea_turtle_full, grp_3$grp, n.pca = 14) # retain 14 PCs automatically; retain 2 DFs

# Visualize DAPC results

pal3 <- seecol(pal_unikn_pref, n = 3) # set up palette

k3_scatter <- scatter(dapc_3, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal3)


# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

# Visualize:
compo3 <- compoplot(dapc_3, col = pal3)

##################################################


# k = 4
# Set up group membership
grp_4 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 4) # evaluated to 25 clusters; retain 50 PCs. at k = 4 clusters

set.seed(1)
xval_4 <- xvalDapc(sea_turtle_full, grp = grp_4$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_4[2:6] # 30 PCs

# 4_a: center xval on 30 PCs:
set.seed(1)
xval_4_a <- xvalDapc(sea_turtle_full, grp = grp_4$grp, n.pca = 20:40, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_4_a[2:6] # 28 PCs

# Result: 4 clusters, 28 PCs

dapc_4 <- dapc(sea_turtle_full, grp_4$grp, n.pca = 28) # retain 12 PCs automatically; retain 3 DFs

pal4 <- seecol(pal_unikn_pref, n = 4)

k4_scatter <- scatter(dapc_4, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal4)


temp90 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

compo4 <- compoplot(dapc_4, col = pal4)

###################################

# k = 5
# Set up group membership
grp_5 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 5) # evaluated to 25 clusters; retain 50 PCs. at k = 5 clusters

set.seed(1)
xval_5 <- xvalDapc(sea_turtle_full, grp = grp_5$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_5[2:6] # 10 PCs

# 5_a: evaluate  as center of optimal distribution:
set.seed(1)
xval_5_a <- xvalDapc(sea_turtle_full, grp = grp_5$grp, n.pca = 1:20, n.rep = 1000, parallel = "multicore", ncpus = 4L) 

xval_5_a[2:6] # 4 is optimal number of PCAs at 5 clusters for correctly predicting subsamples with the lowest error.

# Result: 5 clusters, 4 PCs

dapc_5 <- dapc(sea_turtle_full, grp_5$grp, n.pca = 4) # retain 4 PCs automatically; retain 4 DFs

pal5 <- seecol(pal_unikn_pref, n = 5)

k5_scatter <- scatter(dapc_5, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal5)


temp90 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 5 alleles

temp95 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 7 alleles

temp99 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 9 alleles

compoplot(dapc_5, col = pal5)

###################################

# k = 6

# Set up group membership
grp_6 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 6) # evaluated to 25 clusters; retain 50 PCs. at k = 6 clusters

set.seed(1)
xval_6 <- xvalDapc(sea_turtle_full, grp = grp_6$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_6[2:6] # 15 PCs

# 6_a: evaluate 15 as center of optimal distribution:
set.seed(1)
xval_6_a <- xvalDapc(sea_turtle_full, grp = grp_6$grp, n.pca = 5:25, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_6_a[2:6] # 15 is optimal number of PCAs at 6 clusters for correctly predicting subsamples with the lowest error.

# Result: 6 clusters, 15 PCs

dapc_6 <- dapc(sea_turtle_full, grp_6$grp, n.pca = 15) # retain 15 PCs automatically; retain 5 DFs

pal6 <- seecol(pal_unikn_pref, n = 6)

k6_scatter <- scatter(dapc_6, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = seecol(pal_unikn_pref, n = 6))


temp90 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

compo6 <- compoplot(dapc_6, col = pal6)


###################################

# k = 7

# Set up group membership
grp_7 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 7) # evaluated to 25 clusters; retain 50 PCs. at k = 7 clusters

set.seed(1)
xval_7 <- xvalDapc(sea_turtle_full, grp = grp_7$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_7[2:6] # 20 PCs

# 7_a: evaluate 20 as center of optimal distribution:
set.seed(1)
xval_7_a <- xvalDapc(sea_turtle_full, grp = grp_7$grp, n.pca = 10:30, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_7_a[2:6] # 17 is optimal number of PCAs at 7 clusters for correctly predicting subsamples with the lowest error.

# Result: 7 clusters, 17 PCs

dapc_7 <- dapc(sea_turtle_full, grp_7$grp, n.pca = 17) # retain 17 PCs automatically; retain 6 DFs

pal7 <- seecol(pal_unikn_pref, n = 7)

k7_scatter <- scatter(dapc_7, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal7)

temp90 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

compo7 <- compoplot(dapc_7, col = pal7)


###################################

# k = 8

# Set up group membership
grp_8 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 8) # evaluated to 25 clusters; retain 50 PCs. at k = 8 clusters

set.seed(1)
xval_8 <- xvalDapc(sea_turtle_full, grp = grp_8$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_8[2:6] # 20 PCs

# 8_a: evaluate 20 as center of optimal distribution:
set.seed(1)
xval_8_a <- xvalDapc(sea_turtle_full, grp = grp_8$grp, n.pca = 10:30, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_8_a[2:6] # 20 is optimal number of PCAs at 8 clusters for correctly predicting subsamples with the lowest error.

# Result: 8 clusters, 20 PCs

dapc_8 <- dapc(sea_turtle_full, grp_8$grp, n.pca = 20) # retain 20 PCs automatically; retain 7 DFs

pal8 <- seecol(pal_unikn_pref, n = 8)

k8_scatter <- scatter(dapc_8, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal8)


temp90 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 2 alleles

temp95 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 2 alleles

temp99 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 2 alleles

compo_8 <- compoplot(dapc_8, col = pal8)


###################################

# k = 9

# Set up group membership
grp_9 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 9) # evaluated to 25 clusters; retain 50 PCs. at k = 9 clusters

set.seed(1)
xval_9 <- xvalDapc(sea_turtle_full, grp = grp_9$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_9[2:6] # 20 PCs

# 9_a: evaluate 20 as center of optimal distribution:
set.seed(1)
xval_9_a <- xvalDapc(sea_turtle_full, grp = grp_9$grp, n.pca = 10:30, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_9_a[2:6] # 21 is optimal number of PCAs at 9 clusters for correctly predicting subsamples with the lowest error.

# Result: 9 clusters, 21 PCs

dapc_9 <- dapc(sea_turtle_full, grp_9$grp, n.pca = 21) # retain 21 PCs automatically; retain 8 DFs

pal9 <- seecol(pal_unikn_pref, n = 9)

k9_scatter <- scatter(dapc_9, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal9)


temp90 <- which(apply(dapc_9$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_9$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_9$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

compo_9 <- compoplot(dapc_9, col = pal9)


###################################

# k = 10

# Set up group membership
grp_10 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 10) # evaluated to 25 clusters; retain 50 PCs. at k = 10 clusters
set.seed(1)
xval_10 <- xvalDapc(sea_turtle_full, grp = grp_10$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_10[2:6] # 20 PCs

# 10_a: evaluate 20 as center of optimal distribution:
set.seed(1)
xval_10_a <- xvalDapc(sea_turtle_full, grp = grp_10$grp, n.pca = 10:30, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_10_a[2:6] #  17 is optimal number of PCAs at 10 clusters for correctly predicting subsamples with the lowest error.

# Result: 10 clusters, 17 PCs

dapc_10 <- dapc(sea_turtle_full, grp_10$grp, n.pca = 17) # retain x PCs automatically; retain 9 DFs

pal10 <- seecol(pal_unikn_pref, n = 10)

k10_scatter <- scatter(dapc_10, # data
                       bg = "white", # white background
                       pch = 20, # size and shape of points
                       cstar = 0, # cstar= lines btwn points, 0 for null
                       solid = 0.6,
                       cex = 1.5, #
                       clab = 0.375,
                       leg = TRUE,
                       scree.da = TRUE,
                       scree.pca = TRUE,
                       col = pal10)

temp90 <- which(apply(dapc_10$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 2 alleles

temp95 <- which(apply(dapc_10$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 3 alleles

temp99 <- which(apply(dapc_10$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 5 alleles

compo_10 <- compoplot(dapc_10, col = pal10)

###################################

# k = 11

# Set up group membership
grp_11 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 11) # evaluated to 25 clusters; retain 50 PCs. at k = 11 clusters
set.seed(1)
xval_11 <- xvalDapc(sea_turtle_full, grp = grp_11$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_11[2:6] # 20 PCs

# 11_a: evaluate 20 as center of optimal distribution:
set.seed(1)
xval_11_a <- xvalDapc(sea_turtle_full, grp = grp_11$grp, n.pca = 10:30, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_11_a[2:6] #  17 is optimal number of PCAs at 11 clusters for correctly predicting subsamples with the lowest error.

# Result: 11 clusters, 17 PCs

dapc_11 <- dapc(sea_turtle_full, grp_11$grp, n.pca = 17) # retain 17 PCs automatically; retain 10 DFs

pal11 <- seecol(pal_unikn_pref, n = 11)

k11_scatter <- scatter(dapc_11, # data
                       bg = "white", # white background
                       pch = 20, # size and shape of points
                       cstar = 0, # cstar= lines btwn points, 0 for null
                       solid = 0.6,
                       cex = 1.5, #
                       clab = 0.375,
                       leg = TRUE,
                       scree.da = TRUE,
                       scree.pca = TRUE,
                       col = pal11)

temp90 <- which(apply(dapc_11$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 2 alleles

temp95 <- which(apply(dapc_11$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 2 alleles

temp99 <- which(apply(dapc_11$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 3 alleles

compo_11 <- compoplot(dapc_11, col = pal11)


### K = 3 through 11 has been evaluated and 3 appears to be the optimal cluster number, based on the scatter plots and composition plot information.

# Generate a dataframe of each allele and its supertype
supertypes <- as.data.frame(dapc_3[["grp"]]) # gives supertype membership of each allele.
supertypes

allele_names <- read.csv("MHC_amino_acid_matrix_values_for_DAPC.csv")
allele_names <- allele_names[1]

MHC_alleles_by_supertype <- cbind(supertypes, allele_names)

write.csv(MHC_alleles_by_supertype, "MHC_alleles_by_supertype_final.csv") # lists each allele and its supertype (groups here are named 1,2,3; in manuscript they are named A,B,C)




