## Classification random forest for tumor texture in C. mydas with FP

library(randomForest)
library(caret)
library(dplyr)

# read in data
data_full <- read.csv("classI_juveniles_morpho_FP_v3.csv") # this includes C. mydas and C. caretta records

data_full <- data_full %>% filter(!is.na(Carap_L_sl)) # remove any record without straight carapace length (SCL) measurement (just 1)

# variables: sample ID, location, spp, year, FP, paps_texture, season, allele count, all alleles, supertypes A, B, C, and SCL ("Carap_L_sl" from database download)
data_full <- data_full[c(1, 6, 7, 11, 14, 16, 19:139, 142)] 
colnames(data_full)

data_full <- data_full %>% mutate(location = tolower(location))

# Convert sample ID to character
data_full$sample_ID = as.character(data_full$sample_ID)

# Convert location, species, FP, paps_texture, season, and all alleles, all supertypes to factors
col_names <- names(data_full)
col_names

data_full[col_names[2:3]] <- lapply(data_full[col_names[2:3]], factor) # location and species

data_full[col_names[5:7]] <- lapply(data_full[col_names[5:7]], factor) # FP,paps_texture, season

data_full[col_names[9:127]] <- lapply(data_full[col_names[9:127]], factor) # alleles and supertypes

# Convert allele count, SCL, to numeric
data_full$allele_count <- as.numeric(data_full$allele_count)
data_full$Carap_L_sl <- as.numeric(data_full$Carap_L_sl)

# Just checking class and levels of each variable
sapply(data_full, class)
sapply(data_full, levels)

# Subset to the records of C. mydas with FP, so that just the C. mydas individuals that had tumors (which could either be active or regressing) are present in dataset

data <- data_full %>% filter(species == "Chelonia mydas") %>% filter(FP == 1)
data <- data[,-5] # remove "FP" variable which is part of a separate model ("RandomForest_FP_status.R")

# as a result of filtering to just FP-positive C. mydas, the following alleles are no longer relevant (they were found in just FP negative C. mydas)

data <- subset(data, select = -c(Caca01,Carca09,Carca32,Carca15,Carca19,Caca02,Carca14,Carca11,Carca16,Carca17,Carca27,Carca13,Chmy26,Chmy27,Carca56,Caca03,Carca28,Carca109,Caca04,Chmy59,Caca05,Carca25,Carca81,Caca06,Chmy77,Chmy78,Caca07,Chmy86,Caca08,Chmy88,Chmy92))

# Check for imbalance in response variable (regression)
length(which(data$paps_texture==0)) # 77 C. mydas with rough tumor texture
length(which(data$paps_texture==1)) # 28 C. mydas with smooth tumor texture for 26.7% regression.

# This is unbalanced enough that we need to account for it by over-sampling the under-represented class (regression) and under-sampling the over-represented class (active tumor). Set the sampsize parameter to 2/3 of the under-represented class, and use that sample size in conjunction with the strata option when building the random trees and when training.

# there are 28 regressed tumor individuals; 28*(2/3) = 18.67 == 19 individuals, so we'll sample 19 regressed tumor individuals and 19 active tumor individuals for the training data set so that the resutling tree forest is not biased.

sample_size <- c(19,19) # will use subsequently with strata.

# Random Forest analysis

# Optimize mtry parameter. Default for RF classification is sqrt(p) where p is the number of variables in the dataframe.
# let's also run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, in addition to p, where p is the number of variables.
# Run each mtry at ntree= 100 to 1000 (by increments of 100), looking for plateau where the out of bag error rate (OOB-ER) is minimized while also maximizing the RF algorithm. The mtry value that minimizes the OOB-ER will be chosen for subsequent analyses.

# We are using 95 variables to explain tumor regression status (this excludes the first column, which is sample ID, and regression status which is the response)

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3) # create matrix that the following loop will dump info into
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(9.7, 19.4, 9.4, 18.8, 31.3, 94)){    #values of mtry based on 96 total predictors
    rf_ij <- randomForest(x = data[,2:96], y = as.factor(data$paps_texture), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(data$paps_texture), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 9.7],results_optimization$OOB_ER[results_optimization$mtry == 9.7], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 19.4],results_optimization$OOB_ER[results_optimization$mtry == 19.4], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 9.4],results_optimization$OOB_ER[results_optimization$mtry == 9.4], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 18.8],results_optimization$OOB_ER[results_optimization$mtry == 18.8], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 31.3],results_optimization$OOB_ER[results_optimization$mtry == 31.3], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 94],results_optimization$OOB_ER[results_optimization$mtry == 94], col="red")


# all mtry parameters seem to behave similarly with no discernable plateau (although for mtry = 9.5, improvement after 200 trees), so it may be best to explore further or use defaults, i.e., sqrt(p)

# Create model with default paramters as a baseline, grid search with "oob" method
control <- trainControl(method="oob", number=10, search = "grid")
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(94)
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(paps_texture~., data=data[,2:96], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default) # mtry = 9.695 at 75.2381% accuracy.

# now let's fine tune
## Grid search, where method = "oob"
control <- trainControl(method="oob", number=10, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:94)) # doing full 94
rf_gridsearch <- train(paps_texture~., data=data[,2:96], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch) # most accurate: mtry 7 at 76.19048%

# Now begin the full Random Forest analyses: going to use mtry = 9.695 (default), and mtry = 7; 7  was given from grid search tuning with decently high accuracy.

# run a large number of trees with the above mtry values, in order to know what value of ntree achieves convergence of importance values between forests.
# As a starting point, grow 10,000 trees and increase if necessary.

# forests with default mtry:
rf_default_a <- randomForest(x = data[,2:96], y = as.factor(data$paps_texture), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$paps_texture), sampsize=sample_size)

rf_default_b <- randomForest(x = data[,2:96], y = as.factor(data$paps_texture), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$paps_texture), sampsize=sample_size)

#Check correlation of predictor importance values between forests 
importance_rf_default_a <- data.frame(importance(rf_default_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important predictor)
colnames(importance_rf_default_a)<-c("importance")

importance_rf_default_b <- data.frame(importance(rf_default_b,type=1))
colnames(importance_rf_default_b)<-c("importance")

cor(importance_rf_default_a,importance_rf_default_b) # A correlation of 0.9961736 for predictor importance values between forests when mtry = 9.695 and ntree = 10,000

# forests with mtry = 7

rf_7_a <- randomForest(x = data[,2:96], y = as.factor(data$paps_texture), importance=TRUE ,proximity=TRUE, mtry=7, ntree=10000, strata=as.factor(data$paps_texture), sampsize=sample_size)

rf_7_b <- randomForest(x = data[,2:96], y = as.factor(data$paps_texture), importance=TRUE ,proximity=TRUE, mtry=7, ntree=10000, strata=as.factor(data$paps_texture), sampsize=sample_size)

#Check correlation of locus importance values between forests 
importance_rf_7_a <- data.frame(importance(rf_7_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_7_a)<-c("importance")

importance_rf_7_b <- data.frame(importance(rf_7_b,type=1))
colnames(importance_rf_7_b)<-c("importance")

cor(importance_rf_7_a,importance_rf_7_b) # A correlation of 0.9961842 for predictor importance values between forests when mtry = 7 and ntree = 10,000

# Build final model with mtry = 7, ntree = 10,000, as mtry = 7 had slightly importance value correlation

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model = randomForest(paps_texture ~ .,
                          data = data.train[,-1],
                          importance=TRUE,
                          proximity=TRUE,
                          mtry=7, # based on tuning done above
                          ntree=10000, # based on tuning for ntree above
                          strata=as.factor(data$paps_texture), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model # OOB-ER = 40.62%

# Use the model built above to make predictions on the test data
data.predict = predict(data.model,newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict,Truth = data.test$paps_texture)
sum(diag(table(Prediction = data.predict,Truth = data.test$paps_texture))) / sum(table(Prediction = data.predict,Truth = data.test$paps_texture)) # Accuracy of model with all variables: #75.34247%

# Importance Sampling
importance(data.model)
varImpPlot(data.model, main = "Cm tumor smoothness")

# Let's also build a model based on the default mtry, since it did decently and is very close to the tuned mtry

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model_defaultmtry <- randomForest(paps_texture ~ .,
                                       data = data.train[,-1],
                                       importance=TRUE,
                                       proximity=TRUE,
                                       ntree=10000, # based on tuning for ntree above
                                       strata=as.factor(data$paps_texture), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model_defaultmtry # OOB-ER = 40.62%

# Use the model built above to make predictions on the test data
data.predict_defaultmtry = predict(data.model_defaultmtry, newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict_defaultmtry, Truth = data.test$paps_texture)
sum(diag(table(Prediction = data.predict_defaultmtry,Truth = data.test$paps_texture))) / sum(table(Prediction = data.predict_defaultmtry,Truth = data.test$paps_texture)) # Accuracy of model with all variables: #72.60%

# RF built with mtry = 7 has better accuracy than RF built with default mtry

# Importance Sampling
imp <- importance(data.model_defaultmtry)
varImpPlot(data.model_defaultmtry, main = "Cm tumor smoothness default mtry")