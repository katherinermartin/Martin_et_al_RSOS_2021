## Classification random forest for tumor regression status in C. mydas with FP
## for the sake of clarity with the use of the word "regression": here it refers to the phenomenon of fibropapillomas (tumors) in sea turtles regressing or lessening in severity over time. For this analysis and associated publication, it does not refer to the statistical analysis.

library(randomForest)
library(caret)
library(dplyr)

# read in data

data_full <- read.csv("classI_juveniles_morpho_FP.csv") # this includes C. mydas and C. caretta records

data_full <- data_full %>% filter(!is.na(Carap_L_sl)) # remove any record without straight carapace length (SCL) measurement (just 1)

# variables: sample ID, location, spp, year, FP, paps_regressed, season, allele count, all alleles, supertypes A, B, C, and SCL ("Carap_L_sl" from database download)
data_full <- data_full[c(1, 6, 7, 11, 14, 16, 19:139, 142)] 
colnames(data_full)

data_full <- data_full %>% mutate(location = tolower(location))

# Convert sample ID to character
data_full$sample_ID = as.character(data_full$sample_ID)

# Convert location, species, FP, paps_regressed, season, and all alleles, all supertypes to factors
col_names <- names(data_full)
col_names

data_full[col_names[2:3]] <- lapply(data_full[col_names[2:3]], factor) # location and species

data_full[col_names[5:7]] <- lapply(data_full[col_names[5:7]], factor) # FP,paps_regressed, season

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

data <- subset(data, select = -c(Caca01,Carca09,Carca32,Carca15,Carca19,Caca02,Carca14,Carca11,Carca16,Carca17,Carca27,Carca13,Chmy26,Chmy27,Carca56,Caca03,Carca28,Carca109,Caca04,Chmy59,Caca05,Carca25,Carca81,Caca06,Chmy77,Chmy78,Caca07,Chmy86,Caca08,Chmy88))

# Check for imbalance in response variable (regression)
length(which(data$paps_regressed==0)) # 77 C. mydas with no regression (aka active tumors)
length(which(data$paps_regressed==1)) # 28 C. mydas with regression, for 26.7% regression.

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
  for (j in c(10, 20, 9.5, 19, 32, 95)){    #values of mtry based on 96 total predictors
    rf_ij <- randomForest(x = data[,2:97], y = as.factor(data$paps_regressed), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(data$paps_regressed), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 10],results_optimization$OOB_ER[results_optimization$mtry == 10], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 20],results_optimization$OOB_ER[results_optimization$mtry == 20], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 9.5],results_optimization$OOB_ER[results_optimization$mtry == 9.5], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 19],results_optimization$OOB_ER[results_optimization$mtry == 19], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 32],results_optimization$OOB_ER[results_optimization$mtry == 32], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 95],results_optimization$OOB_ER[results_optimization$mtry == 95], col="red")


# all mtry parameters seem to behave similarly with no discernable plateau (although for mtry = 9.5, improvement after 200 trees), so it may be best to explore further or use defaults, i.e., sqrt(p)

# Create model with default paramters as a baseline, grid search with "oob" method
control <- trainControl(method="oob", number=10, search = "grid")
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(95)
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(paps_regressed~., data=data[,2:97], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default) # mtry = 9.746 at 72.238095% accuracy.

# now let's fine tune
## Grid search, where method = "oob"
control <- trainControl(method="oob", number=10, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:95)) # doing full 95
rf_gridsearch <- train(paps_regressed~., data=data[,2:97], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch) # most accurate: mtry 9 at 76.19048

# Now begin the full Random Forest analyses: going to use mtry = 9.746 (default), and mtry = 9; 9 was given from grid search tuning with decently high accuracy.

# run a large number of trees with the above mtry values, in order to know what value of ntree achieves convergence of importance values between forests.
# As a starting point, grow 10,000 trees and increase if necessary.

# forests with default mtry:
rf_default_a <- randomForest(x = data[,2:97], y = as.factor(data$paps_regressed), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$paps_regressed), sampsize=sample_size)
#save(rf_all_1,file="rf_all_1.Rdata")

rf_default_b <- randomForest(x = data[,2:97], y = as.factor(data$paps_regressed), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$paps_regressed), sampsize=sample_size)
#save(rf_all_2,file="rf_all_2.Rdata")

#Check correlation of predictor importance values between forests 
importance_rf_default_a <- data.frame(importance(rf_default_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important predictor)
colnames(importance_rf_default_a)<-c("importance")

importance_rf_default_b <- data.frame(importance(rf_default_b,type=1))
colnames(importance_rf_default_b)<-c("importance")

cor(importance_rf_default_a,importance_rf_default_b) # A correlation of 0.9965609 for predictor importance values between forests when mtry = 9.746 and ntree = 10,000

# forests with mtry = 9

rf_9_a <- randomForest(x = data[,2:97], y = as.factor(data$paps_regressed), importance=TRUE ,proximity=TRUE, mtry=9, ntree=10000, strata=as.factor(data$paps_regressed), sampsize=sample_size)
#save(rf_all_1,file="rf_all_1.Rdata")

rf_9_b <- randomForest(x = data[,2:97], y = as.factor(data$paps_regressed), importance=TRUE ,proximity=TRUE, mtry=9, ntree=10000, strata=as.factor(data$paps_regressed), sampsize=sample_size)
#save(rf_all_2,file="rf_all_2.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_9_a <- data.frame(importance(rf_9_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_9_a)<-c("importance")

importance_rf_9_b <- data.frame(importance(rf_9_b,type=1))
colnames(importance_rf_9_b)<-c("importance")

cor(importance_rf_9_a,importance_rf_9_b) # A correlation of 0.997668 for predictor importance values between forests when mtry = 9 and ntree = 10,000

# Build final model with mtry = 9, ntree = 10,000, as mtry = 9 had slightly importance value correlation

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model = randomForest(paps_regressed ~ .,
                          data = data.train[,-1],
                          importance=TRUE,
                          proximity=TRUE,
                          mtry=9, # based on tuning done above
                          ntree=10000, # based on tuning for ntree above
                          strata=as.factor(data$paps_regressed), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model # OOB-ER = 40.62%

# Use the model built above to make predictions on the test data
data.predict = predict(data.model,newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict,Truth = data.test$paps_regressed)
sum(diag(table(Prediction = data.predict,Truth = data.test$paps_regressed))) / sum(table(Prediction = data.predict,Truth = data.test$paps_regressed)) # Accuracy of model with all variables: #73.97%

# Importance Sampling
importance(data.model)
varImpPlot(data.model, main = "Cm tumor regression")
dev.copy(pdf, file="~/Documents/UCF/Research/MHC_Class_I/Analysis/Figures/RF/Feb21/RF_importance_tumor_regression.pdf")
dev.off()
