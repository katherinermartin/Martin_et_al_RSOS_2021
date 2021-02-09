## Classification random forest for FP occurence in C. mydas

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

sapply(data_full, class)
sapply(data_full, levels)

## Subset to just FP with just C. mydas.

data <- data_full %>% filter(species == "Chelonia mydas") # keep just the Cm; 268 individuals

# as a result of filtering to just Cm individuals, the following alleles are no longer relevant (they were found just in Cc) and can be removed

data <- subset(data, select = -c(Carca32, Carca14, Carca11, Carca16, Carca17, Carca27, Carca13, Carca56, Caca03, Carca28, Carca109, Caca04, Caca05, Carca25, Caca06, Caca07, Caca08))

# remove paps_regressed, which is part of a separate model ("RandomForest_TumorRegression.R").
data <- data[,-6]


# Check for imbalance in response variable (FP)
length(which(data$FP==0)) #163 Cm without FP
length(which(data$FP==1)) #105 Cm with FP, for ~30.7% FP. This is unbalanced enough that we need to account for it by over-sampling the under-represented class (FP positive) and under-sampling the over-represented class (FP negative). Set the sampsize parameter to 2/3 of the under-represented class, and use that sample size in conjunction with the strata option when building the random trees and when training.

# there are 105 FP positive individuals; 105*(2/3) = 70 individuals, so we'll sample 70 FP positive individuals and 70 FP negative individuals for the training data set so that the resutling tree forest is not biased.

sample_size <- c(70,70) # will use subsequently with strata.

# Random Forest analysis

# Optimize mtry parameter. Default for RF classification is sqrt(p) where p is the number of variables in the dataframe.
# let's also run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, in addition to p, where p is the number of variables.
# Run each mtry at ntree= 100 to 1000 (by increments of 100), looking for plateau where the out of bag error rate (OOB-ER) is minimized while also maximizing the RF algorithm. The mtry value that minimizes the OOB-ER will be chosen for subsequent analyses.

# We are using 108 variables to explain FP occurence (this excludes the first column, which is sample ID, and FP occurence which is the response)

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3) # create matrix that the following loop will dump info into
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(10, 21, 11, 21.6, 36, 108)){    #values of mtry based on 109 total predictors
    rf_ij <- randomForest(x = data[,2:110], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(data$FP), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 10],results_optimization$OOB_ER[results_optimization$mtry == 10], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 21],results_optimization$OOB_ER[results_optimization$mtry == 21], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 11],results_optimization$OOB_ER[results_optimization$mtry == 11], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 21.6],results_optimization$OOB_ER[results_optimization$mtry == 21.6], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 36],results_optimization$OOB_ER[results_optimization$mtry == 36], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 108],results_optimization$OOB_ER[results_optimization$mtry == 108], col="red")

# all mtry parameters seem to behave similarly with no discernable plateau, so it may be best to explore further or use defaults, i.e., sqrt(p)

# Create model with default paramters as a baseline with grid search, method = "oob"
control <- trainControl(method="oob", number=10, search = "grid")
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(108) # square root of total number of variables; this is the default mtry
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(FP~., data=data[,2:110], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default) # mtry = 10.39 at 71.6% accuracy.

# now let's fine tune
## Grid search, where method = "oob"
control <- trainControl(method="oob", number=10, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:108)) # doing full 108
rf_gridsearch <- train(FP~., data=data[,2:110], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch) # most accurate: mtry 82 at 74.6%

# Now begin the full Random Forest analyses: going to use mtry = 10.39 (default), and mtry = 82, which was given from grid search tuning.

# run a large number of trees with the above mtry values, in order to know what value of ntree achieves convergence of importance values between forests.
# As a starting point, grow 10,000 trees and increase if necessary.

# forests with default mtry:

rf_default_a <- randomForest(x = data[,2:110], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_all_1,file="rf_all_1.Rdata")

rf_default_b <- randomForest(x = data[,2:110], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_all_2,file="rf_all_2.Rdata")

#Check correlation of predictor importance values between forests 
importance_rf_default_a <- data.frame(importance(rf_default_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important predictor)
colnames(importance_rf_default_a)<-c("importance")

importance_rf_default_b <- data.frame(importance(rf_default_b,type=1))
colnames(importance_rf_default_b)<-c("importance")

cor(importance_rf_default_a,importance_rf_default_b) # A correlation of 0.9985 for predictor importance values between forests when mtry = 10 and ntree = 10,000

# forest with mtry = 82

rf_82_a <- randomForest(x = data[,2:110], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, mtry=82, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_all_1,file="rf_all_1.Rdata")

rf_82_b <- randomForest(x = data[,2:110], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, mtry=82, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_all_2,file="rf_all_2.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_82_a <- data.frame(importance(rf_82_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_82_a)<-c("importance")

importance_rf_82_b <- data.frame(importance(rf_82_b,type=1))
colnames(importance_rf_82_b)<-c("importance")

cor(importance_rf_82_a,importance_rf_82_b) # A correlation of 0.9999 for predictor importance values between forests when mtry = 82 and ntree = 10,000

# Build final model with mtry = 82, ntree = 10,000, as mtry = 82 had higher importance value correlation

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model = randomForest(FP ~ .,
                          data = data.train[,-1],
                          importance=TRUE,
                          proximity=TRUE,
                          mtry=82, # based on tuning done above
                          ntree=10000, # based on tuning for ntree above
                          strata=as.factor(data$FP), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model # OOB-ER = 25.93%

# Use the model built above to make predictions on the test data
data.predict = predict(data.model,newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict,Truth = data.test$FP)
sum(diag(table(Prediction = data.predict,Truth = data.test$FP))) / sum(table(Prediction = data.predict,Truth = data.test$FP)) # Accuracy of model with all variables: #67.37%

# Importance Sampling
importance(data.model)
varImpPlot(data.model, main = "Cm FP")
dev.copy(pdf, file="~/Documents/UCF/Research/MHC_Class_I/Analysis/Figures/RF/Feb21/RF_importance_FP.pdf")
dev.off()