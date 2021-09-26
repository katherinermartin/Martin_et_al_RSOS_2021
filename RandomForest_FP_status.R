## Classification random forest for FP occurence in C. mydas

library(randomForest)
library(caret)
library(dplyr)

# read in data

data_full <- read.csv("classI_juveniles_morpho_FP_v3.csv") # this includes C. mydas and C. caretta records

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

data <- subset(data, select = -c(Carca32, Carca14, Carca11, Carca16, Carca17, Carca27, Carca13, Carca56, Caca03, Carca28, Carca109, Caca04, Caca05, Carca25, Caca06, Caca07, Caca08, Carca81))

# remove paps_regressed, which is part of a separate model ("RandomForest_TumorRegression.R").
data <- data[,-6]


# Check for imbalance in response variable (FP)
length(which(data$FP==0)) #163 Cm without FP
length(which(data$FP==1)) #105 Cm with FP, for ~39.17% FP. This is unbalanced enough that we need to account for it by over-sampling the under-represented class (FP positive) and under-sampling the over-represented class (FP negative). Set the sampsize parameter to 2/3 of the under-represented class, and use that sample size in conjunction with the strata option when building the random trees and when training.

# there are 105 FP positive individuals; 105*(2/3) = 70 individuals, so we'll sample 70 FP positive individuals and 70 FP negative individuals for the training data set so that the resutling tree forest is not biased.

sample_size <- c(70,70) # will use subsequently with strata.

# Random Forest analysis

# Optimize mtry parameter. Default for RF classification is sqrt(p) where p is the number of variables in the dataframe.
# let's also run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, in addition to p, where p is the number of variables.
# Run each mtry at ntree= 100 to 1000 (by increments of 100), looking for plateau where the out of bag error rate (OOB-ER) is minimized while also maximizing the RF algorithm. The mtry value that minimizes the OOB-ER will be chosen for subsequent analyses.

# We are using 107 variables to explain FP occurence (this excludes the first column, which is sample ID, and FP occurence which is the response)

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3) # create matrix that the following loop will dump info into
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(10.3, 20.7, 10.7, 21.4, 35.7, 107)){    #values of mtry based on 107 total predictors
    rf_ij <- randomForest(x = data[,2:109], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(data$FP), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 10.3],results_optimization$OOB_ER[results_optimization$mtry == 10.3], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 20.7],results_optimization$OOB_ER[results_optimization$mtry == 20.7], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 10.7],results_optimization$OOB_ER[results_optimization$mtry == 10.7], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 21.4],results_optimization$OOB_ER[results_optimization$mtry == 21.4], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 35.7],results_optimization$OOB_ER[results_optimization$mtry == 35.7], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 107],results_optimization$OOB_ER[results_optimization$mtry == 107], col="red")

# all mtry parameters seem to behave similarly with no discernable plateau, so it may be best to explore further or use defaults, i.e., sqrt(p)

# Create model with default paramters as a baseline with grid search, method = "oob"
control <- trainControl(method="oob", number=10, search = "grid")
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(107) # square root of total number of variables; this is the default mtry
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(FP~., data=data[,2:109], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default) # mtry = 10.34 at 73.5% accuracy.

# now let's fine tune
## Grid search, where method = "oob"
control <- trainControl(method="oob", number=10, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:107)) # doing full 107 mtry
rf_gridsearch <- train(FP~., data=data[,2:109], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch) # most accurate: mtry 21 at 74.6%

# Now begin the full Random Forest analyses: going to use mtry = 10.34 (default), and mtry = 21 which was given from grid search tuning.

# run a large number of trees with the above mtry values, in order to know what value of ntree achieves convergence of importance values between forests.
# As a starting point, grow 10,000 trees and increase if necessary.

# forests with default mtry:

rf_default_a <- randomForest(x = data[,2:109], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_default_a,file="rf_default_a.Rdata")

rf_default_b <- randomForest(x = data[,2:109], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_default_b,file="rf_default_b.Rdata")

#Check correlation of predictor importance values between forests 
importance_rf_default_a <- data.frame(importance(rf_default_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important predictor)
colnames(importance_rf_default_a)<-c("importance")

importance_rf_default_b <- data.frame(importance(rf_default_b,type=1))
colnames(importance_rf_default_b)<-c("importance")

cor(importance_rf_default_a,importance_rf_default_b) # A correlation of 0.998319 for predictor importance values between forests when mtry = 10.34 and ntree = 10,000

# forest with mtry = 21

rf_21_a <- randomForest(x = data[,2:109], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, mtry=21, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_21_a,file="rf_21_a.Rdata")

rf_21_b <- randomForest(x = data[,2:109], y = as.factor(data$FP), importance=TRUE ,proximity=TRUE, mtry=21, ntree=10000, strata=as.factor(data$FP), sampsize=sample_size)
#save(rf_21_b,file="rf_21_b.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_21_a <- data.frame(importance(rf_21_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_21_a)<-c("importance")

importance_rf_21_b <- data.frame(importance(rf_21_b,type=1))
colnames(importance_rf_21_b)<-c("importance")

cor(importance_rf_21_a,importance_rf_21_b) # A correlation of 0.9994113 for predictor importance values between forests when mtry = 21 and ntree = 10,000

# Build final model with mtry = 21, ntree = 10,000, as mtry = 21 had higher importance value correlation

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
                          mtry=21, # based on tuning done above
                          ntree=10000, # based on tuning for ntree above
                          strata=as.factor(data$FP), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model # OOB-ER = 29.63%

# Use the model built above to make predictions on the test data
data.predict = predict(data.model,newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict,Truth = data.test$FP)
sum(diag(table(Prediction = data.predict,Truth = data.test$FP))) / sum(table(Prediction = data.predict,Truth = data.test$FP)) # Accuracy of model with all variables: #70.05%

# Importance Sampling
importance(data.model)
varImpPlot(data.model, main = "Cm FP")



# Let's also build a model based on the default mtry, since it also did decently in the tuning.

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model_10 <- randomForest(FP ~ .,
                              data = data.train[,-1],
                              importance=TRUE,
                              proximity=TRUE,
                              ntree=10000, # based on tuning for ntree above
                              strata=as.factor(data$FP), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model_10 # OOB-ER = 37.04%

# Use the model built above to make predictions on the test data
data.predict_10 = predict(data.model_10, newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict_10, Truth = data.test$FP)
sum(diag(table(Prediction = data.predict_10,Truth = data.test$FP))) / sum(table(Prediction = data.predict_10,Truth = data.test$FP)) # Accuracy of model with all variables: #67.914%
# accuracy is comparable between the two mtry values, but the OOB-ER is better with the tuned mtry value