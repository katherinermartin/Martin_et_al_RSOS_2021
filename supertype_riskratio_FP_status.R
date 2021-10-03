# Risk ratio analysis: FP risk by supertype, in C. mydas

library(epiR)
library(tidyverse)
library(ggplot2)
library(svglite)


# Load data
data <- read.csv("/Users/KatieMartin/Documents/UCF/Research/MHC_Class_I/Data/supertype_dataframe_assignment/classI_juveniles_morpho_FP_v3.csv")


# filter to C. mydas records only
Cm_data <- data %>% filter(species == "Chelonia mydas")

Cm_data[Cm_data$FP == 1, "FP"] <- "FP+" # replace 1s with yes and 0s with nos for FP
Cm_data[Cm_data$FP == 0, "FP"] <- "FP-"


# st_A
Cm_data[Cm_data$supertype_A == 1, "supertype_A"] <- "st_A +"
Cm_data[Cm_data$supertype_A == 0, "supertype_A"] <- "st_A -"

tab_st_A <- table(Cm_data$supertype_A, Cm_data$FP) # independent, dependent
tab_st_A

# rearrange so it's in ABCD format
tab2_st_A <- cbind(tab_st_A[,2], tab_st_A[,1]) # rearrange columns
colnames(tab2_st_A) <- c("FP+", "FP-")

tab2_st_A <- rbind(tab2_st_A[2,], tab2_st_A[1,]) # rearrange rows
rownames(tab2_st_A) <- c("st_A +", "st_A -")

tab2_st_A

barplot(tab2_st_A, beside = TRUE, main = "FP in turtles with and without supertype A alleles", legend.text = c("st_A +", "st_A -"))

st_A_risk <- epi.2by2(tab2_st_A, method = "cohort.count", conf.level = 0.95)

st_A_risk

st_A_risk <- as.data.frame(st_A_risk$massoc.summary)

# st_B
Cm_data[Cm_data$supertype_B == 1, "supertype_B"] <- "st_B +"
Cm_data[Cm_data$supertype_B == 0, "supertype_B"] <- "st_B -"

tab_st_B <- table(Cm_data$supertype_B, Cm_data$FP) # independent, dependent
tab_st_B

# rearrange so it's in ABCD format
tab2_st_B <- cbind(tab_st_B[,2], tab_st_B[,1]) # rearrange columns
colnames(tab2_st_B) <- c("FP+", "FP-")

tab2_st_B <- rbind(tab2_st_B[2,], tab2_st_B[1,]) # rearrange rows
rownames(tab2_st_B) <- c("st_B +", "st_B -")

tab2_st_B

barplot(tab2_st_B, beside = TRUE, main = "FP in turtles with and without supertype B alleles", legend.text = c("st_B +", "st_B -"))

st_B_risk <- epi.2by2(tab2_st_B, method = "cohort.count", conf.level = 0.95)

st_B_risk

st_B_risk <- as.data.frame(st_B_risk$massoc.summary)

supertype_FP_risk <- rbind(st_A_risk, st_B_risk)

# st_C
Cm_data[Cm_data$supertype_C == 1, "supertype_C"] <- "st_C +"
Cm_data[Cm_data$supertype_C == 0, "supertype_C"] <- "st_C -"

tab_st_C <- table(Cm_data$supertype_C, Cm_data$FP) # independent, dependent
tab_st_C

# rearrange so it's in ABCD format
tab2_st_C <- cbind(tab_st_C[,2], tab_st_C[,1]) # rearrange columns
colnames(tab2_st_C) <- c("FP+", "FP-")

tab2_st_C <- rbind(tab2_st_C[2,], tab2_st_C[1,]) # rearrange rows
rownames(tab2_st_C) <- c("st_C +", "st_C -")

tab2_st_C

barplot(tab2_st_C, beside = TRUE, main = "FP in turtles with and without supertype C alleles", legend.text = c("st_C +", "st_C -"))

st_C_risk <- epi.2by2(tab2_st_C, method = "cohort.count", conf.level = 0.95)

st_C_risk

st_C_risk <- as.data.frame(st_C_risk$massoc.summary)

supertype_FP_risk <- rbind(supertype_FP_risk, st_C_risk)

supertype_FP_risk <- supertype_FP_risk %>% filter(var == "Inc risk ratio") # filter  to just "inc risk ratio"

supertype_FP_risk$supertype <- c("A", "B", "C")

# graph a forest plot

supertype_FP <- ggplot(data=supertype_FP_risk, aes(x=supertype, y=est, ymin=lower, ymax=upper)) +
  geom_pointrange(shape = 18) +
  geom_hline(yintercept=1, lty=2, color = "red") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("supertype") + ylab("relative risk estimate of FP (95% CI)") +
  theme_bw() # use a white background
supertype_FP

# Need Bonferroni correction for comparing 3 classes

# CI will be 98.3%; familywise error rate is 95% so 1-(alpha/m) where m = 3 comparisons so 1-(0.05/6) = 0.983

library(epiR)
library(tidyverse)
library(ggplot2)
library(svglite)


# Load data
data <- read.csv("/Users/KatieMartin/Documents/UCF/Research/MHC_Class_I/Data/supertype_dataframe_assignment/classI_juveniles_morpho_FP_v3.csv")


# filter to C. mydas records only
Cm_data <- data %>% filter(species == "Chelonia mydas")

Cm_data[Cm_data$FP == 1, "FP"] <- "FP+" # replace 1s with yes and 0s with nos for FP
Cm_data[Cm_data$FP == 0, "FP"] <- "FP-"


# st_A
Cm_data[Cm_data$supertype_A == 1, "supertype_A"] <- "st_A +"
Cm_data[Cm_data$supertype_A == 0, "supertype_A"] <- "st_A -"

tab_st_A <- table(Cm_data$supertype_A, Cm_data$FP) # independent, dependent
tab_st_A

# rearrange so it's in ABCD format
tab2_st_A <- cbind(tab_st_A[,2], tab_st_A[,1]) # rearrange columns
colnames(tab2_st_A) <- c("FP+", "FP-")

tab2_st_A <- rbind(tab2_st_A[2,], tab2_st_A[1,]) # rearrange rows
rownames(tab2_st_A) <- c("st_A +", "st_A -")

tab2_st_A

barplot(tab2_st_A, beside = TRUE, main = "FP in turtles with and without supertype A alleles", legend.text = c("st_A +", "st_A -"))

st_A_risk <- epi.2by2(tab2_st_A, method = "cohort.count", conf.level = 0.983)

st_A_risk

st_A_risk <- as.data.frame(st_A_risk$massoc.summary)

# st_B
Cm_data[Cm_data$supertype_B == 1, "supertype_B"] <- "st_B +"
Cm_data[Cm_data$supertype_B == 0, "supertype_B"] <- "st_B -"

tab_st_B <- table(Cm_data$supertype_B, Cm_data$FP) # independent, dependent
tab_st_B

# rearrange so it's in ABCD format
tab2_st_B <- cbind(tab_st_B[,2], tab_st_B[,1]) # rearrange columns
colnames(tab2_st_B) <- c("FP+", "FP-")

tab2_st_B <- rbind(tab2_st_B[2,], tab2_st_B[1,]) # rearrange rows
rownames(tab2_st_B) <- c("st_B +", "st_B -")

tab2_st_B

barplot(tab2_st_B, beside = TRUE, main = "FP in turtles with and without supertype B alleles", legend.text = c("st_B +", "st_B -"))

st_B_risk <- epi.2by2(tab2_st_B, method = "cohort.count", conf.level = 0.983)

st_B_risk

st_B_risk <- as.data.frame(st_B_risk$massoc.summary)

supertype_FP_risk <- rbind(st_A_risk, st_B_risk)

# st_C
Cm_data[Cm_data$supertype_C == 1, "supertype_C"] <- "st_C +"
Cm_data[Cm_data$supertype_C == 0, "supertype_C"] <- "st_C -"

tab_st_C <- table(Cm_data$supertype_C, Cm_data$FP) # independent, dependent
tab_st_C

# rearrange so it's in ABCD format
tab2_st_C <- cbind(tab_st_C[,2], tab_st_C[,1]) # rearrange columns
colnames(tab2_st_C) <- c("FP+", "FP-")

tab2_st_C <- rbind(tab2_st_C[2,], tab2_st_C[1,]) # rearrange rows
rownames(tab2_st_C) <- c("st_C +", "st_C -")

tab2_st_C

barplot(tab2_st_C, beside = TRUE, main = "FP in turtles with and without supertype C alleles", legend.text = c("st_C +", "st_C -"))

st_C_risk <- epi.2by2(tab2_st_C, method = "cohort.count", conf.level = 0.983)

st_C_risk

st_C_risk <- as.data.frame(st_C_risk$massoc.summary)

supertype_FP_risk <- rbind(supertype_FP_risk, st_C_risk)

supertype_FP_risk <- supertype_FP_risk %>% filter(var == "Inc risk ratio") # filter  to just "inc risk ratio"

supertype_FP_risk$supertype <- c("A", "B", "C")

# graph a forest plot

supertype_FP <- ggplot(data=supertype_FP_risk, aes(x=supertype, y=est, ymin=lower, ymax=upper)) +
  geom_pointrange(shape = 18) +
  geom_hline(yintercept=1, lty=2, color = "red") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("supertype") + ylab("relative risk estimate of FP (98.3% CI)") +
  theme_bw() # use a white background
supertype_FP
