# Risk ratio of top alleles in C. mydas with FP, with Bonferroni correction where the confidence level for the calculation is set to 99.8%

library(epiR)
library(tidyverse)
library(ggplot2)
library(svglite)

# Load data
data <- read.csv("classI_juveniles_morpho_FP.csv")

# filter to C. mydas records only
Cm_data <- data %>% filter(species == "Chelonia mydas")

# Subset to alleles that occur in 10 or more C. mydas
Cm_data <- subset(Cm_data, select = c(sample_ID, FP, Chmy20, Chmy35, Chmy62, Chmy14, Chmy50, Chmy08, Chmy13, Chmy07, Chmy18, Chmy33, Chmy49, Chmy17, Chmy48, Chmy12, Chmy21, Chmy11, Chmy06, Chmy41, Chmy09, Chmy03, Chmy05, Chmy04, Chmy10, Chmy01, Chmy02))


Cm_data[Cm_data$FP == 1, "FP"] <- "FP+" # replace 1s with yes and 0s with nos for FP
Cm_data[Cm_data$FP == 0, "FP"] <- "FP-"

attach(Cm_data)

# Chmy02
Cm_data[Cm_data$Chmy02 == 1, "Chmy02"] <- "Chmy02 +"
Cm_data[Cm_data$Chmy02 == 0, "Chmy02"] <- "Chmy02 -"

tab_Chmy02 <- table(Chmy02, FP) # independent, dependent
tab_Chmy02

# rearrange so it's in ABCD format
tab2_Chmy02 <- cbind(tab_Chmy02[,2], tab_Chmy02[,1]) # rearrange columns
colnames(tab2_Chmy02) <- c("FP+", "FP-")

tab2_Chmy02 <- rbind(tab2_Chmy02[2,], tab2_Chmy02[1,]) # rearrange rows
rownames(tab2_Chmy02) <- c("Chmy02 +", "Chmy02 -")

tab2_Chmy02

barplot(tab2_Chmy02, beside = TRUE, main = "FP in turtles with and without Chmy02 allele", legend.text = c("Chmy02 +", "Chmy02 -"))

Chmy02_risk <- epi.2by2(tab2_Chmy02, method = "cohort.count", conf.level = 0.998)

Chmy02_risk

df_Chmy02 <- as.data.frame(Chmy02_risk$massoc$RR.strata.wald)

# Chmy01
Cm_data[Cm_data$Chmy01 == 1, "Chmy01"] <- "Chmy01 +"
Cm_data[Cm_data$Chmy01 == 0, "Chmy01"] <- "Chmy01 -"

tab_Chmy01 <- table(Chmy01, FP) # independent, dependent
tab_Chmy01

# rearrange so it's in ABCD format
tab2_Chmy01 <- cbind(tab_Chmy01[,2], tab_Chmy01[,1]) # rearrange columns
colnames(tab2_Chmy01) <- c("FP+", "FP-")

tab2_Chmy01 <- rbind(tab2_Chmy01[2,], tab2_Chmy01[1,]) # rearrange rows
rownames(tab2_Chmy01) <- c("Chmy01 +", "Chmy01 -")

tab2_Chmy01

barplot(tab2_Chmy01, beside = TRUE, main = "FP in turtles with and without Chmy01 allele", legend.text = c("Chmy01 +", "Chmy01 -"))

Chmy01_risk <- epi.2by2(tab2_Chmy01, method = "cohort.count", conf.level = 0.998)

Chmy01_risk

df_Chmy01 <- as.data.frame(Chmy01_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(df_Chmy02, df_Chmy01)

# Chmy10
Cm_data[Cm_data$Chmy10 == 1, "Chmy10"] <- "Chmy10 +"
Cm_data[Cm_data$Chmy10 == 0, "Chmy10"] <- "Chmy10 -"

tab_Chmy10 <- table(Chmy10, FP) # independent, dependent
tab_Chmy10

# rearrange so it's in ABCD format
tab2_Chmy10 <- cbind(tab_Chmy10[,2], tab_Chmy10[,1]) # rearrange columns
colnames(tab2_Chmy10) <- c("FP+", "FP-")

tab2_Chmy10 <- rbind(tab2_Chmy10[2,], tab2_Chmy10[1,]) # rearrange rows
rownames(tab2_Chmy10) <- c("Chmy10 +", "Chmy10 -")

tab2_Chmy10

barplot(tab2_Chmy10, beside = TRUE, main = "FP in turtles with and without Chmy10 allele", legend.text = c("Chmy10 +", "Chmy10 -"))

Chmy10_risk <- epi.2by2(tab2_Chmy10, method = "cohort.count", conf.level = 0.998)

Chmy10_risk

df_Chmy10 <- as.data.frame(Chmy10_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy10)

# Chmy04
Cm_data[Cm_data$Chmy04 == 1, "Chmy04"] <- "Chmy04 +"
Cm_data[Cm_data$Chmy04 == 0, "Chmy04"] <- "Chmy04 -"

tab_Chmy04 <- table(Chmy04, FP) # independent, dependent
tab_Chmy04

# rearrange so it's in ABCD format
tab2_Chmy04 <- cbind(tab_Chmy04[,2], tab_Chmy04[,1]) # rearrange columns
colnames(tab2_Chmy04) <- c("FP+", "FP-")

tab2_Chmy04 <- rbind(tab2_Chmy04[2,], tab2_Chmy04[1,]) # rearrange rows
rownames(tab2_Chmy04) <- c("Chmy04 +", "Chmy04 -")

tab2_Chmy04

barplot(tab2_Chmy04, beside = TRUE, main = "FP in turtles with and without Chmy04 allele", legend.text = c("Chmy04 +", "Chmy04 -"))

Chmy04_risk <- epi.2by2(tab2_Chmy04, method = "cohort.count", conf.level = 0.998)

Chmy04_risk

df_Chmy04 <- as.data.frame(Chmy04_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy04)

# Chmy05
Cm_data[Cm_data$Chmy05 == 1, "Chmy05"] <- "Chmy05 +"
Cm_data[Cm_data$Chmy05 == 0, "Chmy05"] <- "Chmy05 -"

tab_Chmy05 <- table(Chmy05, FP) # independent, dependent
tab_Chmy05

# rearrange so it's in ABCD format
tab2_Chmy05 <- cbind(tab_Chmy05[,2], tab_Chmy05[,1]) # rearrange columns
colnames(tab2_Chmy05) <- c("FP+", "FP-")

tab2_Chmy05 <- rbind(tab2_Chmy05[2,], tab2_Chmy05[1,]) # rearrange rows
rownames(tab2_Chmy05) <- c("Chmy05 +", "Chmy05 -")

tab2_Chmy05

barplot(tab2_Chmy05, beside = TRUE, main = "FP in turtles with and without Chmy05 allele", legend.text = c("Chmy05 +", "Chmy05 -"))

Chmy05_risk <- epi.2by2(tab2_Chmy05, method = "cohort.count", conf.level = 0.998)

Chmy05_risk

df_Chmy05 <- as.data.frame(Chmy05_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy05)

# Chmy03
Cm_data[Cm_data$Chmy03 == 1, "Chmy03"] <- "Chmy03 +"
Cm_data[Cm_data$Chmy03 == 0, "Chmy03"] <- "Chmy03 -"

tab_Chmy03 <- table(Chmy03, FP) # independent, dependent
tab_Chmy03

# rearrange so it's in ABCD format
tab2_Chmy03 <- cbind(tab_Chmy03[,2], tab_Chmy03[,1]) # rearrange columns
colnames(tab2_Chmy03) <- c("FP+", "FP-")

tab2_Chmy03 <- rbind(tab2_Chmy03[2,], tab2_Chmy03[1,]) # rearrange rows
rownames(tab2_Chmy03) <- c("Chmy03 +", "Chmy03 -")

tab2_Chmy03

barplot(tab2_Chmy03, beside = TRUE, main = "FP in turtles with and without Chmy03 allele", legend.text = c("Chmy03 +", "Chmy03 -"))

Chmy03_risk <- epi.2by2(tab2_Chmy03, method = "cohort.count", conf.level = 0.998)

Chmy03_risk

df_Chmy03 <- as.data.frame(Chmy03_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy03)

# Chmy09
Cm_data[Cm_data$Chmy09 == 1, "Chmy09"] <- "Chmy09 +"
Cm_data[Cm_data$Chmy09 == 0, "Chmy09"] <- "Chmy09 -"

tab_Chmy09 <- table(Chmy09, FP) # independent, dependent
tab_Chmy09

# rearrange so it's in ABCD format
tab2_Chmy09 <- cbind(tab_Chmy09[,2], tab_Chmy09[,1]) # rearrange columns
colnames(tab2_Chmy09) <- c("FP+", "FP-")

tab2_Chmy09 <- rbind(tab2_Chmy09[2,], tab2_Chmy09[1,]) # rearrange rows
rownames(tab2_Chmy09) <- c("Chmy09 +", "Chmy09 -")

tab2_Chmy09

barplot(tab2_Chmy09, beside = TRUE, main = "FP in turtles with and without Chmy09 allele", legend.text = c("Chmy09 +", "Chmy09 -"))

Chmy09_risk <- epi.2by2(tab2_Chmy09, method = "cohort.count", conf.level = 0.998)

Chmy09_risk

df_Chmy09 <- as.data.frame(Chmy09_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy09)

# Chmy06
Cm_data[Cm_data$Chmy06 == 1, "Chmy06"] <- "Chmy06 +"
Cm_data[Cm_data$Chmy06 == 0, "Chmy06"] <- "Chmy06 -"

tab_Chmy06 <- table(Chmy06, FP) # independent, dependent
tab_Chmy06

# rearrange so it's in ABCD format
tab2_Chmy06 <- cbind(tab_Chmy06[,2], tab_Chmy06[,1]) # rearrange columns
colnames(tab2_Chmy06) <- c("FP+", "FP-")

tab2_Chmy06 <- rbind(tab2_Chmy06[2,], tab2_Chmy06[1,]) # rearrange rows
rownames(tab2_Chmy06) <- c("Chmy06 +", "Chmy06 -")

tab2_Chmy06

barplot(tab2_Chmy06, beside = TRUE, main = "FP in turtles with and without Chmy06 allele", legend.text = c("Chmy06 +", "Chmy06 -"))

Chmy06_risk <- epi.2by2(tab2_Chmy06, method = "cohort.count", conf.level = 0.998)

Chmy06_risk

df_Chmy06 <- as.data.frame(Chmy06_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy06)

# Chmy41
Cm_data[Cm_data$Chmy41 == 1, "Chmy41"] <- "Chmy41 +"
Cm_data[Cm_data$Chmy41 == 0, "Chmy41"] <- "Chmy41 -"

tab_Chmy41 <- table(Chmy41, FP) # independent, dependent
tab_Chmy41

# rearrange so it's in ABCD format
tab2_Chmy41 <- cbind(tab_Chmy41[,2], tab_Chmy41[,1]) # rearrange columns
colnames(tab2_Chmy41) <- c("FP+", "FP-")

tab2_Chmy41 <- rbind(tab2_Chmy41[2,], tab2_Chmy41[1,]) # rearrange rows
rownames(tab2_Chmy41) <- c("Chmy41 +", "Chmy41 -")

tab2_Chmy41

barplot(tab2_Chmy41, beside = TRUE, main = "FP in turtles with and without Chmy41 allele", legend.text = c("Chmy41 +", "Chmy41 -"))

Chmy41_risk <- epi.2by2(tab2_Chmy41, method = "cohort.count", conf.level = 0.998)

Chmy41_risk

df_Chmy41 <- as.data.frame(Chmy41_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy41)

# Chmy11
Cm_data[Cm_data$Chmy11 == 1, "Chmy11"] <- "Chmy11 +"
Cm_data[Cm_data$Chmy11 == 0, "Chmy11"] <- "Chmy11 -"

tab_Chmy11 <- table(Chmy11, FP) # independent, dependent
tab_Chmy11

# rearrange so it's in ABCD format
tab2_Chmy11 <- cbind(tab_Chmy11[,2], tab_Chmy11[,1]) # rearrange columns
colnames(tab2_Chmy11) <- c("FP+", "FP-")

tab2_Chmy11 <- rbind(tab2_Chmy11[2,], tab2_Chmy11[1,]) # rearrange rows
rownames(tab2_Chmy11) <- c("Chmy11 +", "Chmy11 -")

tab2_Chmy11

barplot(tab2_Chmy11, beside = TRUE, main = "FP in turtles with and without Chmy11 allele", legend.text = c("Chmy11 +", "Chmy11 -"))

Chmy11_risk <- epi.2by2(tab2_Chmy11, method = "cohort.count", conf.level = 0.998)

Chmy11_risk

df_Chmy11 <- as.data.frame(Chmy11_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy11)

# Chmy12
Cm_data[Cm_data$Chmy12 == 1, "Chmy12"] <- "Chmy12 +"
Cm_data[Cm_data$Chmy12 == 0, "Chmy12"] <- "Chmy12 -"

tab_Chmy12 <- table(Chmy12, FP) # independent, dependent
tab_Chmy12

# rearrange so it's in ABCD format
tab2_Chmy12 <- cbind(tab_Chmy12[,2], tab_Chmy12[,1]) # rearrange columns
colnames(tab2_Chmy12) <- c("FP+", "FP-")

tab2_Chmy12 <- rbind(tab2_Chmy12[2,], tab2_Chmy12[1,]) # rearrange rows
rownames(tab2_Chmy12) <- c("Chmy12 +", "Chmy12 -")

tab2_Chmy12

barplot(tab2_Chmy12, beside = TRUE, main = "FP in turtles with and without Chmy12 allele", legend.text = c("Chmy12 +", "Chmy12 -"))

Chmy12_risk <- epi.2by2(tab2_Chmy12, method = "cohort.count", conf.level = 0.998)

Chmy12_risk

df_Chmy12 <- as.data.frame(Chmy12_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy12)

# Chmy21
Cm_data[Cm_data$Chmy21 == 1, "Chmy21"] <- "Chmy21 +"
Cm_data[Cm_data$Chmy21 == 0, "Chmy21"] <- "Chmy21 -"

tab_Chmy21 <- table(Chmy21, FP) # independent, dependent
tab_Chmy21

# rearrange so it's in ABCD format
tab2_Chmy21 <- cbind(tab_Chmy21[,2], tab_Chmy21[,1]) # rearrange columns
colnames(tab2_Chmy21) <- c("FP+", "FP-")

tab2_Chmy21 <- rbind(tab2_Chmy21[2,], tab2_Chmy21[1,]) # rearrange rows
rownames(tab2_Chmy21) <- c("Chmy21 +", "Chmy21 -")

tab2_Chmy21

barplot(tab2_Chmy21, beside = TRUE, main = "FP in turtles with and without Chmy21 allele", legend.text = c("Chmy21 +", "Chmy21 -"))

Chmy21_risk <- epi.2by2(tab2_Chmy21, method = "cohort.count", conf.level = 0.998)

Chmy21_risk

df_Chmy21 <- as.data.frame(Chmy21_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy21)

# Chmy17
Cm_data[Cm_data$Chmy17 == 1, "Chmy17"] <- "Chmy17 +"
Cm_data[Cm_data$Chmy17 == 0, "Chmy17"] <- "Chmy17 -"

tab_Chmy17 <- table(Chmy17, FP) # independent, dependent
tab_Chmy17

# rearrange so it's in ABCD format
tab2_Chmy17 <- cbind(tab_Chmy17[,2], tab_Chmy17[,1]) # rearrange columns
colnames(tab2_Chmy17) <- c("FP+", "FP-")

tab2_Chmy17 <- rbind(tab2_Chmy17[2,], tab2_Chmy17[1,]) # rearrange rows
rownames(tab2_Chmy17) <- c("Chmy17 +", "Chmy17 -")

tab2_Chmy17

barplot(tab2_Chmy17, beside = TRUE, main = "FP in turtles with and without Chmy17 allele", legend.text = c("Chmy17 +", "Chmy17 -"))

Chmy17_risk <- epi.2by2(tab2_Chmy17, method = "cohort.count", conf.level = 0.998)

Chmy17_risk

df_Chmy17 <- as.data.frame(Chmy17_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy17)


# Chmy48
Cm_data[Cm_data$Chmy48 == 1, "Chmy48"] <- "Chmy48 +"
Cm_data[Cm_data$Chmy48 == 0, "Chmy48"] <- "Chmy48 -"

tab_Chmy48 <- table(Chmy48, FP) # independent, dependent
tab_Chmy48

# rearrange so it's in ABCD format
tab2_Chmy48 <- cbind(tab_Chmy48[,2], tab_Chmy48[,1]) # rearrange columns
colnames(tab2_Chmy48) <- c("FP+", "FP-")

tab2_Chmy48 <- rbind(tab2_Chmy48[2,], tab2_Chmy48[1,]) # rearrange rows
rownames(tab2_Chmy48) <- c("Chmy48 +", "Chmy48 -")

tab2_Chmy48

barplot(tab2_Chmy48, beside = TRUE, main = "FP in turtles with and without Chmy48 allele", legend.text = c("Chmy48 +", "Chmy48 -"))

Chmy48_risk <- epi.2by2(tab2_Chmy48, method = "cohort.count", conf.level = 0.998)

Chmy48_risk

df_Chmy48 <- as.data.frame(Chmy48_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy48)

# Chmy07
Cm_data[Cm_data$Chmy07 == 1, "Chmy07"] <- "Chmy07 +"
Cm_data[Cm_data$Chmy07 == 0, "Chmy07"] <- "Chmy07 -"

tab_Chmy07 <- table(Chmy07, FP) # independent, dependent
tab_Chmy07

# rearrange so it's in ABCD format
tab2_Chmy07 <- cbind(tab_Chmy07[,2], tab_Chmy07[,1]) # rearrange columns
colnames(tab2_Chmy07) <- c("FP+", "FP-")

tab2_Chmy07 <- rbind(tab2_Chmy07[2,], tab2_Chmy07[1,]) # rearrange rows
rownames(tab2_Chmy07) <- c("Chmy07 +", "Chmy07 -")

tab2_Chmy07

barplot(tab2_Chmy07, beside = TRUE, main = "FP in turtles with and without Chmy07 allele", legend.text = c("Chmy07 +", "Chmy07 -"))

Chmy07_risk <- epi.2by2(tab2_Chmy07, method = "cohort.count", conf.level = 0.998)

Chmy07_risk

df_Chmy07 <- as.data.frame(Chmy07_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy07)

# Chmy18
Cm_data[Cm_data$Chmy18 == 1, "Chmy18"] <- "Chmy18 +"
Cm_data[Cm_data$Chmy18 == 0, "Chmy18"] <- "Chmy18 -"

tab_Chmy18 <- table(Chmy18, FP) # independent, dependent
tab_Chmy18

# rearrange so it's in ABCD format
tab2_Chmy18 <- cbind(tab_Chmy18[,2], tab_Chmy18[,1]) # rearrange columns
colnames(tab2_Chmy18) <- c("FP+", "FP-")

tab2_Chmy18 <- rbind(tab2_Chmy18[2,], tab2_Chmy18[1,]) # rearrange rows
rownames(tab2_Chmy18) <- c("Chmy18 +", "Chmy18 -")

tab2_Chmy18

barplot(tab2_Chmy18, beside = TRUE, main = "FP in turtles with and without Chmy18 allele", legend.text = c("Chmy18 +", "Chmy18 -"))

Chmy18_risk <- epi.2by2(tab2_Chmy18, method = "cohort.count", conf.level = 0.998)

Chmy18_risk

df_Chmy18 <- as.data.frame(Chmy18_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy18)

# Chmy33
Cm_data[Cm_data$Chmy33 == 1, "Chmy33"] <- "Chmy33 +"
Cm_data[Cm_data$Chmy33 == 0, "Chmy33"] <- "Chmy33 -"

tab_Chmy33 <- table(Chmy33, FP) # independent, dependent
tab_Chmy33

# rearrange so it's in ABCD format
tab2_Chmy33 <- cbind(tab_Chmy33[,2], tab_Chmy33[,1]) # rearrange columns
colnames(tab2_Chmy33) <- c("FP+", "FP-")

tab2_Chmy33 <- rbind(tab2_Chmy33[2,], tab2_Chmy33[1,]) # rearrange rows
rownames(tab2_Chmy33) <- c("Chmy33 +", "Chmy33 -")

tab2_Chmy33

barplot(tab2_Chmy33, beside = TRUE, main = "FP in turtles with and without Chmy33 allele", legend.text = c("Chmy33 +", "Chmy33 -"))

Chmy33_risk <- epi.2by2(tab2_Chmy33, method = "cohort.count", conf.level = 0.998)

Chmy33_risk

df_Chmy33 <- as.data.frame(Chmy33_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy33)

# Chmy49
Cm_data[Cm_data$Chmy49 == 1, "Chmy49"] <- "Chmy49 +"
Cm_data[Cm_data$Chmy49 == 0, "Chmy49"] <- "Chmy49 -"

tab_Chmy49 <- table(Chmy49, FP) # independent, dependent
tab_Chmy49

# rearrange so it's in ABCD format
tab2_Chmy49 <- cbind(tab_Chmy49[,2], tab_Chmy49[,1]) # rearrange columns
colnames(tab2_Chmy49) <- c("FP+", "FP-")

tab2_Chmy49 <- rbind(tab2_Chmy49[2,], tab2_Chmy49[1,]) # rearrange rows
rownames(tab2_Chmy49) <- c("Chmy49 +", "Chmy49 -")

tab2_Chmy49

barplot(tab2_Chmy49, beside = TRUE, main = "FP in turtles with and without Chmy49 allele", legend.text = c("Chmy49 +", "Chmy49 -"))

Chmy49_risk <- epi.2by2(tab2_Chmy49, method = "cohort.count", conf.level = 0.998)

Chmy49_risk

df_Chmy49 <- as.data.frame(Chmy49_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy49)

# Chmy08
Cm_data[Cm_data$Chmy08 == 1, "Chmy08"] <- "Chmy08 +"
Cm_data[Cm_data$Chmy08 == 0, "Chmy08"] <- "Chmy08 -"

tab_Chmy08 <- table(Chmy08, FP) # independent, dependent
tab_Chmy08

# rearrange so it's in ABCD format
tab2_Chmy08 <- cbind(tab_Chmy08[,2], tab_Chmy08[,1]) # rearrange columns
colnames(tab2_Chmy08) <- c("FP+", "FP-")

tab2_Chmy08 <- rbind(tab2_Chmy08[2,], tab2_Chmy08[1,]) # rearrange rows
rownames(tab2_Chmy08) <- c("Chmy08 +", "Chmy08 -")

tab2_Chmy08

barplot(tab2_Chmy08, beside = TRUE, main = "FP in turtles with and without Chmy08 allele", legend.text = c("Chmy08 +", "Chmy08 -"))

Chmy08_risk <- epi.2by2(tab2_Chmy08, method = "cohort.count", conf.level = 0.998)

Chmy08_risk

df_Chmy08 <- as.data.frame(Chmy08_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy08)

# Chmy13
Cm_data[Cm_data$Chmy13 == 1, "Chmy13"] <- "Chmy13 +"
Cm_data[Cm_data$Chmy13 == 0, "Chmy13"] <- "Chmy13 -"

tab_Chmy13 <- table(Chmy13, FP) # independent, dependent
tab_Chmy13

# rearrange so it's in ABCD format
tab2_Chmy13 <- cbind(tab_Chmy13[,2], tab_Chmy13[,1]) # rearrange columns
colnames(tab2_Chmy13) <- c("FP+", "FP-")

tab2_Chmy13 <- rbind(tab2_Chmy13[2,], tab2_Chmy13[1,]) # rearrange rows
rownames(tab2_Chmy13) <- c("Chmy13 +", "Chmy13 -")

tab2_Chmy13

barplot(tab2_Chmy13, beside = TRUE, main = "FP in turtles with and without Chmy13 allele", legend.text = c("Chmy13 +", "Chmy13 -"))

Chmy13_risk <- epi.2by2(tab2_Chmy13, method = "cohort.count", conf.level = 0.998)

Chmy13_risk

df_Chmy13 <- as.data.frame(Chmy13_risk$massoc$RR.strata.wald)


Cm_risks <- rbind(Cm_risks, df_Chmy13)

# Chmy14
Cm_data[Cm_data$Chmy14 == 1, "Chmy14"] <- "Chmy14 +"
Cm_data[Cm_data$Chmy14 == 0, "Chmy14"] <- "Chmy14 -"

tab_Chmy14 <- table(Chmy14, FP) # independent, dependent
tab_Chmy14

# rearrange so it's in ABCD format
tab2_Chmy14 <- cbind(tab_Chmy14[,2], tab_Chmy14[,1]) # rearrange columns
colnames(tab2_Chmy14) <- c("FP+", "FP-")

tab2_Chmy14 <- rbind(tab2_Chmy14[2,], tab2_Chmy14[1,]) # rearrange rows
rownames(tab2_Chmy14) <- c("Chmy14 +", "Chmy14 -")

tab2_Chmy14

barplot(tab2_Chmy14, beside = TRUE, main = "FP in turtles with and without Chmy14 allele", legend.text = c("Chmy14 +", "Chmy14 -"))

Chmy14_risk <- epi.2by2(tab2_Chmy14, method = "cohort.count", conf.level = 0.998)

Chmy14_risk

df_Chmy14 <- as.data.frame(Chmy14_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy14)

# Chmy20
Cm_data[Cm_data$Chmy20 == 1, "Chmy20"] <- "Chmy20 +"
Cm_data[Cm_data$Chmy20 == 0, "Chmy20"] <- "Chmy20 -"

tab_Chmy20 <- table(Chmy20, FP) # independent, dependent
tab_Chmy20

# rearrange so it's in ABCD format
tab2_Chmy20 <- cbind(tab_Chmy20[,2], tab_Chmy20[,1]) # rearrange columns
colnames(tab2_Chmy20) <- c("FP+", "FP-")

tab2_Chmy20 <- rbind(tab2_Chmy20[2,], tab2_Chmy20[1,]) # rearrange rows
rownames(tab2_Chmy20) <- c("Chmy20 +", "Chmy20 -")

tab2_Chmy20

barplot(tab2_Chmy20, beside = TRUE, main = "FP in turtles with and without Chmy20 allele", legend.text = c("Chmy20 +", "Chmy20 -"))

Chmy20_risk <- epi.2by2(tab2_Chmy20, method = "cohort.count", conf.level = 0.998)

Chmy20_risk

df_Chmy20 <- as.data.frame(Chmy20_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy20)

# Chmy50
Cm_data[Cm_data$Chmy50 == 1, "Chmy50"] <- "Chmy50 +"
Cm_data[Cm_data$Chmy50 == 0, "Chmy50"] <- "Chmy50 -"

tab_Chmy50 <- table(Chmy50, FP) # independent, dependent
tab_Chmy50

# rearrange so it's in ABCD format
tab2_Chmy50 <- cbind(tab_Chmy50[,2], tab_Chmy50[,1]) # rearrange columns
colnames(tab2_Chmy50) <- c("FP+", "FP-")

tab2_Chmy50 <- rbind(tab2_Chmy50[2,], tab2_Chmy50[1,]) # rearrange rows
rownames(tab2_Chmy50) <- c("Chmy50 +", "Chmy50 -")

tab2_Chmy50

barplot(tab2_Chmy50, beside = TRUE, main = "FP in turtles with and without Chmy50 allele", legend.text = c("Chmy50 +", "Chmy50 -"))

Chmy50_risk <- epi.2by2(tab2_Chmy50, method = "cohort.count", conf.level = 0.998)

Chmy50_risk

df_Chmy50 <- as.data.frame(Chmy50_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy50)

# Chmy35
Cm_data[Cm_data$Chmy35 == 1, "Chmy35"] <- "Chmy35 +"
Cm_data[Cm_data$Chmy35 == 0, "Chmy35"] <- "Chmy35 -"

tab_Chmy35 <- table(Chmy35, FP) # independent, dependent
tab_Chmy35

# rearrange so it's in ABCD format
tab2_Chmy35 <- cbind(tab_Chmy35[,2], tab_Chmy35[,1]) # rearrange columns
colnames(tab2_Chmy35) <- c("FP+", "FP-")

tab2_Chmy35 <- rbind(tab2_Chmy35[2,], tab2_Chmy35[1,]) # rearrange rows
rownames(tab2_Chmy35) <- c("Chmy35 +", "Chmy35 -")

tab2_Chmy35

barplot(tab2_Chmy35, beside = TRUE, main = "FP in turtles with and without Chmy35 allele", legend.text = c("Chmy35 +", "Chmy35 -"))

Chmy35_risk <- epi.2by2(tab2_Chmy35, method = "cohort.count", conf.level = 0.998)

Chmy35_risk

df_Chmy35 <- as.data.frame(Chmy35_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy35)

# Chmy62
Cm_data[Cm_data$Chmy62 == 1, "Chmy62"] <- "Chmy62 +"
Cm_data[Cm_data$Chmy62 == 0, "Chmy62"] <- "Chmy62 -"

tab_Chmy62 <- table(Chmy62, FP) # independent, dependent
tab_Chmy62

# rearrange so it's in ABCD format
tab2_Chmy62 <- cbind(tab_Chmy62[,2], tab_Chmy62[,1]) # rearrange columns
colnames(tab2_Chmy62) <- c("FP+", "FP-")

tab2_Chmy62 <- rbind(tab2_Chmy62[2,], tab2_Chmy62[1,]) # rearrange rows
rownames(tab2_Chmy62) <- c("Chmy62 +", "Chmy62 -")

tab2_Chmy62

barplot(tab2_Chmy62, beside = TRUE, main = "FP in turtles with and without Chmy62 allele", legend.text = c("Chmy62 +", "Chmy62 -"))

Chmy62_risk <- epi.2by2(tab2_Chmy62, method = "cohort.count", conf.level = 0.998)

Chmy62_risk

df_Chmy62 <- as.data.frame(Chmy62_risk$massoc$RR.strata.wald)

Cm_risks <- rbind(Cm_risks, df_Chmy62)

# Read in row names for graph
Cm_allele <- read.csv("Cm_topalleles_rownames.csv")

risk_ratio_top_Cm_alleles <- cbind(Cm_allele, Cm_risks)

# graph a forest plot

ggplot(data=risk_ratio_top_Cm_alleles, aes(x=allele, y=est, ymin=lower, ymax=upper)) +
  geom_pointrange(shape = 18) +
  geom_hline(yintercept=1, lty=2, color = "red") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("allele") + ylab("relative risk estimate of FP (99.8% CI)") +
  theme_bw() # use a white background

