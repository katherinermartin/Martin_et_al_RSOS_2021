# Risk ratio of top alleles and tumor texture in C. mydas with FP, with Bonferroni correction

# CI will be 99.2%; familywise error rate is 95% so 1-(alpha/m) where m = 6 comparisons so 1-(0.05/6) = 0.992

library(epiR)
library(tidyverse)
library(ggplot2)
library(svglite)

# top alleles (those alleles that occur in 10 or more individuals of that species) in Cm with FP so that I can assess relative risk of regression.

# Read in data
data <- read.csv("classI_juveniles_morpho_FP_v3.csv")

# filter to C. mydas records with FP
Cm_data <- data %>% filter(species == "Chelonia mydas") %>% # filter to just Cm
  filter(FP == 1) # filter to just FP + Cms, 105 total


# Subset to the alleles that occur in 10 or more greens with FP and the columns of sample_ID and paps_texture (0s and 1s) but include Chmy33 (only occurs in 4 greens with FP), since it had very high relative importance in the random forest model.
Cm_data <- subset(Cm_data, select = c(sample_ID, paps_texture, Chmy02, Chmy01, Chmy10, Chmy04, Chmy05, Chmy33))


Cm_data[Cm_data$paps_texture == 1, "paps_texture"] <- "smooth" # replace 1s with smooth and 0s with not smooth for pap regression state
Cm_data[Cm_data$paps_texture == 0, "paps_texture"] <- "rough"

attach(Cm_data)

# Chmy01
Cm_data[Cm_data$Chmy01 == 1, "Chmy01"] <- "Chmy01 +"
Cm_data[Cm_data$Chmy01 == 0, "Chmy01"] <- "Chmy01 -"

tab_Chmy01 <- table(Chmy01, paps_texture) # independent, dependent
tab_Chmy01

# rearrange so it's in ABCD format
tab2_Chmy01 <- cbind(tab_Chmy01[,2], tab_Chmy01[,1]) # rearrange columns
colnames(tab2_Chmy01) <- c("smooth", "rough")

tab2_Chmy01 <- rbind(tab2_Chmy01[2,], tab2_Chmy01[1,]) # rearrange rows
rownames(tab2_Chmy01) <- c("Chmy01 +", "Chmy01 -")

tab2_Chmy01

barplot(tab2_Chmy01, beside = TRUE, main = "FP texture in turtles with and without Chmy01 allele", legend.text = c("Chmy01 +", "Chmy01 -"))

Chmy01_risk <- epi.2by2(tab2_Chmy01, method = "cohort.count", conf.level = 0.992)

Chmy01_risk

df_Chmy01 <- as.data.frame(Chmy01_risk$massoc.summary)

# Chmy02
Cm_data[Cm_data$Chmy02 == 1, "Chmy02"] <- "Chmy02 +"
Cm_data[Cm_data$Chmy02 == 0, "Chmy02"] <- "Chmy02 -"

tab_Chmy02 <- table(Chmy02, paps_texture) # independent, dependent
tab_Chmy02

# rearrange so it's in ABCD format
tab2_Chmy02 <- cbind(tab_Chmy02[,2], tab_Chmy02[,1]) # rearrange columns
colnames(tab2_Chmy02) <- c("smooth", "rough")

tab2_Chmy02 <- rbind(tab2_Chmy02[2,], tab2_Chmy02[1,]) # rearrange rows
rownames(tab2_Chmy02) <- c("Chmy02 +", "Chmy02 -")

tab2_Chmy02

barplot(tab2_Chmy02, beside = TRUE, main = "FP texture in turtles with and without Chmy02 allele", legend.text = c("Chmy02 +", "Chmy02 -"))

Chmy02_risk <- epi.2by2(tab2_Chmy02, method = "cohort.count", conf.level = 0.992)

Chmy02_risk

df_Chmy02 <- as.data.frame(Chmy02_risk$massoc.summary)

Cm_risks <- rbind(df_Chmy01, df_Chmy02)

# Chmy04
Cm_data[Cm_data$Chmy04 == 1, "Chmy04"] <- "Chmy04 +"
Cm_data[Cm_data$Chmy04 == 0, "Chmy04"] <- "Chmy04 -"

tab_Chmy04 <- table(Chmy04, paps_texture) # independent, dependent
tab_Chmy04

# rearrange so it's in ABCD format
tab2_Chmy04 <- cbind(tab_Chmy04[,2], tab_Chmy04[,1]) # rearrange columns
colnames(tab2_Chmy04) <- c("smooth", "rough")

tab2_Chmy04 <- rbind(tab2_Chmy04[2,], tab2_Chmy04[1,]) # rearrange rows
rownames(tab2_Chmy04) <- c("Chmy04 +", "Chmy04 -")

tab2_Chmy04

barplot(tab2_Chmy04, beside = TRUE, main = "FP texture in turtles with and without Chmy04 allele", legend.text = c("Chmy04 +", "Chmy04 -"))

Chmy04_risk <- epi.2by2(tab2_Chmy04, method = "cohort.count", conf.level = 0.992)

Chmy04_risk

df_Chmy04 <- as.data.frame(Chmy04_risk$massoc.summary)


Cm_risks <- rbind(Cm_risks, df_Chmy04)

# Chmy05
Cm_data[Cm_data$Chmy05 == 1, "Chmy05"] <- "Chmy05 +"
Cm_data[Cm_data$Chmy05 == 0, "Chmy05"] <- "Chmy05 -"

tab_Chmy05 <- table(Chmy05, paps_texture) # independent, dependent
tab_Chmy05

# rearrange so it's in ABCD format
tab2_Chmy05 <- cbind(tab_Chmy05[,2], tab_Chmy05[,1]) # rearrange columns
colnames(tab2_Chmy05) <- c("smooth", "rough")

tab2_Chmy05 <- rbind(tab2_Chmy05[2,], tab2_Chmy05[1,]) # rearrange rows
rownames(tab2_Chmy05) <- c("Chmy05 +", "Chmy05 -")

tab2_Chmy05

barplot(tab2_Chmy05, beside = TRUE, main = "FP texture in turtles with and without Chmy05 allele", legend.text = c("Chmy05 +", "Chmy05 -"))

Chmy05_risk <- epi.2by2(tab2_Chmy05, method = "cohort.count", conf.level = 0.992)

Chmy05_risk

df_Chmy05 <- as.data.frame(Chmy05_risk$massoc.summary)

Cm_risks <- rbind(Cm_risks, df_Chmy05)

# Chmy10
Cm_data[Cm_data$Chmy10 == 1, "Chmy10"] <- "Chmy10 +"
Cm_data[Cm_data$Chmy10 == 0, "Chmy10"] <- "Chmy10 -"

tab_Chmy10 <- table(Chmy10, paps_texture) # independent, dependent
tab_Chmy10

# rearrange so it's in ABCD format
tab2_Chmy10 <- cbind(tab_Chmy10[,2], tab_Chmy10[,1]) # rearrange columns
colnames(tab2_Chmy10) <- c("smooth", "rough")

tab2_Chmy10 <- rbind(tab2_Chmy10[2,], tab2_Chmy10[1,]) # rearrange rows
rownames(tab2_Chmy10) <- c("Chmy10 +", "Chmy10 -")

tab2_Chmy10

barplot(tab2_Chmy10, beside = TRUE, main = "FP texture in turtles with and without Chmy10 allele", legend.text = c("Chmy10 +", "Chmy10 -"))

Chmy10_risk <- epi.2by2(tab2_Chmy10, method = "cohort.count", conf.level = 0.992)

Chmy10_risk

df_Chmy10 <- as.data.frame(Chmy10_risk$massoc.summary)

Cm_risks <- rbind(Cm_risks, df_Chmy10)

# Chmy33
Cm_data[Cm_data$Chmy33 == 1, "Chmy33"] <- "Chmy33 +"
Cm_data[Cm_data$Chmy33 == 0, "Chmy33"] <- "Chmy33 -"

tab_Chmy33 <- table(Chmy33, paps_texture) # independent, dependent
tab_Chmy33

# rearrange so it's in ABCD format
tab2_Chmy33 <- cbind(tab_Chmy33[,2], tab_Chmy33[,1]) # rearrange columns
colnames(tab2_Chmy33) <- c("smooth", "rough")

tab2_Chmy33 <- rbind(tab2_Chmy33[2,], tab2_Chmy33[1,]) # rearrange rows
rownames(tab2_Chmy33) <- c("Chmy33 +", "Chmy33 -")

tab2_Chmy33

barplot(tab2_Chmy33, beside = TRUE, main = "FP texture in turtles with and without Chmy33 allele", legend.text = c("Chmy33 +", "Chmy33 -"))


# Since one of the categories is zero (Chmy33 - and FP + ), do a Haldane-Anscombe correction: add 1 to all categories

tab2_HA_Chmy33 <- matrix(c(5, 1, 25, 78), ncol =2, byrow = TRUE)
tab2_HA_Chmy33
colnames(tab2_HA_Chmy33) <- c("smooth", "rough")
rownames(tab2_HA_Chmy33) <- c("Chmy33+", "Chmy33-")
tab2_HA_Chmy33

barplot(tab2_HA_Chmy33, beside = TRUE, main = "Regression in turtles with and without Chmy33 allele, HA correction", legend.text = c("Chmy33 +", "Chmy33 -"))

HA_Chmy33_risk <- epi.2by2(tab2_HA_Chmy33, method = "cohort.count", conf.level = 0.992)

HA_Chmy33_risk

df_HA_Chmy33 <- as.data.frame(HA_Chmy33_risk$massoc.summary)


Cm_risks <- rbind(Cm_risks, df_HA_Chmy33)

Cm_risks <- Cm_risks %>% filter(var == "Inc risk ratio") # filter  to just "inc risk ratio"

Cm_risks$allele <- c("Chmy01", "Chmy02", "Chmy04", "Chmy05", "Chmy10", "Chmy33")

# graph a forest plot

ggplot(data=Cm_risks, aes(x=allele, y=est, ymin=lower, ymax=upper)) +
  geom_pointrange(shape = 18) +
  geom_hline(yintercept=1, lty=2, color = "red") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("allele") + ylab("relative risk estimate of FP texture (99.2% CI)") +
  theme_bw() # use a white background
