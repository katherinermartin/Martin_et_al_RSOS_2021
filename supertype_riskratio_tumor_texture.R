# Risk ratio analysis: FP texture risk by supertype, in C. mydas

library(epiR)
library(tidyverse)
library(ggplot2)
library(svglite)

# Load data
data <- read.csv("classI_juveniles_morpho_FP_v3.csv")


# filter to C. mydas records only
Cm_data <- data %>% filter(species == "Chelonia mydas")

# filter to just individuals with FP
Cm_data <- Cm_data %>% filter (FP == 1)

Cm_data[Cm_data$paps_smooth_regressed == 1, "paps_smooth_regressed"] <- "smooth" # replace 1s with regression and 0s with rough
Cm_data[Cm_data$paps_smooth_regressed == 0, "paps_smooth_regressed"] <- "rough"

# remove sea turtles that are truly regressed ("regression inferred from recap records" in column "paps visually categorized as regressed at capture"). Removing from analyses regarding tumor texture. N = 7
Cm_data <- Cm_data %>% filter(paps.visually.categorized.as.regressed.at.capture != c("regression inferred from recap records"))


Cm_data$supertype_A
# st_A
Cm_data[Cm_data$supertype_A == 1, "supertype_A"] <- "st_A +"
Cm_data[Cm_data$supertype_A == 0, "supertype_A"] <- "st_A -"

tab_st_A <- table(Cm_data$supertype_A, Cm_data$paps_smooth_regressed) # independent, dependent
tab_st_A

# rearrange so it's in ABCD format
tab2_st_A <- cbind(tab_st_A[,2], tab_st_A[,1]) # rearrange columns
colnames(tab2_st_A) <- c("smooth", "rough")

tab2_st_A <- rbind(tab2_st_A[2,], tab2_st_A[1,]) # rearrange rows
rownames(tab2_st_A) <- c("st_A +", "st_A -")

tab2_st_A

barplot(tab2_st_A, beside = TRUE, main = "tumor texture in turtles with and without supertype A alleles", legend.text = c("st_A +", "st_A -"))

st_A_risk <- epi.2by2(tab2_st_A, method = "cohort.count", conf.level = 0.95)

st_A_risk$massoc.summary

st_A_risk <- as.data.frame(st_A_risk$massoc.summary)

# st_B
Cm_data[Cm_data$supertype_B == 1, "supertype_B"] <- "st_B +"
Cm_data[Cm_data$supertype_B == 0, "supertype_B"] <- "st_B -"

tab_st_B <- table(Cm_data$supertype_B, Cm_data$paps_smooth_regressed) # independent, dependent
tab_st_B

# rearrange so it's in ABCD format
tab2_st_B <- cbind(tab_st_B[,2], tab_st_B[,1]) # rearrange columns
colnames(tab2_st_B) <- c("smooth", "rough")

tab2_st_B <- rbind(tab2_st_B[2,], tab2_st_B[1,]) # rearrange rows
rownames(tab2_st_B) <- c("st_B +", "st_B -")

tab2_st_B

barplot(tab2_st_B, beside = TRUE, main = "tumor texture in turtles with and without supertype B alleles", legend.text = c("st_B +", "st_B -"))

st_B_risk <- epi.2by2(tab2_st_B, method = "cohort.count", conf.level = 0.95)

st_B_risk

st_B_risk <- as.data.frame(st_B_risk$massoc.summary)

supertype_texture_risk <- rbind(st_A_risk, st_B_risk)

# st_C
Cm_data[Cm_data$supertype_C == 1, "supertype_C"] <- "st_C +"
Cm_data[Cm_data$supertype_C == 0, "supertype_C"] <- "st_C -"

tab_st_C <- table(Cm_data$supertype_C, Cm_data$paps_smooth_regressed) # independent, dependent
tab_st_C

# rearrange so it's in ABCD format
tab2_st_C <- cbind(tab_st_C[,2], tab_st_C[,1]) # rearrange columns
colnames(tab2_st_C) <- c("smooth", "rough")

tab2_st_C <- rbind(tab2_st_C[2,], tab2_st_C[1,]) # rearrange rows
rownames(tab2_st_C) <- c("st_C +", "st_C -")

tab2_st_C

# Since one of the categories is zero (st_C - and smooth), do a Haldane-Anscombe correction: add 1 to all categories

tab2_HA_st_C <- matrix(c(23, 77, 1, 2), ncol =2, byrow = TRUE)
tab2_HA_st_C
colnames(tab2_HA_st_C) <- c("smooth", "rough")
rownames(tab2_HA_st_C) <- c("st_C+", "st_C-")
tab2_HA_st_C

barplot(tab2_HA_st_C, beside = TRUE, main = "tumor texture in turtles with and without supertype C alleles", legend.text = c("st_C +", "st_C -"))

st_C_risk <- epi.2by2(tab2_HA_st_C, method = "cohort.count", conf.level = 0.95)

st_C_risk

st_C_risk <- as.data.frame(st_C_risk$massoc.summary)

supertype_texture_risk <- rbind(supertype_texture_risk, st_C_risk)

supertype_texture_risk <- supertype_texture_risk %>% filter(var == "Inc risk ratio") # filter  to just "inc risk ratio"

supertype_texture_risk$supertype <- c("A", "B", "C")


supertype_texture <- ggplot(data=supertype_texture_risk, aes(x=supertype, y=est, ymin=lower, ymax=upper)) +
  geom_pointrange(shape = 18) +
  geom_hline(yintercept=1, lty=2, color = "red") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("supertype") + ylab("relative risk estimate of smooth texture (95% CI)") +
  theme_bw() # use a white background
supertype_texture



