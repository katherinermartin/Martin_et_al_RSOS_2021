# distribution allele count for C. mydas and C. caretta

library(ggplot2)
library(dplyr)

# read in data:

df <- read.csv("/classI_juveniles_morpho_FP_v3.csv")

mean <- mean(df$allele_count)
mean # avg allele count is 3.57 across all animals
sd(df$allele_count) # standard deviaiton is 1.41

# average allele count per species
# Caretta caretta:
Cc <- df %>% filter(species == "Caretta caretta")
Cc_mean <- mean(Cc$allele_count)
Cc_mean # 3.93
sd(Cc$allele_count) # standard deviation is 1.491

# Chelonia mydas:
Cm <- df %>% filter(species == "Chelonia mydas")
Cm_mean <- mean(Cm$allele_count)
Cm_mean # 3.46
sd(Cm$allele_count) # standard deviation is 1.36

# get frequencies for each category of number of alleles
allele_count <- df %>%
  group_by(allele_count) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n))

# all together:
hist_allele_count <- ggplot(df, aes(x=allele_count)) +
  geom_histogram(bins = 7, fill = "grey", color = "black") +
  geom_vline(aes(xintercept = mean), linetype = "dashed", size = 0.6) +
  theme_bw() +
  scale_x_continuous(breaks=0:7, name ="number of alleles per individual")

hist_allele_count

ggsave(plot = hist_allele_count, filename = "../Analysis/Figures/hist_allele_count.pdf", device = "pdf")
