# distribution allele count for C. mydas and C. caretta

library(ggplot2)
library(dplyr)

# read in data:

df <- read.csv("classI_juveniles_morpho_FP.csv")
mean <- mean(df$allele_count)
mean # avg allele count is 3.57 across all animals
sd(df$allele_count) #1.405

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
