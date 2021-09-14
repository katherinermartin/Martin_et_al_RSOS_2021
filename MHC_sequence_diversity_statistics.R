library(tidyverse)
library(ape)
library(pegas)
library(ape)

# Population genetics stats on any MHC alleles found in Florida C. caretta (i.e., not necessarily exclusive to C. caretta, including alleles found also in C. mydas and Cape Verde turtles)

Cc_mydna <- read.dna("classI_Florida_Ccaretta.fasta", format = "fasta")

Cc_myalign <- as.alignment(Cc_mydna) # convert dna object into alignment object

# Number of segregating sites:
Cc_segregating <- length(seg.sites(Cc_mydna))
Cc_segregating # 95 segregating sites, across 32 alleles for C. caretta

# set length of nucleotide sequence
Cc_length <- 162

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cc_percent_segregating <- Cc_segregating/Cc_length
Cc_percent_segregating # 0.5864198 or 58.6% segregating sites.

# Calculate nucleotide diversity or pi

nuc.div(Cc_mydna, TRUE) # pi = 0.23477947; variance = 0.01359077

# Population genetics stats on MHC alleles found in Florida C. mydas (i.e., not necessarily exclusive to C. mydas, including alleles found also in Florida and Cape Verde C. caretta)

Cm_mydna <- read.dna("classI_Florida_Cmydas.fasta", format = "fasta")

Cm_myalign <- as.alignment(Cm_mydna) # convert dna object into alignment object

# Number of segregating sites:
Cm_segregating <- length(seg.sites(Cm_mydna))
Cm_segregating # 106 segregating sites, across 98 alleles.

# set nucleotide sequence length
Cm_length <- 162

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cm_percent_segregating <- Cm_segregating/Cm_length
Cm_percent_segregating # 0.654321

# Calculate nucleotide diversity or pi
nuc.div(Cm_mydna, TRUE) # pi = 0.201021837; variance = 0.009590092


