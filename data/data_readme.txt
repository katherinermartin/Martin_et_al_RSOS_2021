Information on raw data in repository supporting Martin et al, Adaptive evolution of MHC class I immune genes and disease associations in coastal juvenile sea turtles, submitted to Royal Society Open Science.

1. allele_counts_by_species.csv: input for phylogeny.R script, contains information regarding allele counts by species (C. mydas or C. caretta) and location (this study: Florida; or Stiebens et al. 2013 (doi:10.1098/rspb.2013.0305) Cape Verde)

2. classI_amino_acid_alm_to_humanHLA.fasta: alignment of 124 sea turtle MHC alleles (those from this study and those unique to Stiebens et al. 2013) aligned to human HLA-2 for identification of putative peptide binding residues.

3. classI_Florida_Ccaretta.fasta: alignment of MHC alleles recovered in Florida C. caretta from this study for input into MHC_sequence_diversity_statistics.R

4. classI_Florida_Cmydas.fasta: alignment of MHC alleles recovered in Florida C. mydas from this study for input into MHC_sequence_diversity_statistics.R

5. classI_haplotype_for_popart.nex: nexus file for input to PopArt program, with trait block (individuals per allele), with thanks to Dr. Christophe Eizaguirre for providing individuals per allele for Cape Verde loggerheads sampled in Stiebens et al. 2013 (doi:10.1098/rspb.2013.0305) and Stiebens et al. 2013 (doi:10.1186/1471-2148-13-95)

6. classI_juveniles_morpho_FP_v3.csv: dataframe containing all data relevant to project (i.e., demographic information on turtles, MHC alleles, MHC allele count); input to random forest analyses, risk ratio analyses, and allele count distribution script to generate histogram

7. classI_nucleotide_alm.nex: nexus file of alignment of 124 sea turtle MHC class I alleles and 2 outgroups with block information for input to MrBayes.

8. classI_nucleotide_alm.nex.con.tre: Bayesian phylogeny output from RAxML analysis; used as input to phylogeny.R script.

9. classI_nucleotide_alm.phylip: phylip file of alignment of 124 sea turtle MHC class I alleles and 2 outgroups, for input to RAxML.

10. classI_partition1_for_hyphy.nex: nucleotide alignment of MHC alleles and tree information from MrBayes for input to HyPhy. This partition is of nucleotide positions 1-96.

11. classI_partition1.nex: nucleotide alignment of MHC alleles at positions 1-96 for input into MrBayes, to generate tree information for input to HyPhy.

12. classI_partition2_for_hyphy.nex: nucleotide alignment of MHC alleles and tree information from MrBayes for input to HyPhy. This partition is of nucleotide positions 97-159.

13. classI_partition2.nex: nucleotide alignment of MHC alleles at positions 97-159 for input into MrBayes, to generate tree information for input to HyPhy.

14. classI_SinglePartition_for_hyphy.nex: nucleotide alignment of MHC alleles and tree information from MrBayes for input to HyPhy for model identification and recombination breakpoint test.

15. classI_SinglePartition.nex: nucleotide alignment of MHC alleles for input into MrBayes, to generate tree information for input to HyPhy. Alignment was truncated at 159 bp to provide a correct open reading frame for HyPhy.

16. Martin2021_MHC_alleles_sequencelist.fasta: file of unaligned unique MHC alleles (101 total) recovered in this study. Alleles will be uploaded to GenBank and accession IDs listed in publication upon acceptance.

17. MHC_amino_acid_matrix_values_for_DAPC.csv: contains a matrix of 5 values summarizing biochemical values for each amino acid of 124 MHC alleles used in this study; used as input to supertyping_with_DAPC.R script

18. RAxML_bootstrap.classI_alm_raxml: maximum likelihood phylogeny output from RAxML analysis; used as input to phylogeny.R script.
