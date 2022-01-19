# drosophila-longevity
Fall 2021 Project with Harpak Lab at UT Austin

This project used the data from the preprint "Diet unmasks genetic variants that regulate lifespan in outbred Drosophila" (Pallares et al. 2021).

## Experimental Setup
Following is a brief summary of the Pallares et al. (2021) experimental setup. Please refer to the original supplemental materials for a more 
detailed description.
- Flies came from pools A and B.
- 1,000 flies sampled from each pool at beginning of experiment and recorded for "time 0" (flies were 1-3 days old).
- Four cages filled with 10,000 flies each from pool A, and two cages filled with 10,000 flies each from pool B.
- Three cages were fed a "high sugar" diet (22.5% sugar) and three were fed a standard lab diet (7.5% sugar).
- When the last ~500 surviving flies remain in a given cage, they were collected and recorded for "time N" (flies were 32-59 days old).

## Goal of Analysis 
The goal of the project was twofold: 
  (1) Replicate the authors' original analysis, which used a beta binomial to classify SNPs as having a statistically significant effect on lifespan 
  in the high sugar environment (HS), in the control environment (CTRL), in both environments, or in neither.
  (2) Create a new model that adds an interaction term between time and environment, and compare to the original model results.

## The Data
I used four files from the authors in my analysis.
- "29Jun20_merged_alt_counts_allCHR.txt" - 271,536 rows corresponding to SNPs. Each row has allele counts for the 10,000 flies that were mapped to the alternate gene.
- "29Jun20_merged_both_counts_allCHR.txt" - same as above, but the total count mapped.
- "29Jun20_merged_sample_info.txt" - info on the 10,000 flies (timepoint, cage number, starting pool, sex, HS or CTRL, etc).
- "SummaryTable_allsites_12Nov20.txt" - summary of authors' original analysis, including classification of each site as HS, CTRL, both, or NS.

## Description of R Scripts in Repository
- "beta_binomial.R" - this script recreates the original beta-binomial model from Pallares et al.
- "int_model.R" - this script fits the new model with an interaction term for time:env; it was run on TACC using the associated slurm file after splitting the count data in 272 files of 1000 rows each.
- "int_model.slurm" - slurm file for submitting int_model.R to run on TACC
- "plot_int_model_results.R" - reads in the output from int_model.R and creates a variety of plots for analysis
