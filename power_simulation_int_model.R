library(MASS)
library(dplyr)
library(qvalue)

#' Generate n genotypes with a minor allele frequency
#'
#' @param maf Minor allele frequency. Must be in [0, 1]
#' @param n Total number of genotypes to generate
#'
#' @return Vector of genotypes, where 0 is homozygous major allele,
#' .5 is heterozygous, and 1 is homozygous minor allele
#' @export
#'
#' @examples
generate_genotypes <- function(maf, n) {

  sample(
    x=c(0, 0.5, 1),
    size=n,
    prob=c((1 - maf) ^ 2, 2 * maf * (1 - maf), maf ^ 2),
    replace = TRUE
  )

}

#' Generate Allele Frequency Changes for Simulated Dataset
#'
#' The function assumes that the true allele frequency changes for each site follows one of a few cases:
#' (1) The minor allele frequency change is 0 in both cases.
#' (2) The minor allele frequency changes by \code{cx} percent in the HS case only.
#' (3) The minor allele frequency changes by \code{cx} percent in the C case only.
#' (4) The minor allele frequency changes by \code{cx} percent in both cases.
#' (5) The minor allele frequency changes by \code{cx} percent in C and \code{cx} * \code{hs_amp_coef} in HS.
#' The user is able to input the magnitude of changes as well as the probability of landing in each case.
#'
#' @param sites Number of SNPs for which to generate data
#' @param prob_both_null Probability of being in case 1 in the description
#' @param prob_hs_only Probability of being in case 2 in the description
#' @param prob_c_only Probability of being in case 3 in the description
#' @param prob_same Probability of being in case 4 in the description
#' @param prob_amp Probability of being in case 5 in the description
#' @param cx Percentage change of minor allele frequency at time TN
#' @param hs_amp_coef Amplification constant, see description for more details
#'
#' @return Dataframe with minor allele frequency changes for each site in the HS and C cases.
#' @export
#'
#' @examples
generate_tn_allele_cxs <- function(
  sites,
  prob_both_null,
  prob_hs_only,
  prob_c_only,
  prob_same,
  prob_amp,
  cx,
  hs_amp_coef = 1.0
) {

  # Create matrices to select from so that only one call to sample is made
  # c(both_null, hs_only, c_only, same, amp)
  c_matrix_data <- c(rep(0, 2 * sites), rep(cx, 3 * sites))
  c_matrix <- matrix(data = c_matrix_data, ncol = sites, byrow = TRUE)
  hs_matrix_data <- c(rep(0, sites), rep(cx, sites), rep(0, sites), rep(cx, sites), rep(cx * hs_amp_coef, sites))
  hs_matrix <- matrix(data = hs_matrix_data, ncol = sites, byrow = TRUE)

  # Select from 5 cases with specified probability
  case_gen <- sample(
    x = c(1:5),
    size = sites,
    replace = TRUE,
    prob = c(prob_both_null, prob_hs_only, prob_c_only, prob_same, prob_amp)
  )
  str_possibilities <- c("both_null", "hs_only", "c_only", "same", "hs_amp")
  str_selections <- str_possibilities[case_gen]

  # Find proper indices in matrix
  selected_indices <- c(0:(sites - 1)) * 5 + case_gen
  selected_c_cxs <- c_matrix[selected_indices]
  selected_hs_cxs <- hs_matrix[selected_indices]

  maf_cx_tn_df <- data.frame(c_cx = selected_c_cxs, hs_cx = selected_hs_cxs)
  return(list(cx_df = maf_cx_tn_df, true_signal = str_selections))

}

#' Generate Simulated Drosophila Dataset With Specified Signal Setup
#'
#' The function assumes that the true allele frequency changes for each site follows one of a few cases:
#' (1) The minor allele frequency change is 0 in both cases.
#' (2) The minor allele frequency changes by \code{cx} percent in the HS case only.
#' (3) The minor allele frequency changes by \code{cx} percent in the C case only.
#' (4) The minor allele frequency changes by \code{cx} percent in both cases.
#' (5) The minor allele frequency changes by \code{cx} percent in C and \code{cx} * \code{hs_amp_coef} in HS.
#' The user is able to input the magnitude of changes as well as the probability of landing in each case.
#'
#' @param n Total number of flies. 1/2 are at time 0, 1/4 are high sugar at time N, 1/4 are control at time N
#' @param sites Total number of SNP sites for each fly
#' @param prob_both_null Probability of being in case 1 in the description
#' @param prob_hs_only Probability of being in case 2 in the description
#' @param prob_c_only Probability of being in case 3 in the description
#' @param prob_same Probability of being in case 4 in the description
#' @param prob_amp Probability of being in case 5 in the description
#' @param cx Percentage change of minor allele frequency at time TN
#' @param hs_amp_coef Amplification constant, see description for more details
#' @param nb_params_dist Vector of lists of negative binomial parameters to be sampled from to determine reads for each SNP.
#' For each SNP, one set of these parameters are sampled, and then a negative binomial distribution with these parameters
#' is used to generate \code{n} read counts for each SNP. Each list should contain an element t corresponding to
#' the \code{n} argument in \code{rnbinom} and \code{p} corresponding to the \code{prob} argument in \code{rnbinom}.
#'
#' @return Dataframe containing information from simulations
#' @export
#'
#' @examples
generate_simulated_dataset <- function(
  n,
  sites,
  maf_dist,
  prob_both_null,
  prob_hs_only,
  prob_c_only,
  prob_same,
  prob_amp,
  cx,
  nb_params_dist,
  hs_amp_coef = 1.0
) {

  # Divide fly counts into groups
  n_t0 <- floor(n / 2)
  n_hs_tn <- floor((n - n_t0) / 2)
  n_c_tn <- n_hs_tn

  # Sample parameters before loop for speed
  maf_t0_freqs <- sample(x = maf_dist, size = sites, replace = TRUE)
  freq_change_list <- generate_tn_allele_cxs(
    sites, prob_both_null, prob_hs_only, prob_c_only, prob_same, prob_amp, cx, hs_amp_coef
  )
  freq_change_df <- freq_change_list$cx_df
  true_signal_df <- data.frame(site = c(1:sites), true_signal = freq_change_list$true_signal)

  maf_hs_tn_freqs <- maf_t0_freqs * (
    1 + freq_change_df$hs_cx
  )
  maf_c_tn_freqs <- maf_t0_freqs * (
    1 + freq_change_df$c_cx
  )
  nb_site_params <- sample(x = nb_params_dist, size = sites, replace = TRUE)

  site_vec <- c()
  total_reads_vec <- c()
  alt_reads_vec <- c()
  condition_vec <- c()
  time_vec <- c()

  for(i in c(1:sites)) {

    # Generate genotypes for HS and Control at T0 and TN
    site_t0_maf <- maf_t0_freqs[i]
    genotypes_t0 <- generate_genotypes(site_t0_maf, n_t0)

    site_tn_hs_maf <- maf_hs_tn_freqs[i]
    genotypes_hs_tn <- generate_genotypes(site_tn_hs_maf, n_hs_tn)

    site_tn_c_maf <- maf_c_tn_freqs[i]
    genotypes_c_tn <- generate_genotypes(site_tn_c_maf, n_c_tn)

    # Sample reads based on negative binomial distribution
    # Make sure there are never 0 reads, this should happen very rarely
    total_reads <- pmax(1, rnbinom(n = n, size = nb_site_params[i][[1]]$t, prob = nb_site_params[i][[1]]$p))
    alt_reads <- rbinom(n = n, size = total_reads, prob = c(genotypes_t0, genotypes_hs_tn, genotypes_c_tn))

    site_vec <- c(site_vec, rep(i, n))
    total_reads_vec <- c(total_reads_vec, total_reads)
    alt_reads_vec <- c(alt_reads_vec, alt_reads)
    condition_vec <- c(condition_vec, rep("both", n_t0), rep("HS", n_hs_tn), rep("C", n_c_tn))
    time_vec <- c(time_vec, rep("T0", n_t0), rep("TN", n_hs_tn + n_c_tn))

  }

  sim_df <- data.frame(
    site = site_vec,
    total_reads = total_reads_vec,
    alt_reads = alt_reads_vec,
    condition = condition_vec,
    time = time_vec
  )

  return(list(sim_df = sim_df, true_signal_df = true_signal_df))

}

test_generate_simulated_dataset <- function() {

  # The point of this function is to build some tests that will allow me to verify the correctness of the function
  # At some point it would probably be better to build a testing suite if this is moved to an R package
  sim_df <- generate_simulated_dataset(
    n = 20,
    sites = 10,
    maf_dist = c(.25, .4),
    prob_both_null = .5,
    prob_hs_only = 0.1,
    prob_c_only = 0.1,
    prob_same = 0.1,
    prob_amp = 0.2,
    cx = -0.15,
    hs_amp_coef = 1.5,
    nb_params_dist = list(list(t = 25, p = .4), list(t = 10, p = .3))
  )$sim_df

  # Test basic expected features of simulated dataset
  num_rows_t0 <- sim_df %>% filter(time == 'T0') %>% nrow()
  stopifnot(num_rows_t0 == 100)

  num_rows_tn <- sim_df %>% filter(time == 'TN') %>% nrow()
  stopifnot(num_rows_tn == 100)

  num_rows_both <- sim_df %>% filter(condition == 'both') %>% nrow()
  stopifnot(num_rows_both == 100)

  num_rows_c <- sim_df %>% filter(condition == 'C') %>% nrow()
  stopifnot(num_rows_c == 50)

  num_rows_hs <- sim_df %>% filter(condition == 'HS') %>% nrow()
  stopifnot(num_rows_hs == 50)

  # Some of the tests below might have some degree of stochasticity
  # But I think they'll be necessary to at least document in order to show the things that I'm calculating
  # In order to make sure that the power simulation is correct

  sim_df2 <- generate_simulated_dataset(
    n = 500,
    sites = 250,
    maf_dist = c(.25, .4),
    prob_both_null = 1.0,
    prob_hs_only = 0.0,
    prob_c_only = 0.0,
    prob_same = 0.0,
    prob_amp = 0.0,
    cx = -0.15,
    nb_params_dist = list(list(t = 25, p = .4), list(t = 10, p = .3))
  )$sim_df

  # Make sure allele frequencies don't vary accross time and condition in the all null case
  sim_df2 <- sim_df2 %>% mutate(allele_prop = alt_reads / total_reads)
  sim_df2_summary <- sim_df2 %>% group_by(time, condition) %>% summarize(maf = mean(allele_prop))
  max_diff <- max(sim_df2_summary$maf) - min(sim_df2_summary$maf)
  stopifnot(max_diff < 0.01)

  sim_df3 <- generate_simulated_dataset(
    n = 500,
    sites = 250,
    maf_dist = c(.25, .4),
    prob_both_null = 0.0,
    prob_hs_only = 0.5,
    prob_c_only = 0.5,
    prob_same = 0.0,
    prob_amp = 0.0,
    cx = -0.15,
    nb_params_dist = list(list(t = 25, p = .4), list(t = 10, p = .3))
  )$sim_df

  sim_df3 <- sim_df3 %>% mutate(allele_prop = alt_reads / total_reads)
  sim_df3_summary <- sim_df3 %>%
    group_by(condition) %>%
    summarize(maf = mean(allele_prop)) %>%
    filter(condition != 'both')
  maf_tn_diff <- max(sim_df3_summary$maf) - min(sim_df3_summary$maf)
  stopifnot(maf_tn_diff < 0.01)

  sim_df4 <- generate_simulated_dataset(
    n = 500,
    sites = 250,
    maf_dist = c(.25, .4),
    prob_both_null = 0.0,
    prob_hs_only = 0.0,
    prob_c_only = 0.0,
    prob_same = 1.0,
    prob_amp = 0.0,
    cx = -0.15,
    nb_params_dist = list(list(t = 25, p = .4), list(t = 10, p = .3))
  )$sim_df

  sim_df4 <- sim_df4 %>% mutate(allele_prop = alt_reads / total_reads)
  sim_df4_summary <- sim_df4 %>%
    group_by(time) %>%
    summarize(maf = mean(allele_prop))
  maf_t0 <- (sim_df4_summary %>% filter(time == 'T0'))$maf[1]
  maf_tn <- (sim_df4_summary %>% filter(time == 'TN'))$maf[1]
  maf_cx <- (maf_tn / maf_t0) - 1
  maf_diff_vs_expected <- maf_cx - .15
  stopifnot(maf_diff_vs_expected < 0.01)

  # Finally, I want to do a test of the amplification
  # Again, check for correct MAFs in this case
  sim_df5 <- generate_simulated_dataset(
    n = 2000,
    sites = 1000,
    maf_dist = c(.25, .4),
    prob_both_null = 0.5,
    prob_hs_only = 0.0,
    prob_c_only = 0.0,
    prob_same = 0.0,
    prob_amp = 0.5,
    hs_amp_coef = 2.0,
    cx = -0.1,
    nb_params_dist = list(list(t = 25, p = .4), list(t = 10, p = .3))
  )$sim_df
  sim_df5 <- sim_df5 %>% mutate(allele_prop = alt_reads / total_reads)
  sim_df5_summary <- sim_df5 %>%
    group_by(time, condition) %>%
    summarize(maf = mean(allele_prop))
  t0_maf <- (sim_df5_summary %>% filter(condition == "both"))$maf[1]
  stopifnot(abs(t0_maf - .325) < .01)
  tn_maf_c <- (sim_df5_summary %>% filter(condition == "C"))$maf[1]
  stopifnot(abs(tn_maf_c - .95 * .325) < .01)
  tn_maf_hs <- (sim_df5_summary %>% filter(condition == "HS"))$maf[1]
  stopifnot(abs(tn_maf_hs - .9 * .325) < .01)

}

#test_generate_simulated_dataset()

sim_df3 <- generate_simulated_dataset(
  n = 500,
  sites = 250,
  maf_dist = c(.25, .4),
  prob_both_null = 0.0,
  prob_hs_only = 0.5,
  prob_c_only = 0.5,
  prob_same = 0.0,
  prob_amp = 0.0,
  cx = -0.15,
  nb_params_dist = list(list(t = 25, p = .4), list(t = 10, p = .3))
)

print(0)

calc_pallares_results <- function(sim_df, n, sites, pval_thresh = .05, qval_thresh = .1) {

  site_vec <- c()
  test_result_vec <- c()

  # Make sure that the sim_df is sorted
  sim_df <- sim_df %>% arrange(site)
  cond_c_pval_vec <- numeric(sites)
  cond_hs_pval_vec <- numeric(sites)

  for(i in c(1:sites)) {

    # slice instead of filter for speed
    site_df <- sim_df %>% slice(c((1 + (i - 1) * n):(i * n)))
    bbn_model <- betabin(cbind(alt_reads, total_reads - alt_reads) ~ condition, ~1, data=site_df)
    mod_sum <- summary(bbn_model)
    cond_c_pval_vec[i] <- attributes(mod_sum)$Coef['conditionC', 4]
    cond_hs_pval_vec[i] <- attributes(mod_sum)$Coef['conditionHS', 4]

  }

  hs_pval_test <- cond_hs_pval_vec < pval_thresh
  c_pval_test <- cond_c_pval_vec < pval_thresh

  hs_qval_test <- qvalue(p = cond_hs_pval_vec, fdr.level = qval_thresh)$significant
  c_qval_test <- qvalue(p = cond_c_pval_vec, fdr.level = qval_thresh)$significant

  # Would like to classify signals into 4 categories
  # (1) both null
  # (2) hs only
  # (3) c only
  # (4) same (or shared)

  # shared signal if the p-value < .05 in one environment and < 10% FDR in the other
  shared_test <- (hs_pval_test & c_qval_test) | (c_pval_test & hs_qval_test)

  # hs only if the fdr test passes for hs but the pval test fails for c
  hs_only_test <- !c_pval_test & hs_qval_test

  # c only if the fdr test passes for c but the pval test fails for hs
  c_only_test <- !hs_pval_test & c_qval_test

  results_df <- data.frame(
    site = c(1:sites),
    shared_test = shared_test,
    hs_only_test = hs_only_test,
    c_only_test = c_only_test
  )

  results_df <- results_df %>%
    mutate(
      classification = case_when(
        shared_test ~ "shared",
        hs_only_test ~ "hs_only",
        c_only_test ~ "c_only",
        TRUE ~ "null"
      )
    )

  results_df <- results_df %>% select(c(site, classification))
  return(results_df)

}

# Now, I have a function that will get me the dataframe as well as the true signals
# The next step is implementing our method as well as the method of the authors
# Once we do that, we'll be able to figure out the relative power / error of the two approaches
# I think the first step is to build a method that will do all the regressions based on the site
# and then go from there

# I feel pretty confident that things are working correctly in the generation of simulated datasets
# But I need to do additional testing in order to make sure that the properties I expect the data to have will actually
# hold. I'm a bit worried that there's a bug somewhere and that the power simulations may be incorrect if so


# Now, it would be good to go through and document the code
# Once that's done I should run it in debug to make sure that things look like I expect they should
# I'm also happy to add some default parameters to make things simpler

# The above procedure works very well to simulate a variety of alternative distributions.
# However, it does not give us any ability to specifically calculate a correlation in signals
# It would be particularly useful to simulate an amplification effect, as that is what we're interested in

# Some are null, some are C only, some are HS only, some are equal, some are amplified in the HS environment
# This code generates MAF changes from this distribution


# I'm a bit worried about the computational complexity here just because of the requirement of permutations to be run 10 times


# Now, I think that what I need to do is create functions that will run the relevant regression on each dataset
# This will allow me to run the power simulation

# One difficulty I am imagining is the various p-value thresholds that one might use. The difficult part here is that
# the authors have somewhat of a piecemeal method for determining different effects. I'm not really sure how well
# they control for type I error, or how well we will. It seems strange in both cases, so I want to investigate that
# It seems like simulation should be a very good way to do this.

# There are also clearly other dimensions of the simulations that I haven't been thinking about.
# What percentage of the signals should be null?
# What percentage of the non-null signals for each environment should be shared?
# I think that I could add a few additional parameters here that should make things much easier
# And then as a result I'll have a much more flexible simulation that will allow me to test more parts of the models
# This would allow us to understand if Type I error is calibrated as well as power.

# There are many dimensions of complexity I can add to the simulation function above
# Basically, the distribution of effects can have a complex structure, and I need to figure out how exactly I want to build that
# structure. I want to be able to sample from some distribution of effects. I'll have to think pretty extensively
# Through the range of possible alternative distributions. I think there are many possibilities here. I'm not sure
# exactly what is the most relevant.

# The next steps will be to (1) Think more about the alternative distributions and (2) Write out the specific procedure for
# both the paper (easier) and our method (more thinking involved)


library(data.table)
library(dplyr)
set.seed(37)

# take samples of 10k just to get frequencies
alt <- fread('data/29Jun20_merged_alt_counts_allCHR.txt') %>% distinct(site, .keep_all = TRUE) %>% sample_n(10000)
both <- fread('data/29Jun20_merged_both_counts_allCHR.txt') %>% filter(site %in% alt$site)

alt_counts <- alt %>% select(-site) %>% rowSums()
both_counts <- both %>% select(-site) %>% rowSums()

alt_count_df <- data.frame(site=alt$site, alt_count = alt_counts)
both_count_df <- data.frame(site=both$site, total_count = both_counts)

count_df <- alt_count_df %>% inner_join(both_count_df, by = c("site"))
count_df <- count_df %>% mutate(alt_freq = alt_count / total_count)
