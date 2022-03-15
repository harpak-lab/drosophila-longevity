# Here, I'd like to reproduce MASH estimates on the Drosophila variants
'%>%' <- dplyr::'%>%'
set.seed(2022)

# First, I will try to fit this on a much smaller sample of the data

# read in data and select relevant columns
summary_table <- read.delim('data/SummaryTable_allsites_12Nov20.txt')
summary_table <- summary_table %>%
  dplyr::sample_n(1000) %>%
  dplyr::select(c(site, pval_CTRL, pval_HS, coef_CTRL, coef_HS, sig_cat))

# replace 0 p-values with small numbers
summary_table <- summary_table %>%
  dplyr::mutate(
    pval_CTRL = pmax(.0001, pval_CTRL),
    pval_HS = pmax(.0001, pval_HS)
  )

# construct std error estimates from coefficients and p-values
summary_table <- summary_table %>%
  dplyr::mutate(
    std_error_ctrl = abs(coef_CTRL) / qnorm((2 - pval_CTRL) / 2),
    std_error_hs = abs(coef_HS) / qnorm((2 - pval_HS) / 2)
  )

reg_fx_mat <- t(matrix(
  data = c(summary_table$coef_CTRL, summary_table$coef_HS),
  nrow = 2
))
colnames(reg_fx_mat) <- c("ctrl", "hs")

reg_se_mat <- t(matrix(
  data = c(summary_table$std_error_ctrl, summary_table$std_error_hs),
  nrow = 2
))
colnames(reg_se_mat) <- c("ctrl", "hs")

mash_data <- mashr::mash_set_data(reg_fx_mat, reg_se_mat)

# Now, want to construct covariance matrices to feed into mash
cov_mat_list <- list()

cov_mat_list[['no_effect']] <- matrix(
  data = rep(0, 4), nrow = 2, dimnames = list(
    rows = c("ctrl", "hs"), cols = c("ctrl", "hs")
  )
)

cov_mat_list[['hs_spec']] <- matrix(
  data = c(0, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(
    rows = c("ctrl", "hs"), cols = c("ctrl", "hs")
  )
)

cov_mat_list[['ctrl_spec']] <- matrix(
  data = c(1, 0, 0, 0), nrow = 2, byrow = TRUE, dimnames = list(
    rows = c("ctrl", "hs"), cols = c("ctrl", "hs")
  )
)

make_amp_cov_mat <- function(
  desired_corr, amp_coef = 1, amp = TRUE, amp_hs = TRUE
) {

  if (amp_hs && amp) {

    ctrl_sd <- 1
    hs_sd <- ctrl_sd * amp_coef

  } else if(!amp_hs && amp) {

    hs_sd <- 1
    ctrl_sd <- hs_sd * amp_coef

  } else {

    hs_sd <- 1
    ctrl_sd <- 1

  }

  # now, have the standard deviations
  cov_hs_ctrl <- desired_corr * hs_sd * ctrl_sd

  cov_mat <- matrix(
    data = c(ctrl_sd ^ 2, cov_hs_ctrl, cov_hs_ctrl, hs_sd ^ 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      rows = c("ctrl", "hs"), cols = c("ctrl", "hs")
    )
  )

  return(cov_mat)

}

desired_corrs <- seq(from = -1, to = 1, by = .25)
desired_amp <- c(3, 2, 1.5)

for(corr in desired_corrs) {

  cov_mat_list[[glue::glue('equal_corr_{corr}')]] <- make_amp_cov_mat(
    desired_corr = corr, amp = FALSE
  )

  for(cond in c("hs", "ctrl")) {

    for(amp in desired_amp) {

      cov_mat_list[[glue::glue('{cond}_amp_{amp}_corr_{corr}')]] <- make_amp_cov_mat(
        desired_corr = corr, amp_hs = (cond == "hs"), amp_coef = amp
      )

    }

  }

}

mash_out <- mashr::mash(
  data = mash_data,
  Ulist = cov_mat_list,
  algorithm.version = "Rcpp",
  outputlevel = 1
)

cov_mat_ests <- mashr::get_estimated_pi(mash_out)

