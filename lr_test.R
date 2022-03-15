library(dplyr)
library(ggplot2)
library(aod)
library(data.table)

set.seed(37)

summary_table <- read.delim('data/SummaryTable_allsites_12Nov20.txt')
summary_table <- summary_table %>%
  select(c(site, pval_CTRL, pval_HS, coef_CTRL, coef_HS, q_perm_CTRL, q_perm_HS, sig_cat)) %>%
  filter(sig_cat != 'NS')

shared_signif_df <- summary_table %>%
  filter(sig_cat == "shared")

alt <- readRDS("rds_data/alt_shared_sig.rds")

both <- readRDS("rds_data/both_shared_sig.rds")

info <- read.delim('data/29Jun20_merged_sample_info.txt')

variant_vec <- c()
pval_vec <- c()
hs_coef_vec <- c()
hs_pval_vec <- c()
c_coef_vec <- c()
c_pval_vec <- c()


for(i in c(1:nrow(alt))) {

  print(i)
  info$tmp_x <- as.vector(t(alt[i,-1]))
  info$tmp_y <- as.vector(t(both[i,-1]))

  tmp_info <- subset(info, tmp_y>0 & sex!='unknown')
  int_model <- betabin(cbind(tmp_x, tmp_y - tmp_x) ~ condition + sequencing_batch + meta_cage + sex, ~1, data=tmp_info, control = list(maxit  = 1000000))
  time_model <- betabin(cbind(tmp_x, tmp_y - tmp_x) ~ timepoint + sequencing_batch + meta_cage + sex, ~1, data=tmp_info, control = list(maxit  = 1000000))

  int_model_ll <- as.numeric(logLik(int_model))
  time_model_ll <- as.numeric(logLik(time_model))

  lr_stat <- -2 * (time_model_ll - int_model_ll)
  lr_pval <- pchisq(lr_stat, df = 1, lower.tail = FALSE)

  variant_vec <- c(variant_vec, alt$site[i])
  pval_vec <- c(pval_vec, lr_pval)
  hs_coef_vec <- c(hs_coef_vec, attributes(summary(int_model))$Coef['conditionHS', 'Estimate'])
  hs_pval_vec <- c(hs_pval_vec, attributes(summary(int_model))$Coef['conditionHS', 'Pr(> |z|)'])
  c_coef_vec <- c(c_coef_vec, attributes(summary(int_model))$Coef['conditionC', 'Estimate'])
  c_pval_vec <- c(c_pval_vec, attributes(summary(int_model))$Coef['conditionC', 'Pr(> |z|)'])

}

pval_df <- data.frame(
  variant = variant_vec,
  lr_pval = pval_vec,
  c_coef = c_coef_vec,
  hs_coef = hs_coef_vec,
  c_coef_pval = c_pval_vec,
  hs_pval = hs_pval_vec
)

saveRDS(pval_df, "rds_data/pval_shared.rds")

library(qqconf)

qq_conf_plot(pval_df$pval, distribution = qunif, log10 = TRUE)

summary_table <- summary_table %>%
  inner_join(pval_df, by = c('site' = 'variant'))




