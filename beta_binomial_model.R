# this script recreates the original beta-binomial model from Pallares et al.

library(data.table)
library(aod)
library(ggplot2)

# count data for beta binomial model
alt <- fread('29Jun20_merged_alt_counts_allCHR.txt')
both <- fread('29Jun20_merged_both_counts_allCHR.txt')
info <- read.delim('29Jun20_merged_sample_info.txt')+

# check that sample/individual IDs match up properly
identical(names(alt[,-1]), names(both[,-1]))
identical(as.character(info$sample_ID), names(both[,-1]))

# create variables to save results
pvals <- matrix(nrow = dim(alt)[1], ncol=6)
coefs <- matrix(nrow = dim(alt)[1], ncol=6)
se <- matrix(nrow = dim(alt)[1], ncol=6)
lik <- matrix(nrow = dim(alt)[1], ncol=2)

# fit model and save pvals, coefs, se, lik at each of the sites
for (i in 1:nrow(alt)) {
  info$tmp_x <- as.vector(t(alt[i,-1]))
  info$tmp_y <- as.vector(t(both[i,-1]))
  # exclude samples with no reads at a given site, and samples that couldn't be sexed
  tmp_info <- subset(info, tmp_y>0 & sex!='unknown')
  mod1 <- tryCatch(betabin(cbind(tmp_x, tmp_y - tmp_x) ~ condition + sequencing_batch + meta_cage + sex, ~1, data=tmp_info), error=function(x){})
  mod2 <- tryCatch(betabin(cbind(tmp_x, tmp_y - tmp_x) ~ sequencing_batch + meta_cage + sex, ~1, data=tmp_info), error=function(x){})
  tryCatch(lik[i,1] <- logLik(mod1), error = function(x){})
  tryCatch(lik[i,2] <- logLik(mod2), error = function(x){})
  tryCatch(pvals[i,] <- attributes(summary(mod1))$Coef[,4], error=function(x){})
  tryCatch(se[i,] <- attributes(summary(mod1))$Coef[,2], error=function(x){})
  tryCatch(coefs[i,] <- attributes(summary(mod1))$Coef[,1], error=function(x){})
}

# summary data from authors to compare pvals and coefs
summary <- read.delim('SummaryTable_allsites_12Nov20.txt')

# vector to save author classifications
author_sig_cat <- rep("NS", nrow(alt))

# look for mismatches
for (i in 1:nrow(alt)) {
  summary_idx <- which(summary$site == alt$site[i], arr.ind = TRUE)
  if (length(summary_idx) != 0 ) {
    # save author_sig_cat
    author_sig_cat[i] <- summary$sig_cat[summary_idx]
    # check if p values mismatch
    if (round(pvals[i,2],3) != round(summary$pval_CTRL[summary_idx],3)) {
      print("mismatching pval_CTRL at:") 
      print(alt$site[i])
      print(i)
    }
    if (round(pvals[i,3],3) != round(summary$pval_HS[summary_idx],3)) {
      print("mismatching pval_HS at:")
      print(alt$site[i])
      print(i)
    }
    # check if coefficients mismatch
    if (round(coefs[i,2],3) != round(summary$coef_CTRL[summary_idx],3)) {
      print("mismatching coef_CTRL at:")
      print(alt$site[i])
      print(i)
    }
    if (round(coefs[i,3],3) != round(summary$coef_HS[summary_idx],3)) {
      print("mismatching coef_HS at:")
      print(alt$site[i])
      print(i)
    }
  }
}
# p values and coefficients matched up to 3 decimal places

# plot regression coefficients (recreate figure 2C in Pallares et al.)
df.coefs <- data.frame(ctrl = coefs[,2],
                       hs = coefs[,3],
                       sigcat = author_sig_cat)
ggplot(df.coefs, aes(x = ctrl, y = hs, color = sigcat, shape = sigcat)) +
  geom_point() + 
  theme_bw() +
  ylab("Genetic effect: HS diet") +
  xlab("Genetic effect: CTRL diet") +
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "lightgray",
                                "HS" = "orange",
                                "CTRL" = "skyblue")) +
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 1,
                                "HS" = 16,
                                "CTRL" = 16)) +
  xlim(-1.5,1) + 
  ylim(-1.5,1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
