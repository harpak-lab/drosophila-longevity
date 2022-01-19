# this script reads in the TACC output for the interaction model and creates
# several plots comparing these results to the original author classifications

library(data.table)
library(ggplot2)
library(gridExtra)

# read in coefs, pvals, se from the 272 output files
int_model_coefs <- matrix(nrow = 271536, ncol = 6)
int_model_pvals <- matrix(nrow = 271536, ncol = 6)
int_model_se <- matrix(nrow = 271536, ncol = 6)
for (i in 1:272) {
  coefs_tmp <- as.matrix(suppressWarnings(fread(input = paste0("int_model_coefs_",i,".txt"))[,2:7]))
  pvals_tmp <- as.matrix(suppressWarnings(fread(input = paste0("int_model_pvals_",i,".txt"))[,2:7]))
  se_tmp <- as.matrix(suppressWarnings(fread(input = paste0("int_model_se_",i,".txt"))[,2:7]))
  int_model_coefs[(1000 * (i - 1) + 1):(1000 * (i - 1) + nrow(coefs_tmp)),] <- coefs_tmp
  int_model_pvals[(1000 * (i - 1) + 1):(1000 * (i - 1) + nrow(coefs_tmp)),] <- pvals_tmp
  int_model_se[(1000 * (i - 1) + 1):(1000 * (i - 1) + nrow(coefs_tmp)),] <- se_tmp
}

# save full version (all SNPs) of output files
write.table(int_model_coefs, file = "ALL_int_model_coefs.txt")
write.table(int_model_pvals, file = "ALL_int_model_pvals.txt")
write.table(int_model_se, file = "ALL_int_model_se.txt")

# look at -log10 of pvals from the interaction model
my_pvals <- -log10(int_model_pvals[,6])
my_pvals[which(is.na(my_pvals))] <- 0
min(my_pvals)
max(my_pvals)

# iterate through 10,000 possible pval significance thresholds 
pval_threshold <- seq(from = 2, to = 15.2, length = 10000)
# for each pval threshold, count how many SNPs had a significant interaction term
num_int_sig <- rep(NA, 10000)
for (i in 1:length(num_int_sig)) {
  num_int_sig[i] <- length(which(my_pvals > pval_threshold[i]))
}
plot(x = pval_threshold, y = num_int_sig)

# read in author classifications
author_classifications <- rep(NA, length = nrow(alt))
summary_info <- read.delim('SummaryTable_allsites_12Nov20.txt')
alt <- fread('29Jun20_merged_alt_counts_allCHR.txt')
for (i in 1:nrow(alt)) {
  summary_idx <- which(summary_info$site == alt$site[i], arr.ind = TRUE)
  if (length(summary_idx) != 0 ) {
    author_classifications[i] <- summary_info$sig_cat[summary_idx]  
  }
}

# find how many of the significant interaction terms at each pval threshold
# were originally classified as HS, both, or CTRL
num_HS_sig <- rep(0, length(pval_threshold))
num_both_sig <- rep(0, length(pval_threshold))
num_CTRL_sig <- rep(0, length(pval_threshold))
for (i in 1:length(pval_threshold)) {
  num_HS_sig[i] = length(which(author_classifications[which(my_pvals > pval_threshold[i])] == "HS"))
  num_both_sig[i] = length(which(author_classifications[which(my_pvals > pval_threshold[i])] == "both"))
  num_both_sig[i] = length(which(author_classifications[which(my_pvals > pval_threshold[i])] == "CTRL"))
}

# what proportion of the significant interaction terms fell into each category?
prop_HS_sig <- num_HS_sig / num_int_sig
prop_both_sig <- num_both_sig / num_int_sig
prop_CTRL_sig <- num_CTRL_sig / num_int_sig

# create dataframe for plotting
plot_data = data.frame(pval_threshold = pval_threshold,
                       num_int_sig = num_int_sig,
                       num_HS_sig = num_HS_sig,
                       num_both_sig = num_both_sig,
                       num_CTRL_sig = num_CTRL_sig,
                       prop_HS_sig = prop_HS_sig,
                       prop_both_sig = prop_both_sig,
                       prop_CTRL_sig = prop_CTRL_sig)

# plot proportion classified as CTRL, HS, or shared at each significance threshold
ggplot(data = plot_data) +
  ylab("Proportion") +
  xlab("-log10(p) Significance Threshold")+
  geom_line(aes(x = pval_threshold, y = prop_HS_sig, color = "HS"), lwd = 1.5) +
  geom_line(aes(x = pval_threshold, y = prop_both_sig, color = "shared"), lwd = 1.5) +
  geom_line(aes(x = pval_threshold, y = prop_CTRL_sig, color = "CTRL"), lwd = 1.5) + 
  theme_bw() +
  scale_color_manual(name = "Author Classification", 
                     values = c("HS" = "cornflowerblue", 
                                "shared" = "deeppink2",
                                "CTRL" = "goldenrod1"))

# what do the interaction model p values look like for sites the authors classified as significant?
log10_int_pvals <- -log10(int_model_pvals[,6])
log10_time_pvals <- -log10(int_model_pvals[,2])
log10_int_pvals[which(is.na(log10_int_pvals))] <- 0
log10_time_pvals[which(is.na(log10_time_pvals))] <- 0
plot_data2 <- data.frame(log10_int_pvals = log10_int_pvals,
                         log10_time_pvals = log10_time_pvals,
                         author_classifications = author_classifications)

ggplot(data = plot_data2, aes(x = log10_int_pvals, y = log10_time_pvals, 
                              color = author_classifications, shape = author_classifications)) +
  geom_point() + 
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "red",
                                "HS" = "purple",
                                "CTRL" = "green3")) +
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 16,
                                "HS" = 16,
                                "CTRL" = 16)) +
  theme_bw() 

# want to change order of levels for author_classifications so that NS is plotted first
plot_data3 <- data.frame(log10_int_pvals = log10_int_pvals,
                         log10_time_pvals = log10_time_pvals,
                         author_classifications = author_classifications)

plot_data3$author_classifications[which(plot_data3$author_classifications == "NS")] <- "a NS"
plot_data3$author_classifications[which(plot_data3$author_classifications == "shared")] <- "b shared"
plot_data3$author_classifications[which(plot_data3$author_classifications == "HS")] <- "c HS"
plot_data3$author_classifications[which(plot_data3$author_classifications == "CTRL")] <- "d CTRL"

plot_data3 <- plot_data3[order(plot_data3$author_classifications),]

plot_data3$author_classifications[which(plot_data3$author_classifications == "a NS")] <- "NS"
plot_data3$author_classifications[which(plot_data3$author_classifications == "b shared")] <- "shared"
plot_data3$author_classifications[which(plot_data3$author_classifications == "c HS")] <- "HS"
plot_data3$author_classifications[which(plot_data3$author_classifications == "d CTRL")] <- "CTRL"


ggplot(data = plot_data3, aes(x = log10_int_pvals, y = log10_time_pvals, 
                              color = author_classifications, shape = author_classifications)) +
  geom_point(alpha = 0.35) + 
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "deeppink2",
                                "HS" = "cornflowerblue",
                                "CTRL" = "goldenrod2")) +
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 16,
                                "HS" = 16,
                                "CTRL" = 16)) +
  theme_bw() +
  ylab("-log10(p) [time]") +
  xlab("-log10(p) [time:environment]") +
  labs(color = "Author Classification", shape = "Author Classification") 

# plot coefficients
int_coefs <- int_model_coefs[,6]
time_coefs <- int_model_coefs[,2]
plot_data4 <- data.frame(int_coefs = int_coefs,
                         time_coefs = time_coefs,
                         author_classifications = author_classifications)

plot_data4$author_classifications[which(plot_data4$author_classifications == "NS")] <- "a NS"
plot_data4$author_classifications[which(plot_data4$author_classifications == "shared")] <- "b shared"
plot_data4$author_classifications[which(plot_data4$author_classifications == "HS")] <- "c HS"
plot_data4$author_classifications[which(plot_data4$author_classifications == "CTRL")] <- "d CTRL"

plot_data4 <- plot_data4[order(plot_data4$author_classifications),]

plot_data4$author_classifications[which(plot_data4$author_classifications == "a NS")] <- "NS"
plot_data4$author_classifications[which(plot_data4$author_classifications == "b shared")] <- "shared"
plot_data4$author_classifications[which(plot_data4$author_classifications == "c HS")] <- "HS"
plot_data4$author_classifications[which(plot_data4$author_classifications == "d CTRL")] <- "CTRL"


ggplot(data = plot_data4, aes(x = int_coefs, y = time_coefs, 
                              color = author_classifications, shape = author_classifications)) +
  geom_point(alpha = 0.35) + 
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "deeppink2",
                                "HS" = "cornflowerblue",
                                "CTRL" = "goldenrod2")) +
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 16,
                                "HS" = 16,
                                "CTRL" = 16)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  xlab("time:environment coefficient") +
  ylab("time coefficient")+
  labs(color = "Author Classification", shape = "Author Classification") 


# look at int_coef + time_coef for the HS sites
plot_data5 <- data.frame(author_classifications = author_classifications,
                         idx = 1:length(author_classifications),
                         time_coefs = time_coefs,
                         int_coefs = int_coefs,
                         sum_coefs = int_coefs + time_coefs)

plot_data5$author_classifications[which(plot_data5$author_classifications == "NS")] <- "a NS"
plot_data5$author_classifications[which(plot_data5$author_classifications == "shared")] <- "b shared"
plot_data5$author_classifications[which(plot_data5$author_classifications == "HS")] <- "c HS"
plot_data5$author_classifications[which(plot_data5$author_classifications == "CTRL")] <- "d CTRL"

plot_data5 <- plot_data5[order(plot_data5$author_classifications),]

plot_data5$author_classifications[which(plot_data5$author_classifications == "a NS")] <- "NS"
plot_data5$author_classifications[which(plot_data5$author_classifications == "b shared")] <- "shared"
plot_data5$author_classifications[which(plot_data5$author_classifications == "c HS")] <- "HS"
plot_data5$author_classifications[which(plot_data5$author_classifications == "d CTRL")] <- "CTRL"

ggplot(plot_data5, aes(x = time_coefs, y = sum_coefs, color = author_classifications, shape = author_classifications)) +
  geom_point(alpha = 0.35) +
  theme_bw() +
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "deeppink2",
                                "HS" = "cornflowerblue",
                                "CTRL" = "goldenrod2")) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 16,
                                "HS" = 16,
                                "CTRL" = 16)) +
  labs(color = "Author Classification", shape = "Author Classification")

# want to split plot into version where they had same signs, and where they had opposite signs
plot_data5_samesigns <- plot_data5[which(plot_data5$time_coefs * plot_data5$int_coefs > 0),]
plot_data5_diffsigns <- plot_data5[which(plot_data5$time_coefs * plot_data5$int_coefs < 0),]

ggplot(plot_data5_samesigns, aes(x = idx, y = sum_coefs, color = author_classifications, shape = author_classifications)) +
  geom_point(alpha = 0.35) +
  theme_bw() +
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "deeppink2",
                                "HS" = "cornflowerblue",
                                "CTRL" = "goldenrod2")) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 16,
                                "HS" = 16,
                                "CTRL" = 16)) +
  labs(color = "Author Classification", shape = "Author Classification") +
  ylab("beta_t + beta_te") +
  xlab("") +
  ggtitle("Same Signs") +
  theme(axis.text.x = element_blank())


ggplot(plot_data5_diffsigns, aes(x = idx, y = sum_coefs, color = author_classifications, shape = author_classifications)) +
  geom_point(alpha = 0.35) +
  theme_bw() +
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "deeppink2",
                                "HS" = "cornflowerblue",
                                "CTRL" = "goldenrod2")) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 16,
                                "HS" = 16,
                                "CTRL" = 16)) +
  labs(color = "Author Classification", shape = "Author Classification") +
  ylab("beta_t + beta_te") +
  xlab("") +
  ggtitle("Opposite Signs") +
  theme(axis.text.x = element_blank())

# z scores for coefs
int_z_scores <- int_model_coefs[,6] / int_model_se[,6]

time_z_scores <- int_model_coefs[,2] / int_model_se[,2]

mean(int_z_scores)
int_z_scores[which(is.na(int_z_scores))] <- 0

mean(time_z_scores)
time_z_scores[which(is.na(time_z_scores))] <- 0

plot_data6 <- data.frame(int_z_scores = int_z_scores,
                         time_z_scores = time_z_scores,
                         author_classifications = author_classifications)

plot_data6$author_classifications[which(plot_data6$author_classifications == "NS")] <- "a NS"
plot_data6$author_classifications[which(plot_data6$author_classifications == "shared")] <- "b shared"
plot_data6$author_classifications[which(plot_data6$author_classifications == "HS")] <- "c HS"
plot_data6$author_classifications[which(plot_data6$author_classifications == "CTRL")] <- "d CTRL"

plot_data6 <- plot_data6[order(plot_data6$author_classifications),]

plot_data6$author_classifications[which(plot_data6$author_classifications == "a NS")] <- "NS"
plot_data6$author_classifications[which(plot_data6$author_classifications == "b shared")] <- "shared"
plot_data6$author_classifications[which(plot_data6$author_classifications == "c HS")] <- "HS"
plot_data6$author_classifications[which(plot_data6$author_classifications == "d CTRL")] <- "CTRL"


ggplot(data = plot_data6, aes(x = int_z_scores, y = time_z_scores, 
                              color = author_classifications, shape = author_classifications)) +
  geom_point(alpha = 0.35) + 
  scale_color_manual(values = c("NS" = "lightgray",
                                "shared" = "deeppink2",
                                "HS" = "cornflowerblue",
                                "CTRL" = "goldenrod2")) +
  scale_shape_manual(values = c("NS" = 1,
                                "shared" = 16,
                                "HS" = 16,
                                "CTRL" = 16)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(color = "Author Classification", shape = "Author Classification") +
  xlab("time:environment coef/se") +
  ylab("time coef/se")

