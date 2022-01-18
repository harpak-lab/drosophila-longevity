# this script fits the new model with an interaction term for time:env;
# I ran it on TACC (using the associated slurm file) after splitting the 
# count data into 272 files of 1000 rows each

#!/usr/bin/env Rscript

library(aod, lib.loc = '/home1/08397/rkknow/R/x86_64-pc-linux-gnu-library/4.0/')
library(data.table, lib.loc = '/home1/08397/rkknow/R/x86_64-pc-linux-gnu-library/4.0/')
library(dplyr, lib.loc = '/home1/08397/rkknow/R/x86_64-pc-linux-gnu-library/4.0/')

# get which portion of data we're running from slurm file
args <- commandArgs(trailingOnly = TRUE)
file_num <- args[1]

# read in count data
alt <- fread(file = paste0('Data/29Jun20_merged_alt_counts_allCHR_',file_num,'.txt'))
both <- fread(file = paste0('Data/29Jun20_merged_both_counts_allCHR_',file_num,'.txt'))
info <- read.delim(file = 'Data/29Jun20_merged_sample_info.txt')

# create variables to save pvals, coefs, se
pvals <- matrix(nrow=dim(alt)[1], ncol = 6)
coefs <- matrix(nrow=dim(alt)[1], ncol = 6)
se <- matrix(nrow=dim(alt)[1], ncol = 6)

# fit the model and save results
for (i in 1:dim(alt)[1]) {
	info$tmp_x <- as.vector(t(alt[i,-1]))
	info$tmp_y <- as.vector(t(both[i,-1]))

	# exclude samples with no reads at a given site, and samples that couldn't be sexed
	tmp_info <- subset(info, tmp_y>0 & sex!='unknown')
	
	# declare new time variable and env
	tmp_info$time <- case_when(tmp_info$timepoint == "T0" ~ 0,
				   tmp_info$timepoint == "TN" ~ 1)
	tmp_info$env <- case_when(tmp_info$condition == "both" ~ 0,
				  tmp_info$condition == "C" ~ 0,
				  tmp_info$condition == "HS" ~ 1)

	# fit model with interaction term
	int_model <- tryCatch(betabin(cbind(tmp_x, tmp_y - tmp_x) ~ time + env:time + sequencing_batch + meta_cage + sex,
	                              ~1, data = tmp_info), error = function(x){})

	# save results
	tryCatch(pvals[i,] <- attributes(summary(int_model))$Coef[,4], error = function(x){})
	tryCatch(se[i,] <- attributes(summary(int_model))$Coef[,2], error = function(x){})
	tryCatch(coefs[i,] <- attributes(summary(int_model))$Coef[,1], error = function(x){})
}

# write model results to txt files
write.table(pvals, file = paste0("ModelResults/int_model_pvals_",file_num,".txt"), quote = FALSE)
write.table(se, file = paste0("ModelResults/int_model_se_",file_num,".txt"), quote = FALSE)
write.table(coefs, file = paste0("ModelResults/int_model_coefs_",file_num,".txt"), quote = FALSE)
