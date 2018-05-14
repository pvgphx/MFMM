main = "C:/Users/bingsi/Dropbox/MyWorkNew/2016_12_26_migraine_data_covariate_Modality_Selection_new_idea_2018_5_11_JASA" 
setwd(main)

#*******************************************
# Install a few R packages needed in the code
#*******************************************
install.packages("MASS")
install.packages("mvtnorm")
install.packages("Matrix")
install.packages("gglasso")
install.packages("psych")

#********************************************
# Initialization parameters
#********************************************
# Number of modes
m = 4 
# Number of informative features
p0 = 10; 
# Number of subjects in subtype 1
num0 = 48; 
# Number of subjects, including subtypes 1 and 2
n = 120; 
# Number of experiments
replicate = 20 

# tuning parameter range for number of clusters
k_seq = 1:3
# tuning parameter range for number of factors 
r_seq = 1:3
# tuning parameter range for beta
beta_seq = seq(5, 15, 1)
# tuning parameter range for gamma
gamma_seq = seq(1500,3500, 100)	


#*********************************************************************************
# Fit MFMM for each of the 20 runs
# Compute mode and feature selection accuracy for each run
#*********************************************************************************
mode_our = vector(length = replicate)
specificity_H_our_each = matrix(nrow = 4, ncol = replicate)
sensitivity_H_our_each = matrix(nrow = 4, ncol = replicate)

for(num_r in 1:replicate){
	setwd(paste(main, "/data/", num_r, sep = ""))	
	source("main_mfmm.R")
	
	#*******************************************************************
	# Load the simulated data 
	# Load the  
	#*******************************************************************
	
	load("best.RData")
	load("X.RData")
	load("Z.RData")
	#print(best)

	lambda_k = k_seq[best[1]]
	lambda_r = r_seq[best[2]]
	lambda_beta = beta_seq[best[3]]
	lambda_gamma = gamma_seq[best[4]]

	load("fit.RData")
	#it = 50; eps = 0.0001;scaling = FALSE; init = NULL; 
	#fit = mfmm(X,Z,lambda_k, r_list = rep(lambda_r, m), lambda_beta, lambda_gamma, x.z = NULL, x.w = NULL, it = it, eps = 0.0001, seed = 4, scaling = FALSE)
	#print(fit$beta)
	#save(fit, file = "fit.RData")
	
	#*******************************************************************
	# Facilitate the computation of mode and feature selection accuracy 
	#*******************************************************************
	mode_our[num_r] = best[1]
	for(num_m in 1:m){specificity_H_our_each[num_m, num_r] = sum(fit$H[1:p0,,num_m] != 0)/(p0 * lambda_r)}
	for(num_m in 1:m){sensitivity_H_our_each[num_m, num_r] = sum(fit$H[(p0 + 1):(4 * p0),,num_m] == 0)/(3 *p0 * lambda_r)}	
}

#*******************************************************
# Compute averaged mode selection accuracy for 20 runs
#*******************************************************
sum(mode_our != 1)/replicate

#**********************************************************
# Compute averaged feature selection accuracy for 20 runs
# Reproduce Table 2 in the paper
#******************************************8***************
Table_2 = matrix(nrow = 2, ncol = 4)
rownames(Table_2) = c("Feature selection sensitivity(%)", "Feature selection specificity(%)")
colnames(Table_2) = c("Mode 1", "Mode 2", "Mode 3", "Mode 4")
for(num_m in 1:4){
	Table_2[1, num_m] = paste(round(mean(specificity_H_our_each[num_m,])*100, digits = 1), "+/-", round(sd(specificity_H_our_each[num_m, ])*100, digits = 1), sep = "") 
	Table_2[2, num_m] = paste(round(mean(sensitivity_H_our_each[num_m, ])*100, digits = 1), "+/-", round(sd(sensitivity_H_our_each[num_m, ])*100, digits = 1), sep = "") 
}
print(Table_2)





