#***********************************************************
# It takes hours to search through the space and 
# find the optimal tuning parameters with minimal BIC. 
# So we conducted the search first and saved the optimal
# parameters we found by fitting MFMM across the space
# of the parameters. 
# 
# To show that the optimal tuning 
# parameters we provided can be reproduced from the 
# code, we take the 3rd dataset as an example.  
# 
# First, we load the optimal parameters from saved
# object "best.RData",fit MFMM with these parameters
# and print the BIC (-24941.14) of the fitted model. 
# 
# Then, we encourage reviewers to try another combination 
# of parameters to double check if there is any other 
# combination resulting a BIC lower than -24941.14. 
# If not, it supports that the parameters saved in 
# "best.RData" is valid and can be reproduced. 
#**************************************************************

#*******************************************
# Install a few R packages needed in the code
#*******************************************
install.packages("MASS")
install.packages("mvtnorm")
install.packages("Matrix")
install.packages("gglasso")
install.packages("psych")

#********************************************
# Range for tuning parameters
#********************************************
# tuning parameter range for number of clusters
k_seq = 1:3
# tuning parameter range for number of factors 
r_seq = 1:3
# tuning parameter range for beta
beta_seq = seq(5, 15, 1)
# tuning parameter range for gamma
gamma_seq = seq(1500,3500, 100)	
# number of modes
m = 4

#*****************************************************************
# Load the 3rd dataset as an example
# Feel free to try other datasets
#******************************************************************
num_r = 3
setwd(paste(main, "/data/", num_r, sep = ""))	
source("main_mfmm.R")	
load("best.RData")
load("X.RData")
load("Z.RData")

#*****************************************************************
# Load the optimal parameters with minimal BIC from "best.RData"
#******************************************************************
lambda_k = k_seq[best[1]]
lambda_r = r_seq[best[2]]
lambda_beta = beta_seq[best[3]]
lambda_gamma = gamma_seq[best[4]]
optimal_parameters = c(lambda_k, lambda_r, lambda_beta, lambda_gamma)
names(optimal_parameters) = c("#clusters", "#factors", "lambda1", "lambda2")
print(optimal_parameters)
# Fit a MFMM with the optimal tuning parameters
it = 50; eps = 0.0001;scaling = FALSE; init = NULL; 
fit = mfmm(X,Z,lambda_k, r_list = rep(lambda_r, m), lambda_beta, lambda_gamma, x.z = NULL, x.w = NULL, it = it, eps = 0.0001, seed = 4, scaling = FALSE)
print(fit$bic)

#*****************************************************************
# Try another combination of the tuning parameters
#*****************************************************************
#k_seq = 1:3
#r_seq = 1:3
#beta_seq = seq(5, 15, 1)
#gamma_seq = seq(1500,3500, 100)	
lambda_k = 1
lambda_r = 2
lambda_beta = 15
lambda_gamma = 3000
other_parameters = c(lambda_k, lambda_r, lambda_beta, lambda_gamma)
names(other_parameters) = c("#clusters", "#factors", "lambda1", "lambda2")
print(other_parameters)
# Fit a MFMM with other parameters
it = 50; eps = 0.0001;scaling = FALSE; init = NULL; 
fit = mfmm(X,Z,lambda_k, r_list = rep(lambda_r, m), lambda_beta, lambda_gamma, x.z = NULL, x.w = NULL, it = it, eps = 0.0001, seed = 4, scaling = FALSE)
print(fit$bic)


