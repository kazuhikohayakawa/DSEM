# Multilevel Bayes

# I. load data
rm(list=ls())   # clear all variables
cat("\014")     # clear console

# getwd()
# setwd("C:/Users/econ02/Desktop/real_data")

# sink("Estimation_results_all_2024_1111.txt",split =TRUE)

library(cmdstanr)
library(coda)
library(VAR.etp)
library(brms)
library(freqdom)
library(brms)


Data <-read.csv("Data_McNeish.csv", header = T)
str(Data)

T <- 50
N <- 100

y_TN <- Data$Urge
x_TN <- Data$Dep
id   <- matrix(1:N,nrow=N,ncol=1)
Job  <- Data$Job     
Home <- Data$Home 
indices<-seq_along(Job) %%50 ==0
Job <- Job[indices]
Home <- Home[indices]


T1     <- T-1
Y_TN   <- matrix(y_TN, T, N)
X_TN   <- matrix(x_TN, T, N)

#centering
y       <- matrix(Y_TN[2:T,],  T1*N, 1)
y_lag1  <- matrix(Y_TN[1:T1,], T1*N, 1)
x       <- matrix(X_TN[2:T,],  T1*N, 1)
ID      <- kronecker(id,matrix(1,nrow=T1,ncol=1))

Y_lag1_c   <- sweep(Y_TN[1:T1,], 2, apply(Y_TN[1:T1,],2,mean))
y_lag1_c   <-  matrix(Y_lag1_c, T1*N, 1)
X_TN_c    <- sweep( X_TN[2:T,], 2, apply(X_TN[2:T,],2,mean) )
x_c        <- matrix(X_TN_c, T1*N, 1)
  
Job_c  <- Job - mean(Job)  # N*1 vector
Home_c <- Home - mean(Home)   # N*1 vector
# Dep_b  <- apply(X_TN[2:T,],2,mean)    # N*1 vector
Dep_b <- apply(X_TN, 2, mean)

Job_C     <- matrix(1, nrow = T1, ncol=1)%*%matrix(Job_c,1,N)
Job_C     <- matrix(Job_C, T1*N,1)   # (N*T1)*1 vector
Home_C    <- matrix(1, nrow = T1, ncol=1)%*%matrix(Home_c,1,N)
Home_C    <- matrix(Home_C, T1*N,1)  # (N*T1)*1 vector
Dep_B     <- matrix(1, nrow = T1, ncol=1)%*%matrix(Dep_b,1,N)
Dep_B     <- matrix(Dep_b, T1*N,1)   # (N*T1)*1 vector

options(warn = -1)


##########################################
source("proc_MG_real.R")
source("proc_lmer_real.R")
source("proc_nmle_real.R")
source("proc_brm_real.R")
source("proc_stan_real.R")
source("proc_reg_OLS_JK_BC_real.R")

# n_ite    = 2000
# n_burnin = 500
n_ite    = 500
n_burnin = 100

name_head <- c("Estimate",  "CI95.lower", "CI95.upper")
name_para <- c("alpha", "phi", "beta", "tau","gamma_alpha1", "gamma_alpha2", "gamma_alpha3", "gamma_phi1", "gamma_phi2", "gamma_beta1", "gamma_beta2", "gamma_tau1", "gamma_tau2", "ome2_alpha", "ome2_phi", "ome2_beta", "ome2_tau")


#########################################################

MG_results      <- proc_MG_real(Y_TN, X_TN, Job_c, Home_c, Dep_b, ID, N )
REML1_results   <- proc_lmer_real(y, y_lag1_c, x_c, Job_C, Home_C, Dep_B, ID, N, T1 )
REML2_results   <- proc_nmle_real(y, y_lag1_c, x_c, Job_C, Home_C, Dep_B, Job_c, Home_c, Dep_b, ID, N, T1 )
Bayes1_results  <- proc_brm_real(y, y_lag1_c, x_c, Job_C, Home_C, Dep_B, ID, n_ite, n_burnin )
Bayes2_results  <- proc_stan_real(y, y_lag1_c, x_c, Job_c, Home_c, Dep_b, ID, N, T1, n_ite, n_burnin )

##################################################

theta_gamma_est_MG    <- MG_results[[1]]
theta_gamma_CI_MG     <- MG_results[[2]]
theta_gamma_est_JKMG  <- MG_results[[3]]
theta_gamma_CI_JKMG   <- MG_results[[4]]
theta_gamma_est_BCMG1 <- MG_results[[5]]
theta_gamma_CI_BCMG1  <- MG_results[[6]]
theta_gamma_est_BCMG2 <- MG_results[[7]]
theta_gamma_CI_BCMG2  <- MG_results[[8]]

# omega_MG_est        <- MG_results[[9]]
# omega_JKMG_est      <- MG_results[[10]]
# omega_BCMG1_est     <- MG_results[[11]]
# omega_BCMG2_est     <- MG_results[[12]]
omega_bc_MG_est     <- MG_results[[13]]
omega_bc_JKMG_est   <- MG_results[[14]]
omega_bc_BCMG1_est  <- MG_results[[15]]
omega_bc_BCMG2_est  <- MG_results[[16]]
time_MGs            <- MG_results[[17]]



theta_gamma_est_REML1    <- REML1_results[[1]]
theta_gamma_CI_REML1     <- REML1_results[[2]]
omega_REML1_est          <- REML1_results[[3]]
time_REML1               <- REML1_results[[4]]


theta_gamma_est_REML2    <- REML2_results[[1]]
theta_gamma_CI_REML2     <- REML2_results[[2]]
omega_REML2_est          <- REML2_results[[3]]
time_REML2               <- REML2_results[[4]]


theta_gamma_est_Bayes1    <- Bayes1_results[[1]]
theta_gamma_se_Bayes1     <- Bayes1_results[[2]]
theta_gamma_Bayes1_CI_L   <- Bayes1_results[[3]]
theta_gamma_Bayes1_CI_U   <- Bayes1_results[[4]]
Bayes1_Rhat               <- Bayes1_results[[5]]
omega_Bayes1_est          <- Bayes1_results[[6]]
omega_Bayes1_CI_L         <- Bayes1_results[[7]]
omega_Bayes1_CI_U         <- Bayes1_results[[8]]
time_Bayes1               <- Bayes1_results[[9]]


fixed_effect_Bayes2       <- Bayes2_results[[1]]
omega_Bayes2              <- Bayes2_results[[2]]
time_Bayes2               <- Bayes2_results[[3]]

theta_gamma_est_Bayes2    <- as.matrix(fixed_effect_Bayes2[,1])
theta_gamma_se_Bayes2     <- as.matrix(fixed_effect_Bayes2[,2])
theta_gamma_Bayes2_CI_L   <- as.matrix(fixed_effect_Bayes2[,3])
theta_gamma_Bayes2_CI_U   <- as.matrix(fixed_effect_Bayes2[,4])
Bayes2_Rhat               <- as.matrix(fixed_effect_Bayes2[,5])
omega_Bayes2_est          <- as.matrix(omega_Bayes2[,1])
omega_Bayes2_CI_L         <- as.matrix(omega_Bayes2[,3])
omega_Bayes2_CI_U         <- as.matrix(omega_Bayes2[,4])


########################################################################
### MG ###
theta_gamma_est_CI_MG    <- cbind(theta_gamma_est_MG,    theta_gamma_CI_MG)
theta_gamma_est_CI_JKMG  <- cbind(theta_gamma_est_JKMG,  theta_gamma_CI_JKMG)
theta_gamma_est_CI_BCMG1 <- cbind(theta_gamma_est_BCMG1, theta_gamma_CI_BCMG1)
theta_gamma_est_CI_BCMG2 <- cbind(theta_gamma_est_BCMG2, theta_gamma_CI_BCMG2)

omega_est_CI_MG    <- cbind(omega_bc_MG_est,    matrix(99,4,2))
omega_est_CI_JKMG  <- cbind(omega_bc_JKMG_est,  matrix(99,4,2))
omega_est_CI_BCMG1 <- cbind(omega_bc_BCMG1_est, matrix(99,4,2))
omega_est_CI_BCMG2 <- cbind(omega_bc_BCMG2_est, matrix(99,4,2))

theta_gamma_omega_est_CI_MG    <- rbind(theta_gamma_est_CI_MG,    omega_est_CI_MG)
theta_gamma_omega_est_CI_JKMG  <- rbind(theta_gamma_est_CI_JKMG,  omega_est_CI_JKMG)
theta_gamma_omega_est_CI_BCMG1 <- rbind(theta_gamma_est_CI_BCMG1, omega_est_CI_BCMG1)
theta_gamma_omega_est_CI_BCMG2 <- rbind(theta_gamma_est_CI_BCMG2, omega_est_CI_BCMG2)

colnames(theta_gamma_omega_est_CI_MG)    <- name_head
colnames(theta_gamma_omega_est_CI_JKMG)  <- name_head
colnames(theta_gamma_omega_est_CI_BCMG1) <- name_head
colnames(theta_gamma_omega_est_CI_BCMG2) <- name_head

rownames(theta_gamma_omega_est_CI_MG)    <- name_para
rownames(theta_gamma_omega_est_CI_JKMG)  <- name_para
rownames(theta_gamma_omega_est_CI_BCMG1) <- name_para
rownames(theta_gamma_omega_est_CI_BCMG2) <- name_para


########################################################################
### REML1 ###
theta_gamma_est_CI_REML1 <- cbind(theta_gamma_est_REML1, theta_gamma_CI_REML1)
omega_est_CI_REML1       <- cbind(omega_REML1_est,    matrix(99,4,2))

theta_gamma_omega_est_CI_REML1  <- rbind(theta_gamma_est_CI_REML1,  omega_est_CI_REML1)

colnames(theta_gamma_omega_est_CI_REML1)    <- name_head
rownames(theta_gamma_omega_est_CI_REML1)    <- name_para

########################################################################
### REML2 ###
theta_gamma_est_CI_REML2 <- cbind(theta_gamma_est_REML2, theta_gamma_CI_REML2)
omega_est_CI_REML2       <- cbind(omega_REML2_est,    matrix(99,4,2))

theta_gamma_omega_est_CI_REML2  <- rbind(theta_gamma_est_CI_REML2, omega_est_CI_REML2)

colnames(theta_gamma_omega_est_CI_REML2)    <- name_head
rownames(theta_gamma_omega_est_CI_REML2)    <- name_para

########################################################################
### Bayes ###

theta_gamma_omega_est_Bayes1   <- rbind(theta_gamma_est_Bayes1, omega_Bayes1_est)
theta_gamma_omega_est_Bayes2   <- rbind(theta_gamma_est_Bayes2, omega_Bayes2_est)

theta_gamma_omega_Bayes1_CI_L  <- rbind( theta_gamma_Bayes1_CI_L, omega_Bayes1_CI_L)
theta_gamma_omega_Bayes1_CI_U  <- rbind( theta_gamma_Bayes1_CI_U, omega_Bayes1_CI_U)
theta_gamma_omega_Bayes2_CI_L  <- rbind( theta_gamma_Bayes2_CI_L, omega_Bayes2_CI_L)
theta_gamma_omega_Bayes2_CI_U  <- rbind( theta_gamma_Bayes2_CI_U, omega_Bayes2_CI_U)

theta_gamma_omega_est_CI_Bayes1 <- cbind( theta_gamma_omega_est_Bayes1,  theta_gamma_omega_Bayes1_CI_L, theta_gamma_omega_Bayes1_CI_U, rbind(Bayes1_Rhat, matrix(99,4,1) ) )
theta_gamma_omega_est_CI_Bayes2 <- cbind( theta_gamma_omega_est_Bayes2,  theta_gamma_omega_Bayes2_CI_L, theta_gamma_omega_Bayes2_CI_U, rbind(Bayes2_Rhat, matrix(99,4,1) ) )

colnames(theta_gamma_omega_est_CI_Bayes1)    <- c(name_head, "Rhat")
rownames(theta_gamma_omega_est_CI_Bayes1)    <- name_para

colnames(theta_gamma_omega_est_CI_Bayes2)    <- c(name_head, "Rhat")
rownames(theta_gamma_omega_est_CI_Bayes2)    <- name_para



cat(" \n **********************************************************")
cat(" \n **********************************************************")
cat(" \n  Estimation results of the real data model:  \n")
cat("  \n") 
cat(" MG  \n")
print(round(theta_gamma_omega_est_CI_MG, 4))
cat("  \n")
cat(" JKMG \n")
print(round(theta_gamma_omega_est_CI_JKMG, 4))
cat("  \n")
cat(" BCMG1  \n")
print(round(theta_gamma_omega_est_CI_BCMG1, 4))
cat("  \n")
cat(" BCMG2  \n")
print(round(theta_gamma_omega_est_CI_BCMG2, 4))
cat("  \n")
cat(" REML1  \n")
print(round(theta_gamma_omega_est_CI_REML1, 4))
cat("  \n")
cat(" REML2  \n")
print(round(theta_gamma_omega_est_CI_REML2, 4))
cat("  \n")
cat(" Bayes1  \n")
print(round(theta_gamma_omega_est_CI_Bayes1, 4))
cat("  \n")
cat(" Bayes2 \n")
print(round(theta_gamma_omega_est_CI_Bayes2, 4))
cat("  \n")

########################
comp_time <- cbind(time_MGs, time_REML1, time_REML2, time_Bayes1, time_Bayes2)
colnames(comp_time) <- c("MG","REML1", "REML2", "Bayes(brm)", "Bayes(stan)")
cat("**************************************************************************************************************  \n")
cat(" comp_time   \n")
print(comp_time)
cat("  \n")
cat("**************************************************************************************************************  \n")




# sink()
