# proc_reg_lmer_nmle_ARX1
proc_lmer_real <- function (y, y_lag1, x, Job_C, Home_C, Dep_B, ID, N, T ) { 

  
  time_lmer1 <- Sys.time()
  
  library(VAR.etp)
  library(matrixcalc)  
  library(MASS)
  library(sandwich)
  library(e1071)
  
  library(lme4)
  library(nlme)
  
  dim_theta <-4
  dim_gamma <-9
  
  # Job_C     <- matrix(1, nrow = T, ncol=1)%*%matrix(Job_C,1,N)
  # Job_C     <- matrix(Job_C, T*N,1)
  # Home_C    <- matrix(1, nrow = T, ncol=1)%*%matrix(Home_C,1,N)
  # Home_C    <- matrix(Home_C, T*N,1)
  # Dep_b     <- matrix(1, nrow = T, ncol=1)%*%matrix(Dep_b,1,N)
  # Dep_b     <- matrix(Dep_b, T*N,1)
  
  # suppressWarnings()
  
  ####### lmer-wg
  reg_lmer  <- lmer( y ~ (1 + y_lag1 + x) + ( Job_C + Home_C + Dep_B) + ( y_lag1:Job_C + y_lag1:Home_C )  +( x:Job_C + x:Home_C )  + (  1 + y_lag1 + x | ID))
  coef_lmer <- summary(reg_lmer)$coefficients
  # Vcov_lmer <- as.matrix(summary(reg_lmer)$vcov)  # variance and covariance of fixed effect 
  RE_lmer   <- data.frame(VarCorr(reg_lmer))  # random effect
  
  
  alpha_lmer_est <- coef_lmer[1,1]
  alpha_lmer_se  <- coef_lmer[1,2]
  phi_lmer_est   <- coef_lmer[2,1]
  phi_lmer_se    <- coef_lmer[2,2]
  beta_lmer_est  <- coef_lmer[3,1]
  beta_lmer_se   <- coef_lmer[3,2]
  tau2_lmer_est  <- log(RE_lmer[7,4])
  tau2_lmer_se   <- 99.99
  
  alpha_lmer_Job_est     <- coef_lmer[4,1]
  alpha_lmer_Job_se      <- coef_lmer[4,2]
  alpha_lmer_Home_est    <- coef_lmer[5,1]
  alpha_lmer_Home_se     <- coef_lmer[5,2]  
  alpha_lmer_Dep_est     <- coef_lmer[6,1]
  alpha_lmer_Dep_se      <- coef_lmer[6,2]  
  
  phi_lmer_Job_est    <- coef_lmer[7,1]
  phi_lmer_Job_se     <- coef_lmer[7,2]  
  phi_lmer_Home_est   <- coef_lmer[8,1]
  phi_lmer_Home_se    <- coef_lmer[8,2]
  
  beta_lmer_Job_est    <- coef_lmer[9,1]
  beta_lmer_Job_se     <- coef_lmer[9,2]  
  beta_lmer_Home_est   <- coef_lmer[10,1]
  beta_lmer_Home_se    <- coef_lmer[10,2]
  
  tau2_lmer_Job_est    <- 99.99
  tau2_lmer_Job_se     <- 999.99 
  tau2_lmer_Home_est   <- 99.99
  tau2_lmer_Home_se    <- 999.99
  

  theta_gam_lmer_est    <- matrix(c(alpha_lmer_est, phi_lmer_est, beta_lmer_est, tau2_lmer_est,alpha_lmer_Job_est, alpha_lmer_Home_est,alpha_lmer_Dep_est, phi_lmer_Job_est, phi_lmer_Home_est,beta_lmer_Job_est, beta_lmer_Home_est, tau2_lmer_Job_est, tau2_lmer_Home_est), dim_theta+dim_gamma,1)
  theta_gam_lmer_se     <- matrix(c(alpha_lmer_se, phi_lmer_se, beta_lmer_se, tau2_lmer_se,alpha_lmer_Job_se, alpha_lmer_Home_se, alpha_lmer_Dep_se, phi_lmer_Job_se, phi_lmer_Home_se, beta_lmer_Job_se, beta_lmer_Home_se, tau2_lmer_Job_se, tau2_lmer_Home_se), dim_theta+dim_gamma,1)
  omega_lmer            <- matrix( c(RE_lmer[1,4], RE_lmer[2,4], RE_lmer[3,4], 99) , dim_theta,1 )## RE_lmer_nc[7,4]??

  theta_gam_lmer_CI <- cbind(theta_gam_lmer_est - 1.96*theta_gam_lmer_se, theta_gam_lmer_est + 1.96*theta_gam_lmer_se)
  ########################################################################### 
#
  
  time_lmer2 <- Sys.time()
  time_lmer  <- difftime(time_lmer2,  time_lmer1,  units = 'sec')
  
  
  return(list(
    theta_gamma_REML_coef    = theta_gam_lmer_est,    
    theta_gamma_REML_CI      = theta_gam_lmer_CI,
    omega_REML       = omega_lmer,
    time_REML        = time_lmer
  ) ) 
}