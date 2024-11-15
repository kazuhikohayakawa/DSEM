proc_nmle_real <- function (y, y_lag1, x, Job_C, Home_C, Dep_B, Job_c, Home_c, Dep_b, ID, N, T ) {
  
  time_nlme1 <- Sys.time()
  
  library(VAR.etp)
  library(matrixcalc)  
  library(MASS)
  library(sandwich)
  library(e1071)
  
  library(lme4)
  library(nlme)
  
  time_lmer1 <- Sys.time()
  
  dim_theta<- 4
  dim_gamma<- 9
  

  
  ####### nlme_nc
  Data <- data.frame( y, y_lag1, x, Job_C, Home_C, Dep_B, ID, N, T ) # using for nlme.
  reg_nlme <- lme(
    fixed   = y ~ (1 + y_lag1 + x) + ( Job_C + Home_C + Dep_B) + ( y_lag1:Job_C + y_lag1:Home_C )  +( x:Job_C + x:Home_C ) ,
    random  = list(ID = pdSymm(form = ~ 1 + y_lag1 + x)),
    data    = Data,
    weights = varIdent(form = ~ 1 +  Job_C + Home_C |  ID),  # 不组别有不同的方差
    control = list(msMaxIter = 10000, msMaxEval = 10000),
    method  = 'REML' 
  )

  
  
  coef_nlme     <- summary(reg_nlme)$coefficients
  Vcov_nlme     <- as.matrix(vcov(reg_nlme)) # variance of estimated fixed effect
  RE_nlme       <- as.matrix(VarCorr(reg_nlme))  # random effect
  
  weight         <- as.vector(coef(reg_nlme$modelStruct$varStruct, unconstrained=FALSE))
  sig2_1_hat     <- (summary(reg_nlme)$sigma)^2
  sig2_N_hat     <- 1:N
  sig2_N_hat[1]  <- sig2_1_hat
  for (i in 2:N){
    sig2_N_hat[i] <- (summary(reg_nlme)$sigma* weight[i-1])^2
  }
  tau_N_hat <- log(sig2_N_hat)
  # for the estimation of fix and random effect of sig2, and gamma

  
  tau_reg <- lm(tau_N_hat ~ 1 + Job_c + Home_c )
  
  alpha_nlme_est    <- as.vector(coef_nlme$fixed[1])
  alpha_nlme_se     <- sqrt(Vcov_nlme[1,1])
  phi_nlme_est      <- as.vector(coef_nlme$fixed[2])
  phi_nlme_se       <- sqrt(Vcov_nlme[2,2])
  beta_nlme_est     <- as.vector(coef_nlme$fixed[3])
  beta_nlme_se      <- sqrt(Vcov_nlme[3,3])
  tau_nlme_est      <- summary(tau_reg)$coefficients[1,1]
  tau_nlme_se       <- summary(tau_reg)$coefficients[1,2]
  
  
  alpha_nlme_Job_est    <- as.vector(coef_nlme$fixed[4])
  alpha_nlme_Job_se     <- sqrt(Vcov_nlme[4,4])
  alpha_nlme_Home_est   <- as.vector(coef_nlme$fixed[5])
  alpha_nlme_Home_se    <- sqrt(Vcov_nlme[5,5])
  alpha_nlme_Dep_est    <- as.vector(coef_nlme$fixed[6])
  alpha_nlme_Dep_se     <- sqrt(Vcov_nlme[6,6])
  
  phi_nlme_Job_est    <- as.vector(coef_nlme$fixed[7])
  phi_nlme_Job_se     <- sqrt(Vcov_nlme[7,7])
  phi_nlme_Home_est   <- as.vector(coef_nlme$fixed[8])
  phi_nlme_Home_se    <- sqrt(Vcov_nlme[8,8])
  
  beta_nlme_Job_est    <- as.vector(coef_nlme$fixed[9])
  beta_nlme_Job_se     <- sqrt(Vcov_nlme[9,9])
  beta_nlme_Home_est   <- as.vector(coef_nlme$fixed[10])
  beta_nlme_Home_se    <- sqrt(Vcov_nlme[10,10])
  
  tau_nlme_Job_est   <- as.vector(tau_reg$coefficients[2]) 
  tau_nlme_Job_se    <- summary(tau_reg)$coefficients[2,2]
  tau_nlme_Home_est  <- as.vector(tau_reg$coefficients[3]) 
  tau_nlme_Home_se   <- summary(tau_reg)$coefficients[3,2]
  
  theta_gam_nlme_est  <- matrix( c(alpha_nlme_est, phi_nlme_est, beta_nlme_est, tau_nlme_est,alpha_nlme_Job_est,alpha_nlme_Home_est, alpha_nlme_Dep_est,phi_nlme_Job_est, phi_nlme_Home_est, beta_nlme_Job_est, beta_nlme_Home_est, tau_nlme_Job_est, tau_nlme_Home_est  ),dim_theta+dim_gamma,1)
  theta_gam_nlme_se   <- matrix( c(alpha_nlme_se, phi_nlme_se, beta_nlme_se, tau_nlme_se, alpha_nlme_Job_se, alpha_nlme_Home_se, alpha_nlme_Dep_se, phi_nlme_Job_se, phi_nlme_Home_se, beta_nlme_Job_se, beta_nlme_Home_se, tau_nlme_Job_se, tau_nlme_Home_se),dim_theta+dim_gamma,1)
  omega_nlme          <- matrix( as.numeric(c(RE_nlme[1,1], RE_nlme[2,1], RE_nlme[3,1], RE_nlme[4,1] ), dim_theta, 1))
  
  theta_gam_nlme_CI <- cbind(theta_gam_nlme_est - 1.96*theta_gam_nlme_se, theta_gam_nlme_est + 1.96*theta_gam_nlme_se)
  ########################################################################### 
  time_nlme2 <- Sys.time()
  time_nlme  <- difftime(time_nlme2,  time_nlme1,  units = 'sec')
  
  return(list(
    theta_gamma_REML_coef    = theta_gam_nlme_est,    
    theta_gamma_REML_CI      = theta_gam_nlme_CI,
    omega            = omega_nlme,
    time_REML        = time_nlme
  ) ) 
}