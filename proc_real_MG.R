# proc_ARX1_OLS_KP1_KP2
proc_MG_real <- function (Y_TN, X_TN, Job_c, Home_c, Dep_b, ID, N ) {
  
  library(VAR.etp)
  library(matrixcalc)  
  library(MASS)
  library(sandwich)
  library(e1071)
  
  source("proc_reg_OLS_JK_BC_real.R")
  
  time_MG1 <- Sys.time()
  
  dim <- dim(Y_TN)
  TT<- dim[1]
  T <- TT -1

  # dim_theta1 <- 3 + nrow(X_TN)
  dim_Z <- 2
  
  dim_theta <- 4
  dim_gamma <- 9
  dim_reg <- 3  # number of regressors (cons + y(-1) + X)

  
  ########################################################################### 
  centering <- 0

  theta_hat_N <- proc_reg_OLS_JK_BC_real(Y_TN, X_TN, N )
  
  theta_OLS_N <- theta_hat_N[[1]] # N x dim(theta)
  theta_JK_N  <- theta_hat_N[[2]] # N x dim(theta)
  theta_KP1_N <- theta_hat_N[[3]] # N x dim(theta)
  theta_KP2_N <- theta_hat_N[[4]] # N x dim(theta)
  var_OLS_N   <- theta_hat_N[[5]] # N x dim(theta)
  var_JK_N    <- theta_hat_N[[6]] # N x dim(theta)
  var_KP1_N   <- theta_hat_N[[7]] # N x dim(theta)
  var_KP2_N   <- theta_hat_N[[8]] # N x dim(theta)
  
  
  
########################################################################### 
  theta_gamma_MG_coef_0    <- matrix(0, dim_Z+2, 1)
  theta_gamma_JKMG_coef_0  <- matrix(0, dim_Z+2, 1)
  theta_gamma_BCMG1_coef_0 <- matrix(0, dim_Z+2, 1)
  theta_gamma_BCMG2_coef_0 <- matrix(0, dim_Z+2, 1)
  
  theta_gamma_MG_se_0    <- matrix(0, dim_Z+2, 1)
  theta_gamma_JKMG_se_0  <- matrix(0, dim_Z+2, 1)
  theta_gamma_BCMG1_se_0 <- matrix(0, dim_Z+2, 1)
  theta_gamma_BCMG2_se_0 <- matrix(0, dim_Z+2, 1)
  
  theta_gamma_MG_coef    <- matrix(0, dim_Z+1, dim_theta-1)
  theta_gamma_JKMG_coef  <- matrix(0, dim_Z+1, dim_theta-1)
  theta_gamma_BCMG1_coef <- matrix(0, dim_Z+1, dim_theta-1)
  theta_gamma_BCMG2_coef <- matrix(0, dim_Z+1, dim_theta-1)

  theta_gamma_MG_se    <- matrix(0, dim_Z+1, dim_theta-1)
  theta_gamma_JKMG_se  <- matrix(0, dim_Z+1, dim_theta-1)
  theta_gamma_BCMG1_se <- matrix(0, dim_Z+1, dim_theta-1)
  theta_gamma_BCMG2_se <- matrix(0, dim_Z+1, dim_theta-1)


  resid_MG     <- matrix(0, N, dim_theta) 
  resid_JKMG   <- matrix(0, N, dim_theta) 
  resid_BCMG1  <- matrix(0, N, dim_theta) 
  resid_BCMG2  <- matrix(0, N, dim_theta)
  
  # k=1
  
  reg_MG  <- lm(theta_OLS_N[,1] ~ 1 + Job_c + Home_c + Dep_b)
  theta_gamma_MG_coef_0 <- summary(reg_MG)$coef[,1]
  theta_gamma_MG_se_0   <- summary(reg_MG)$coef[,2]
  resid_MG[,1]          <- reg_MG$resid
  
  reg_JKMG  <- lm(theta_JK_N[,1] ~ 1 + Job_c + Home_c + Dep_b)
  theta_gamma_JKMG_coef_0 <- summary(reg_JKMG)$coef[,1]
  theta_gamma_JKMG_se_0   <- summary(reg_JKMG)$coef[,2]
  resid_JKMG[,1]          <- reg_JKMG$resid
  
  reg_BCMG1  <- lm(theta_KP1_N[,1] ~ 1 + Job_c + Home_c + Dep_b)
  theta_gamma_BCMG1_coef_0 <- summary(reg_BCMG1)$coef[,1]
  theta_gamma_BCMG1_se_0   <- summary(reg_BCMG1)$coef[,2]
  resid_BCMG1[,1]           <- reg_BCMG1$resid
  
  reg_BCMG2  <- lm(theta_KP2_N[,1] ~ 1 + Job_c + Home_c + Dep_b)
  theta_gamma_BCMG2_coef_0 <- summary(reg_BCMG2)$coef[,1]
  theta_gamma_BCMG2_se_0   <- summary(reg_BCMG2)$coef[,2]
  resid_BCMG2[,1]            <- reg_BCMG2$resid
  
  
  for ( k in (2:dim_theta) ){
      reg_MG  <- lm(theta_OLS_N[,k] ~ 1 + Job_c + Home_c)
      theta_gamma_MG_coef[,k-1] <- summary(reg_MG)$coef[,1]
      theta_gamma_MG_se[,k-1]   <- summary(reg_MG)$coef[,2]
      resid_MG[,k]            <- reg_MG$resid

      reg_JKMG  <- lm(theta_JK_N[,k] ~ 1 + Job_c + Home_c)
      theta_gamma_JKMG_coef[,k-1] <- summary(reg_JKMG)$coef[,1]
      theta_gamma_JKMG_se[,k-1]   <- summary(reg_JKMG)$coef[,2]
      resid_JKMG[,k]            <- reg_JKMG$resid
  
      reg_BCMG1  <- lm(theta_KP1_N[,k] ~ 1 + Job_c + Home_c)
      theta_gamma_BCMG1_coef[,k-1] <- summary(reg_BCMG1)$coef[,1]
      theta_gamma_BCMG1_se[,k-1]   <- summary(reg_BCMG1)$coef[,2]
      resid_BCMG1[,k]            <- reg_BCMG1$resid
  
      reg_BCMG2  <- lm(theta_KP2_N[,k] ~ 1 + Job_c + Home_c)
      theta_gamma_BCMG2_coef[,k-1] <- summary(reg_BCMG2)$coef[,1]
      theta_gamma_BCMG2_se[,k-1]   <- summary(reg_BCMG2)$coef[,2]
      resid_BCMG2[,k]            <- reg_BCMG2$resid
  }


  
  theta_MG_coef    <- matrix( c( theta_gamma_MG_coef_0[1], theta_gamma_MG_coef[1,]),    dim_theta, 1)
  theta_JKMG_coef  <- matrix( c( theta_gamma_JKMG_coef_0[1], theta_gamma_JKMG_coef[1,]),  dim_theta, 1)
  theta_BCMG1_coef <- matrix( c( theta_gamma_BCMG1_coef_0[1],theta_gamma_BCMG1_coef[1,]), dim_theta, 1)
  theta_BCMG2_coef <- matrix( c( theta_gamma_BCMG2_coef_0[1],theta_gamma_BCMG2_coef[1,]), dim_theta, 1)
  
  gamma_MG_coef    <- matrix( c( theta_gamma_MG_coef_0[2:(dim_Z+2)], theta_gamma_MG_coef[2:(dim_Z+1),]),  dim_gamma  , 1 )
  gamma_JKMG_coef  <- matrix( c( theta_gamma_JKMG_coef_0[2:(dim_Z+2)], theta_gamma_JKMG_coef[2:(dim_Z+1),]),  dim_gamma, 1 )
  gamma_BCMG1_coef <- matrix( c( theta_gamma_BCMG1_coef_0[2:(dim_Z+2)], theta_gamma_BCMG1_coef[2:(dim_Z+1),]), dim_gamma, 1 )
  gamma_BCMG2_coef <- matrix( c( theta_gamma_BCMG2_coef_0[2:(dim_Z+2)], theta_gamma_BCMG2_coef[2:(dim_Z+1),]), dim_gamma, 1 )
  
  theta_gamma_MG_coef     <- rbind(theta_MG_coef,    gamma_MG_coef)
  theta_gamma_JKMG_coef   <- rbind(theta_JKMG_coef,  gamma_JKMG_coef)
  theta_gamma_BCMG1_coef  <- rbind(theta_BCMG1_coef, gamma_BCMG1_coef)
  theta_gamma_BCMG2_coef  <- rbind(theta_BCMG2_coef, gamma_BCMG2_coef)
  
  theta_MG_se    <- matrix( c( theta_gamma_MG_se_0[1], theta_gamma_MG_se[1,]),    dim_theta, 1)
  theta_JKMG_se  <- matrix( c( theta_gamma_JKMG_se_0[1], theta_gamma_JKMG_se[1,]),  dim_theta, 1)
  theta_BCMG1_se <- matrix( c( theta_gamma_BCMG1_se_0[1],theta_gamma_BCMG1_se[1,]), dim_theta, 1)
  theta_BCMG2_se <- matrix( c( theta_gamma_BCMG2_se_0[1],theta_gamma_BCMG2_se[1,]), dim_theta, 1)
  
  
  gamma_MG_se    <- matrix( c( theta_gamma_MG_se_0[2:(dim_Z+2)], theta_gamma_MG_se[2:(dim_Z+1),]),  dim_gamma  , 1 )
  gamma_JKMG_se  <- matrix( c( theta_gamma_JKMG_se_0[2:(dim_Z+2)], theta_gamma_JKMG_se[2:(dim_Z+1),]),  dim_gamma, 1 )
  gamma_BCMG1_se <- matrix( c( theta_gamma_BCMG1_se_0[2:(dim_Z+2)], theta_gamma_BCMG1_se[2:(dim_Z+1),]), dim_gamma, 1 )
  gamma_BCMG2_se <- matrix( c( theta_gamma_BCMG2_se_0[2:(dim_Z+2)], theta_gamma_BCMG2_se[2:(dim_Z+1),]), dim_gamma, 1 )
  
  theta_gamma_MG_se     <- rbind(theta_MG_se,    gamma_MG_se)
  theta_gamma_JKMG_se   <- rbind(theta_JKMG_se,  gamma_JKMG_se)
  theta_gamma_BCMG1_se  <- rbind(theta_BCMG1_se, gamma_BCMG1_se)
  theta_gamma_BCMG2_se  <- rbind(theta_BCMG2_se, gamma_BCMG2_se)
  

  results_coef_se <- rbind(theta_gamma_MG_coef,    theta_gamma_MG_se,    theta_gamma_JKMG_coef,  theta_gamma_JKMG_se,
                           theta_gamma_BCMG1_coef, theta_gamma_BCMG1_se, theta_gamma_BCMG2_coef, theta_gamma_BCMG2_se )
  
  # print(results_coef_se)
  
###############################################################################
# Estimation of random effects

  
  Omega_MG    <- cov(resid_MG)
  Omega_JKMG  <- cov(resid_JKMG)
  Omega_BCMG1 <- cov(resid_BCMG1)
  Omega_BCMG2 <- cov(resid_BCMG2)

  # Omega_bc_MG     <- Omega_MG    
  # Omega_bc_JKMG   <- Omega_JKMG  
  # Omega_bc_BCMG1  <- Omega_BCMG1 
  # Omega_bc_BCMG2  <- Omega_BCMG2 
  
  
  Omega_bc_MG     <- Omega_MG    - diag(c(apply(var_OLS_N[,], 2, mean) ))
  Omega_bc_JKMG   <- Omega_JKMG  - diag(c(apply(var_JK_N[,],  2, mean) ))
  Omega_bc_BCMG1  <- Omega_BCMG1 - diag(c(apply(var_KP1_N[,], 2, mean) ))
  Omega_bc_BCMG2  <- Omega_BCMG2 - diag(c(apply(var_KP2_N[,], 2, mean) ))
   

  omega_MG    <- matrix(diag(Omega_MG),    dim_theta, 1)
  omega_JKMG  <- matrix(diag(Omega_JKMG),  dim_theta, 1)
  omega_BCMG1 <- matrix(diag(Omega_BCMG1), dim_theta, 1)
  omega_BCMG2 <- matrix(diag(Omega_BCMG2), dim_theta, 1)
  
  omega_bc_MG    <- matrix(diag(Omega_bc_MG),    dim_theta, 1)
  omega_bc_JKMG  <- matrix(diag(Omega_bc_JKMG),  dim_theta, 1)
  omega_bc_BCMG1 <- matrix(diag(Omega_bc_BCMG1), dim_theta, 1)
  omega_bc_BCMG2 <- matrix(diag(Omega_bc_BCMG2), dim_theta, 1)
  

  # print(Omega_MG)
  
  # n_reg2 <- dim_theta*(dim_theta+1)/2
  # omega_MG    <- matrix(0, n_reg2, 1) 
  # omega_JKMG  <- matrix(0, n_reg2, 1) 
  # omega_BCMG1 <- matrix(0, n_reg2, 1) 
  # omega_BCMG2 <- matrix(0, n_reg2, 1) 
  # 
  # for (j in (0:(dim_theta-2)) ) {
  #   omega_MG[((j*dim_theta)-j*(j-1)/2+1)    : ((j+1)*dim_theta-j*(j+1)/2)] <- diag(Omega_MG[(j+1)   :dim_theta, 1:(dim_theta-j)])
  #   omega_JKMG[((j*dim_theta)-j*(j-1)/2+1)  : ((j+1)*dim_theta-j*(j+1)/2)] <- diag(Omega_JKMG[(j+1) :dim_theta, 1:(dim_theta-j)])
  #   omega_BCMG1[((j*dim_theta)-j*(j-1)/2+1) : ((j+1)*dim_theta-j*(j+1)/2)] <- diag(Omega_BCMG1[(j+1):dim_theta, 1:(dim_theta-j)])
  #   omega_BCMG2[((j*dim_theta)-j*(j-1)/2+1) : ((j+1)*dim_theta-j*(j+1)/2)] <- diag(Omega_BCMG2[(j+1):dim_theta, 1:(dim_theta-j)])
  # }
  # 
  # omega_MG[n_reg2]    <- Omega_MG[dim_theta, 1]
  # omega_JKMG[n_reg2]  <- Omega_JKMG[dim_theta, 1]
  # omega_BCMG1[n_reg2] <- Omega_BCMG1[dim_theta, 1]
  # omega_BCMG2[n_reg2] <- Omega_BCMG2[dim_theta, 1]

  
  results_omega     <- rbind(omega_MG,     omega_JKMG,     omega_BCMG1,     omega_BCMG2)
  results_omega_bc  <- rbind(omega_bc_MG, omega_bc_JKMG,   omega_bc_BCMG1,  omega_bc_BCMG2)

  time_MG2 <- Sys.time()
  time_MG  <- difftime(time_MG2,  time_MG1,  units = 'sec')
  

  # return(list(theta_gamma_MG_coef,    theta_gamma_MG_se,  theta_gamma_JKMG_coef,  theta_gamma_JKMG_se))
  return(list(
    theta_gamma_MG_coef    = theta_gamma_MG_coef,    
    theta_gamma_MG_se      = theta_gamma_MG_se,
    theta_gamma_JKMG_coef  = theta_gamma_JKMG_coef,  
    theta_gamma_JKMG_se    = theta_gamma_JKMG_se,
    theta_gamma_BCMG1_coef = theta_gamma_BCMG1_coef, 
    theta_gamma_BCMG1_se   = theta_gamma_BCMG1_se, 
    theta_gamma_BCMG2_coef = theta_gamma_BCMG2_coef, 
    theta_gamma_BCMG2_se   = theta_gamma_BCMG2_se, 
    omega_MG       = omega_MG,
    omega_JKMG     = omega_JKMG,
    omega_BCMG1    = omega_BCMG1,
    omega_BCMG2    = omega_BCMG2,
    omega_bc_MG    = omega_bc_MG,
    omega_bc_JKMG  = omega_bc_JKMG,
    omega_bc_BCMG1 = omega_bc_BCMG1,  
    omega_bc_BCMG2 = omega_bc_BCMG2,
    time_MG        = time_MG
    ) ) 
}