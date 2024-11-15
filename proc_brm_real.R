
proc_brm_real <- function  (y, y_lag1, x, Job_C, Home_C, Dep_B, ID, n_ite, n_burnin )  {
  
  time_brms1 <- Sys.time()

  
  library(cmdstanr)
  library(coda)
  library(brms)
  
  dim_theta <- 4
  dim_gamma <- 9
  
  y_lag1 <- as.matrix(y_lag1)
  Job_C  <- as.matrix(Job_C)
  Home_C <- as.matrix(Home_C)
  Dep_B  <- as.matrix(Dep_B)
  ID     <- as.matrix(ID)
  # suppressWarnings()
  
  ####### brm
  Data <- data.frame( y, y_lag1, x, Job_C, Home_C, Dep_B, ID, N, T ) # using for nlme.
  formula <- bf(
    y ~ (1 + y_lag1 + x) + ( Job_C + Home_C + Dep_B) + ( y_lag1:Job_C + y_lag1:Home_C )  +( x:Job_C + x:Home_C ) + (  1 + y_lag1 + x | ID) ,
    sigma ~ 1 + Job_C + Home_C + (1|ID)
  )
  
  
  fit <- brm( formula = formula,
              data = Data,
              family = gaussian(), 
              chains = 4, 
              cores  = 4, 
              iter = n_ite, 
              warmup = n_burnin
              )
 
  model_summary <- summary(fit)
  model_fixed_effects  <- as.matrix(model_summary$fixed)
  model_random_effects <- as.matrix(model_summary$random$ID)
  Rhat                 <- summary(fit)$fixed[,"Rhat"]
  Rhat                 <- matrix( c( Rhat[1], Rhat[3:4], Rhat[2], Rhat[5:13]  ) ,13,1  )
 
  ## fixed effect
  alpha_brm_est <- model_fixed_effects[1,1]
  alpha_brm_se  <- model_fixed_effects[1,2]
  alpha_brm_CI  <- model_fixed_effects[1,c(3,4)]
  phi_brm_est   <- model_fixed_effects[3,1]
  phi_brm_se    <- model_fixed_effects[3,2]
  phi_brm_CI    <- model_fixed_effects[3,c(3,4)]
  beta_brm_est  <- model_fixed_effects[4,1]
  beta_brm_se   <- model_fixed_effects[4,2]
  beta_brm_CI   <- model_fixed_effects[4,c(3,4)]
  tau_brm_est  <- model_fixed_effects[2,1]
  tau_brm_se   <- model_fixed_effects[2,2]
  tau_brm_CI   <- model_fixed_effects[2,c(3,4)]
  
  alpha_brm_Job_est    <- model_fixed_effects[5,1]
  alpha_brm_Job_se     <- model_fixed_effects[5,2]
  alpha_brm_Job_CI     <- model_fixed_effects[5,c(3,4)]
  alpha_brm_Home_est   <- model_fixed_effects[6,1]
  alpha_brm_Home_se    <- model_fixed_effects[6,2]
  alpha_brm_Home_CI    <- model_fixed_effects[6,c(3,4)]
  alpha_brm_Dep_est    <- model_fixed_effects[7,1]
  alpha_brm_Dep_se     <- model_fixed_effects[7,2]
  alpha_brm_Dep_CI     <- model_fixed_effects[7,c(3,4)]
  
  phi_brm_Job_est     <- model_fixed_effects[8,1]
  phi_brm_Job_se      <- model_fixed_effects[8,2]
  phi_brm_Job_CI      <- model_fixed_effects[8,c(3,4)]
  phi_brm_Home_est    <- model_fixed_effects[9,1]
  phi_brm_Home_se     <- model_fixed_effects[9,2]
  phi_brm_Home_CI     <- model_fixed_effects[9,c(3,4)]
 

  beta_brm_Job_est     <- model_fixed_effects[10,1]
  beta_brm_Job_se      <- model_fixed_effects[10,2]
  beta_brm_Job_CI      <- model_fixed_effects[10,c(3,4)]
  beta_brm_Home_est    <- model_fixed_effects[11,1]
  beta_brm_Home_se     <- model_fixed_effects[11,2]
  beta_brm_Home_CI     <- model_fixed_effects[11,c(3,4)]
  
  tau_brm_Job_est     <- model_fixed_effects[12,1]
  tau_brm_Job_se      <- model_fixed_effects[12,2]
  tau_brm_Job_CI      <- model_fixed_effects[12,c(3,4)]
  tau_brm_Home_est    <- model_fixed_effects[13,1]
  tau_brm_Home_se     <- model_fixed_effects[13,2]
  tau_brm_Home_CI     <- model_fixed_effects[13,c(3,4)]
  
  theta_gamma_CI <- rbind(alpha_brm_CI, phi_brm_CI, beta_brm_CI, tau_brm_CI, alpha_brm_Job_CI, alpha_brm_Home_CI, alpha_brm_Dep_CI, phi_brm_Job_CI, phi_brm_Home_CI, beta_brm_Job_CI, beta_brm_Home_CI, tau_brm_Job_CI, tau_brm_Home_CI  )
  
  ### random effect
  omega_brm <- matrix(c(model_random_effects[1,1]^2, model_random_effects[2,1]^2, model_random_effects[3,1]^2,model_random_effects[4,1]^2),dim_theta,1)
  omega_brm_CI <- model_random_effects[c(1,2,3,4),c(3,4)]
  #####
  
  theta_gam_brm_est <- matrix( c(alpha_brm_est, phi_brm_est, beta_brm_est, tau_brm_est, alpha_brm_Job_est, alpha_brm_Home_est, alpha_brm_Dep_est, phi_brm_Job_est, phi_brm_Home_est, beta_brm_Job_est, beta_brm_Home_est, tau_brm_Job_est, tau_brm_Home_est  ), dim_theta + dim_gamma, 1)
  theta_gam_brm_se  <- matrix( c(alpha_brm_se, phi_brm_se, beta_brm_se, tau_brm_se,  alpha_brm_Job_se, alpha_brm_Home_se, alpha_brm_Dep_se, phi_brm_Job_se, phi_brm_Home_se, beta_brm_Job_se, beta_brm_Home_se, tau_brm_Job_se, tau_brm_Home_se   ),  dim_theta + dim_gamma, 1)
  

  
  time_brms2 <- Sys.time()
  time_brm <- difftime(time_brms2, time_brms1, units = 'sec')
  
  ########################################################################### 
  
  return(list(
    theta_gamma_Bayes_coef   = as.matrix(theta_gam_brm_est),    
    theta_gamma_Bayes_se     = as.matrix(theta_gam_brm_se),
    theta_gamma_Bayes_CI_L   = as.matrix(theta_gamma_CI[,1], (dim_theta+dim_gamma),1 ), 
    theta_gamma_Bayes_CI_U   = as.matrix(theta_gamma_CI[,2], (dim_theta+dim_gamma),1 ),
    theta_Rhat = Rhat,
    omega_Bayes      = omega_brm,
    omega_Bayes_CI_L = as.matrix(omega_brm_CI[,1], dim_theta,1 ),
    omega_Bayes_CI_U = as.matrix(omega_brm_CI[,2], dim_theta,1 ),
    time_Bayes      = time_brm
  ) ) 
}
