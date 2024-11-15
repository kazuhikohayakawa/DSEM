proc_stan_real <- function  (y, y_lag1, x, Job_c, Home_c, Dep_b, ID, N, T, n_ite, n_burnin )  {
  
  time_stan1 <- Sys.time()

  
  library(cmdstanr)
  library(coda)
  
  dim_theta <- 4
  dim_gamma <- 9


  # suppressWarnings()
  dat <- list(
    NT = N*(T),
    T  =  T, 
    N  = N,   
    y  = as.vector(y),
    x  = as.vector(x),
    y_lag  =  as.vector(y_lag1),
    Job_c  = as.vector(Job_c),
    Home_c = as.vector(Home_c),
    Dep_b  = as.vector(Dep_b)
  )
  
  # n_ite    <- 500
  # n_burnin <- 100
  
  model  <- cmdstan_model('real_data2020.stan')
  fit_CmdStanR1 <- model$sample(data=dat,
                                # init           = init2,
                                # chains          = n_chains,
                                parallel_chains = 4,
                                iter_warmup     = n_burnin,
                                iter_sampling   = n_ite - n_burnin,
                                seed=1234,
                                refresh = 0)
  
  params_for_summary <- c("alpha",  "phi",   "beta", "tau", "gam", "ome2_alpha", "ome2_phi", "ome2_beta", "ome2_tau" )
  summary_CmdStanR1 <- fit_CmdStanR1$summary(params_for_summary, "mean", "sd", function(x){  quantile(x, probs = c(0.025, 0.975), names = TRUE)}, "rhat")
  
  ######### result
  result_stan <- as.matrix(summary_CmdStanR1[,2:6])
  rownames(result_stan) <- c("alpha",  "phi",   "beta", "tau", "gam","gam","gam","gam","gam","gam","gam","gam","gam", "ome2_alpha", "ome2_phi", "ome2_beta", "ome2_tau"  )
  
 
  
  time_stan2 <- Sys.time()
  time_stan <- difftime(time_stan2, time_stan1, units = 'sec')
  
  ########################################################################### 
  
  return(list(
    theta_gamma_Bayes_fixed     = result_stan[1:(dim_theta+dim_gamma),],   
    theta_gamma_Bayes_random    = result_stan[(dim_theta+dim_gamma+1):(dim_theta+dim_gamma+4) ,],                  
    time_Bayes      = time_stan
  ) ) 
}