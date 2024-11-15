# proc_ARX1_OLS_KP1_KP2
proc_reg_OLS_JK_BC_real <- function(Y_TN, X_TN, N) {
  
  library(VAR.etp)
  library(matrixcalc)  
  library(MASS)
  library(sandwich)
  library(e1071)
  
  dim <- dim(Y_TN)
  TT<- dim[1]
  T <- TT -1
  dim_theta <-4
  
  theta_OLS_N  <- matrix(0, N, dim_theta)
  theta_KP1_N  <- matrix(0, N, dim_theta)
  theta_KP2_N  <- matrix(0, N, dim_theta)
  theta_JK_N   <- matrix(0, N, dim_theta)
  
  var_OLS_N   <- matrix(0, N, dim_theta)
  var_JK_N    <- matrix(0, N, dim_theta)
  var_KP1_N   <- matrix(0, N, dim_theta)
  var_KP2_N   <- matrix(0, N, dim_theta)
  
  
  for (i in 1:N){
    Y_i  <- as.matrix(Y_TN[2:TT,i],T,1)
    Y1_i <- as.matrix(Y_TN[1:T,i], T,1)
    X_i  <- as.matrix(X_TN[2:TT,i],T,1)
    
    # # within-centering
    Y1_i <- Y1_i - mean(Y1_i)
    X_i  <- X_i  - mean(X_i)
    
    W_i  <- cbind(matrix(1,T,1), Y1_i, X_i)
    
    ### OLS ###
    est_OLS_i <- lm(Y_i ~ 1 + Y1_i + X_i )
    
    alpha_OLS_i <- est_OLS_i$coef[1] # intercept
    phi_OLS_i   <- est_OLS_i$coef[2] # AR(1) coefficient
    beta_OLS_i  <- est_OLS_i$coef[3] # coef. of x_it
    u_OLS_i     <- as.vector(est_OLS_i $ residuals) 
    sig2_OLS_i  <- sum(u_OLS_i^2)/ (est_OLS_i$ df.residual)
    tau2_OLS_i  <- log(sig2_OLS_i)
    theta_OLS_i <- c(alpha_OLS_i, phi_OLS_i, beta_OLS_i, tau2_OLS_i)
    
    delta_OLS_i <- as.vector(rbind(phi_OLS_i, beta_OLS_i, alpha_OLS_i))
    
    
    theta_OLS_N[i,] <- theta_OLS_i
    
    var_OLS_i      <- sig2_OLS_i*solve(t(W_i) %*% W_i)
    var_OLS_tau2_i <- ( (2 + kurtosis(u_OLS_i, type=2)))/T 
    var_OLS_N[i,] <- c(diag(var_OLS_i), var_OLS_tau2_i)
    
    
    ############################
    ### KP1 ###
    ####  bias-correction1  ####
    ones <- matrix(rep(1,T),T,1)
    W    <-  cbind(Y1_i, X_i, ones )  
    F    <- matrix(0,T,1)
    C    <- matrix(0,T,T)
    e1   <- matrix(c(1,0,0),3,1)
    for (t in 1:T){
      F[t] <- phi_OLS_i^(t-1)
    }
    for (t in 1:(T-1) ) {
      C[,t]  <- rbind(matrix(rep(0,t),t,1), matrix(F[1: (T-t) ],T-t,1)) 
    }
    Z <- cbind(Y1_i[1] *F + C %*% cbind( X_i, ones ) %*% matrix(c(beta_OLS_i, alpha_OLS_i),2,1), cbind( X_i, ones ) )
    D <- t(Z)%*%Z + sig2_OLS_i* sum( diag( t(C) %*% C)) * (e1 %*% t(e1))  
    invD <- ginv(D)
    bias <- - sig2_OLS_i *invD %*%( t(Z)%*% C %*% Z %*% invD %*% e1 +  sum(diag(t(Z) %*% C %*% Z %*% invD ))*e1 
                                    + as.vector(2*sig2_OLS_i* t(e1) %*% invD %*% e1 * sum(diag( (C %*% t(C) %*% C)))) * e1 )  
    #######
    delta_KP1_i <- delta_OLS_i - bias
    u_KP1_i     <- Y_i - W %*% delta_KP1_i
    sig2_KP1_i  <- as.vector(t(u_KP1_i) %*% u_KP1_i/ (T- 3)) 
    tau2_KP1_i  <- log(sig2_KP1_i)
    
    theta_KP1_i <- c(delta_KP1_i[3], delta_KP1_i[1], delta_KP1_i[2], tau2_KP1_i)
    theta_KP1_N[i,]  <- theta_KP1_i
    
    var_KP1_i      <- sig2_KP1_i*solve(t(W_i) %*% W_i)
    var_KP1_tau2_i <- ( (2 + kurtosis(u_KP1_i, type=2)))/T 
    var_KP1_N[i,]  <- c(diag(var_KP1_i), var_KP1_tau2_i)
    
    
    ##################
    # KP2
    Lam <- diag(rep(1,T))
    
    for (t in 1:(T-1) ){
      Lam[t+1,t] = -phi_OLS_i
    }
    ome = 0
    Omega = diag(rep(T+1))
    Omega[1,1] = ome
    G = cbind( matrix(rep(0,T),T,1), C )
    
    Q = solve(D)
    q1 = Q %*% e1
    q11 = as.vector(t(e1) %*% Q %*% e1)
    b2a = sum(diag(Q %*% t(Z) %*% C %*% Z  )) * q1
    b2b = Q %*% t(Z) %*% C %*% Z %*% q1
    
    b4a = -2*q11*sum (diag(G %*% t(G) %*% C))
    b4b = 2*q11*sum(diag(Q %*% t(Z) %*% G %*% t(G) %*% C %*% Z))
    b4c = 2*q11*sum( diag(Q %*% t(Z) %*% G %*% t(G) %*% t(C) %*% Z))
    b4d = -2*q11*sum(diag(Q %*% t(Z) %*% G %*% t(G) %*% Z %*% Q %*% t(Z) %*% C %*% Z))
    b4e = -q11*sum(diag(Q %*% t(Z) %*% C%*% Z)) %*% sum(diag((Q %*% t(Z) %*%G %*% t(G)%*%Z)))
    b4f =  4 * t(q1) %*% t(Z) %*% G %*% t(G) %*% C %*% Z %*% q1
    b4g = 2*t(q1) %*% t(Z) %*% G %*% t(G) %*% t(C) %*% Z %*% q1
    b4h = - 4*t(q1) %*% t(Z) %*% G %*% t(G) %*% Z %*% Q %*% t(Z) %*% C %*% Z %*% q1
    b4i = -2*t(q1) %*% t(Z) %*% G %*% t(G) %*% Z %*% Q %*% t(Z) %*% t(C) %*% Z %*% q1
    
    b4j = -t(q1) %*% t(Z) %*% G %*% t(G) %*% Z %*% q1 %*% sum(diag(Q %*% t(Z) %*%C %*% Z))
    b4k = - 2*t(q1) %*% t(Z) %*% C %*% Z %*% q1 %*% sum(diag(Q %*% t(Z) %*% G %*% t(G) %*% Z))
    b4l = -as.vector( q11*sum(diag(Q %*% t(Z) %*% G %*% t(G) %*% Z)) + t(q1) %*% t(Z) %*%  G %*% t(G) %*% Z %*% q1 ) * Q %*% t(Z) %*% C %*% Z %*% q1
    b4m = -2*as.vector( (q11*sum(diag((Q %*% t(Z) %*% C %*% Z))) + t(q1) %*% t(Z) %*% C %*% Z %*% q1)) * Q %*% t(Z) %*% G %*% t(G) %*% Z %*% q1
    b4n = 2 * q11 * Q %*% t(Z) %*% (G %*% t(G) %*% C + C %*% G %*% t(G) + G %*% t(G) %*% t(C)) %*% Z %*% q1
    b4o = -2*q11 * Q %*% t(Z) %*% G %*% t(G) %*% Z %*% Q %*% t(Z) %*% (C+t(C)) %*% Z %*% q1
    b4p = -2*q11 * Q %*% t(Z) %*% C %*% Z %*% Q %*% t(Z) %*% G %*% t(G) %*% Z %*% q1
    
    b6a =  8*(q11^2) *sum(diag(G %*% t(G) %*% G %*% t(G) %*%C))
    b6b = -2*q11^2 *sum(diag(G %*% t(G) %*% G %*% t(G))) *sum(diag(Q %*% t(Z) %*% C %*% Z))
    b6c = -4 * q11^2 * sum(diag(G %*% t(G) %*% C)) %*% sum(diag(Q %*% t(Z) %*% G %*% t(G) %*% Z))
    b6d = -12 *q11* as.vector(t(q1) %*% t(Z) %*% G %*% t(G) %*% Z %*% q1) *sum(diag(G %*% t(G)%*% C))
    b6e = -8*q11 * (t(q1) %*% t(Z) %*% C %*% Z %*% q1) * sum(diag(G %*% t(G) %*% G %*% t(G) ))
    
    b6f = -2* q11^2 *sum(diag(G %*% t(G) %*% G %*% t(G))) * Q %*% t(Z) %*% C %*% Z %*% q1
    b6g = -8* q11^2 *sum(diag(G %*% t(G) %*% C )) * Q %*% t(Z) %*% G %*% t(G) %*% Z %*% q1
    
    b8a = 12*q11^3* sum(diag(G %*% t(G) %*%C)) *sum(diag( G %*% t(G) %*% G %*% t(G) )) * q1
    
    s2 = sig2_OLS_i
    
    bias = -s2*( b2a + b2b )
    +s2^2*(  as.vector(b4a + b4b + b4c + b4d + b4e + b4f + b4g + b4h + b4i + b4j + b4k) * q1 + b4l + b4m + b4n + b4o + b4p)
    +s2^3* (as.vector( (b6a + b6b + b6c + b6d + b6e))*q1 + b6f+b6g ) 
    - s2^4*(b8a)
    
    #######
    delta_KP2_i  <- delta_OLS_i - as.vector(bias)
    
    u_KP2_i  <- Y_i - W %*% delta_KP2_i
    sig2_KP2_i <- as.vector(t(u_KP2_i) %*% u_KP2_i/(T-3))
    tau2_KP2_i <- log(sig2_KP2_i)
    
    theta_KP2_i<- c(delta_KP2_i[3], delta_KP2_i[1], delta_KP2_i[2], tau2_KP2_i )
    theta_KP2_N[i,]  <- theta_KP2_i
    
    var_KP2_i      <- sig2_KP2_i*solve(t(W_i) %*% W_i)
    var_KP2_tau2_i <- ( (2 + kurtosis(u_KP2_i, type=2)))/T 
    var_KP2_N[i,] <- c(diag(var_KP2_i), var_KP2_tau2_i)
    
    
    ######### JKMG ########
    T_half <- floor(T/2)
    Y_i_1  <- as.vector( Y_i[1:T_half] )
    Y1_i_1 <- as.vector( Y1_i[1:T_half] )
    X_i_1  <- as.vector( X_i[1:T_half] )
    
    Y_i_2  <- as.vector( Y_i[(T_half+1):T] )
    Y1_i_2 <- as.vector( Y1_i[(T_half+1):T] )
    X_i_2  <- as.vector( X_i[(T_half+1):T]  )
    
    lm1 <- lm(Y_i_1 ~ 1+ Y1_i_1 + X_i_1 )
    alpha_OLS1_i <- lm1$coef[1]
    phi_OLS1_i   <- lm1$coef[2]
    beta_OLS1_i  <- lm1$coef[3]
    u_OLS1_i     <- as.vector(lm1 $ residuals) 
    sig2_OLS1_i  <- sum(u_OLS1_i^2)/ (lm1$ df.residual)
    tau2_OLS1_i  <- log(sig2_OLS1_i)
    theta_OLS1_i <- c(alpha_OLS1_i, phi_OLS1_i, beta_OLS1_i, tau2_OLS1_i)
    
    lm2 <- lm(Y_i_2 ~ 1+ Y1_i_2 + X_i_2 )
    alpha_OLS2_i <- lm2$coef[1]
    phi_OLS2_i   <- lm2$coef[2]
    beta_OLS2_i  <- lm2$coef[3]
    u_OLS2_i     <- as.vector(lm2 $ residuals) 
    sig2_OLS2_i  <- sum(u_OLS2_i^2)/ (lm2$ df.residual)
    tau2_OLS2_i  <- log(sig2_OLS2_i)
    theta_OLS2_i <- c(alpha_OLS2_i, phi_OLS2_i, beta_OLS2_i, tau2_OLS2_i)
    
    theta_JK_i     <- 2*theta_OLS_i - (theta_OLS1_i + theta_OLS2_i)/2
    theta_JK_N[i,] <- theta_JK_i
    sig2_JK_i      <- exp(theta_JK_i[4])
    
    # avg_W_a <- matrix( rep(apply(W_i,2,mean), each=T_half), T_half, 3)
    # avg_W_b <- matrix( rep(apply(W_i,2,mean), each=T-T_half), T-T_half, 3)
    # avg_W_1 <- matrix( rep(apply(W_i[1:T_half, 1:3],2,mean), each=T_half), T_half, 3)
    # avg_W_2 <- matrix( rep(apply(W_i[(T_half+1):T, 1:3],2,mean), each=T-T_half), T-T_half, 3)
    # 
    # W_JK_i_a <- W_i[1:T_half,1:3]     - (2*avg_W_a - avg_W_1)
    # W_JK_i_b <- W_i[(T_half+1):T,1:3] - (2*avg_W_b - avg_W_2)
    # 
    # W_JK_i <- as.matrix(rbind(W_JK_i_a, W_JK_i_b), T,3)
    # W_JK_i <- cbind(matrix(1,T,1), W_JK_i[,2:3])
    
    # var_JK_i  <- sig2_JK_i*solve(t(W_i)%*%W_i)%*%(t(W_JK_i)%*%(W_JK_i))%*%solve(t(W_i)%*%W_i);
    var_JK_i  <- sig2_JK_i*solve(t(W_i)%*%W_i)
    
    u_JK_i    <- Y_i - W_i%*%theta_JK_i[1:3]
    var_JK_tau2_i <- ( (2 + kurtosis(u_JK_i, type=2)))/T 
    var_JK_N[i,] <- c(diag(var_JK_i), var_JK_tau2_i)
    
    
  }
  
  ########################################################################### 
  
  
  return(list(theta_OLS_N, theta_JK_N,  theta_KP1_N,  theta_KP2_N, var_OLS_N, var_JK_N,  var_KP1_N,  var_KP2_N  ))
}