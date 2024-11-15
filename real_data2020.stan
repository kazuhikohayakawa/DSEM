 data {
 int<lower=0> NT;
 int<lower=0> T;
 int<lower=0> N;
 vector[NT] y;   
 vector[NT] x;   
 vector[NT] y_lag;   
 vector[N] Job_c;
 vector[N] Home_c;
 vector[N] Dep_b;
 }
 
 parameters {
   real alpha;
   real<lower = -1, upper = 1> phi;
   real beta;
   // real<lower=0> tau2;
   real tau;
   vector[9] gam;

   vector[N] alpha_N;
   vector<lower = -1, upper = 1>[N] phi_N;
   vector[N] beta_N;
   vector[N] tau_N;
   // vector<lower=0>[N]  tau2_N;

   real<lower=0> ome2_alpha;
   real<lower=0> ome2_phi;
   real<lower=0> ome2_beta;
   real<lower=0> ome2_tau;
 }
 

 model {
     for (n in 1:N){
            alpha_N[n] ~ normal( alpha + gam[1]* Job_c[n] + gam[2]* Home_c[n] + gam[3]* Dep_b[n],  sqrt(ome2_alpha) );
            phi_N[n]   ~ normal( phi   + gam[4]* Job_c[n] + gam[5]* Home_c[n] , sqrt(ome2_phi) );
            beta_N[n]  ~ normal( beta  + gam[6]* Job_c[n] + gam[7]* Home_c[n] , sqrt(ome2_beta) );
            tau_N[n]   ~ normal( tau   + gam[8]* Job_c[n] + gam[9]* Home_c[n] , sqrt(ome2_tau) );
}
     for (n in 1:N){
      for (t in 1:T){
        y[  t+ T*(n-1)  ]~ normal(   alpha_N[n]  + y_lag[  t+ T*(n-1)  ] * phi_N[n]  +  beta_N[n] *x[  t+ T*(n-1)  ]  ,  sqrt( exp(tau_N[n]))  );
         }
        }
}
