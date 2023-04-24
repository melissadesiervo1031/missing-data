// 
//
/*----------------------- Data --------------------------*/
  data {
      int<lower = 1> N;   // Length of state and observation time series
      vector[N] P_obs;    // observed GPP
      vector[N] light;    // covariate matrix
      vector[N] Q;        // covariate matrix
      vector[N] sdo;      // standard deviation of observation error
      vector[N] miss_vec; // vector of 0's (missing) and 1's (data) describing
                          // locaions of missing GPP data

  }
 
/*----------------------- Parameters --------------------------*/
  parameters {
      vector[3] beta;                // intercept and covariate coefficients
      real<lower = 0, upper=1 > phi; // Auto-regressive parameter
      real<lower=0> sdp;             // Standard deviation of the process equation
      vector[N] mu;                  // underlying process mean
  }
 


/*----------------------- Model --------------------------*/
  model {
      // Prior distributions
      phi ~ beta(1,1);
      beta ~ normal(0, 5);
      sdp ~ normal(0, 1);
      
      // distribution of first state
      mu[1] ~ normal(P_obs[1], sdp);
      
      for(i in 2:N){
          mu[i] ~ normal(beta[1]*(1-phi) + beta[2] * light[i] + beta[3] * Q[i] + phi * P_obs[i-1], sdp);
          
          if(miss_vec[i]==1){
              target += normal_lpdf(P_obs[i] | mu[i], sdo[i]);   
          }
      }
 
      // likelihood
 
  }
 
  generated quantities {
      // Posterior predictive distribution
      real y_rep[N-1] = normal_rng(mu[2:N], sdo[2:N]);
      
  }
 
  