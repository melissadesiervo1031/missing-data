// 
//
/*----------------------- Data --------------------------*/
  data {
      int<lower = 1> N_obs;         // Number of observed data points
      int<lower = 0> N_mis;         // Number of missing observations
      int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs]; // indices of observations
      int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis]; // indices of missing obs
      vector[N_obs] P_obs;          // observed GPP
      vector[N_obs + N_mis] light; // covariate 
      vector[N_obs + N_mis] discharge;     // covariate 
  }
  
/*----------------------- Transformed data---------------*/
  transformed data{
      int<lower = 0> N = N_obs + N_mis;
  }
  
/*----------------------- Parameters --------------------*/
  parameters {
      vector[N_mis] P_mis;          // estimates of missing observations
      vector[3] beta;               // intercept and covariate coefficients
      real<lower = 0, upper=1> phi; // Auto-regressive parameter
      real<lower=0> sigma;          // Standard deviation of the random innovations
  }
  
/*----------------------- Transformed parameters---------*/
  transformed parameters {
      vector[N] y; 
      vector[N] y2;
      vector[N] mu;
      
      y[ii_obs] = P_obs;
      y[ii_mis] = P_mis;
      
      mu = beta[1] + beta[2]*light + beta[3]*discharge;
      
      y2 = y - mu;
  }
  
/*----------------------- Model -------------------------*/
  model {
      // Prior distributions
      phi ~ beta(1,1);
      beta ~ normal(0, 5);
      sigma ~ normal(0, 1);
      
      // likelihood
      for(i in 2:N){
          y2[i] ~ normal(phi * y2[i-1], sigma);
      }
 
  }
  