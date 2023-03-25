

/*----------------------- Data --------------------------*/
  data {
      int N;              // Length of state and observation time series
      vector[N] P_obs;    // observed GPP
      vector[N] light;    // covariate: light
      vector[N] Q;        // covariate: discharge
      vector[N] sdo;      // standard deviation of the daily GPP estimates (fixed)
      vector[N] miss_vec; // vector of 0's (missing) and 1's (data) describing
                             // locaions of missing GPP data
  }

/*----------------------- Parameters --------------------------*/
  parameters {
      vector[3] beta;                // intercept and covariate coefficients
      real<lower = 0, upper=1 > phi; // Auto-regressive parameter
      real<lower=0> sdp;             // Standard deviation of the process equation
      vector <lower = 0> [N] GPP;    // latent variable - true GPP
  }

/*----------------------- Model --------------------------*/
  model {
      // Prior distributions
      phi ~ beta(1,1);
      beta ~ normal(0, 5);
      sdp ~ normal(0, 1);
  
      // Distribution for the zero and first state
      GPP[1] ~ normal(P_obs[1], sdp);
  
      // Distributions for states
      for(i in 2:N){
          GPP[i] ~ normal(beta[1] + GPP[i-1]*phi + light[i]*beta[2] + Q[i]*beta[3], sdp); // process model with error
         
          if(miss_vec[i]==1){
              target += normal_lpdf(P_obs[i]| GPP[i], sdo[i]); // observation model with fixed observation error
          }
      }

  }
  
  generated quantities {
      real GPP_rep[N];
      real mu[N];
      mu[1] = normal_rng(P_obs[1], sdp);
  
      // Distributions for states
      for(i in 2:N){
          mu[i] = normal_rng(beta[1]+ mu[i-1]*phi + light[i]*beta[2] + Q[i]*beta[3], sdp); // process model with error
      }
      
      GPP_rep = normal_rng(mu, sdo); // observation model with fixed observation error

  }
    
    
