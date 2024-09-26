

/*----------------------- Data --------------------------*/
  data {
    int N;           // Length of state and observation time series
    vector[N] P_obs; // observed GPP
    vector[N] light; // covariate
    vector[N] sdo;   // standard deviation of the observations (fixed)
  }

/*----------------------- Parameters --------------------------*/

  parameters {
    vector[2] beta;                // intercept and slope
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real<lower=0> sdp;             // Standard deviation of the process equation
    vector <lower = 0> [N] GPP;    // latent variable - true GPP
  }

/*----------------------- Model --------------------------*/

  model {
    // Prior distributions
    // phi ~ beta(1,1);
    beta ~ normal(0, 5);
    sdp ~ normal(0,1);

    // Distribution for the zero and first state
    GPP[1] ~ normal(P_obs[1], sdp);

    // Distributions for states
    for(i in 2:N){
       GPP[i] ~ normal(beta[1]*(1-phi) + GPP[i-1]*phi + light[i]*beta[2], sdp); // process model with error
    }

    // Observation error model
    P_obs ~ normal(GPP, sdo);

  }
  
  generated quantities {
    real GPP_rep[N];
    real mu[N];
    mu[1] = normal_rng(P_obs[1], sdp);

    // Distributions for states
    for(i in 2:N){
       mu[i] = normal_rng(beta[1]*(1-phi) + mu[i-1]*phi + light[i]*beta[2], sdp); // process model with error
    }

    // Observation error model
    GPP_rep = normal_rng(mu, sdo);
  
  }
    
    
