

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] logP_obs; // log observations
    vector[N] light;
    real logP_0; //initial underlying state
  }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    vector <lower = 0> [N] logP; //log latent state variable
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0; // intercept
    real b1; // light coefficient 
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
  }

/*----------------------- Model --------------------------*/

  model {
    // Prior distributions
    phi ~ beta(1,1);
    b0 ~ normal(0, 5);
    b1 ~ normal(0, 5);
    sdp ~ normal(0,1);
    sdo ~ normal(0.05,0.01);
    
    // Distribution for the zero and first state
    logP[1] ~ normal(logP_0,sdp);
    
    // Distributions for states
    for(i in 2:N){
       logP[i] ~ normal(b0+ logP[i-1]*phi + light[i]*b1, sdp);// process model with error
        }
   
    // Observation error model       
    logP_obs ~ normal(logP, sdo); 
  }
