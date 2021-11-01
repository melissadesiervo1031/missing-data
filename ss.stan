
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    vector[N] light;
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    vector[N] z; //latent state variable
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0; // intercept
    real b1; // light coefficient 
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
      }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    sdo ~ normal(0, 1);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    
    // Distribution for the first state
    z[1] ~ normal(y[1], 0.1);
  
    // Distributions for all other states
    for(t in 2:N){
       z[t] ~ normal(b0+ z[t-1]*phi+light[t]*b1, sdp);// process model with error
    }
       y ~ normal(exp(z), sdo); // observation model with fixed observation error
     }
  
  
    
