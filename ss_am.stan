
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    real z0; //intital state value
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    vector[N] z; //latent state variable
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real b0; 
    
    }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    sdo ~ normal(0, 1);
    sdp ~ normal(0, 5);
    phi ~ beta(1,1);
    b0 ~ normal(0, 5);
    
    
    // Distribution for the first state
    z[1] ~ normal(z0,sdp);
  
    // Distributions for all other states
    for(i in 2:N){
       z[i] ~ normal(z[i-1]*phi+b0, sdp);// process model with error
    
       y[i] ~ normal(z[i], sdo); // observation model with fixed observation error
    }

 }

