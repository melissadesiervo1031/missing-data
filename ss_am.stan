
 

/*----------------------- Data --------------------------*/
  data {
    int<lower=1> N; // Length of state and observation time series
    vector[N] y; // Observations
    real m0; //intital mean of prior state value
    real<lower=0> sdp0; //variance of prior distribution
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    real z0; //state [0];
    vector[N] z; //latent state variable
    real<lower =0> sdp; // Standard deviation of the process equation
    real<lower =0> sdo; // Standard deviation of the process equation
    real<lower = -1, upper=1 > phi; // Auto-regressive parameter
    real b0; 
    
    }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    phi ~ uniform(-1,1);
    b0 ~ normal(0, 5);
    sdp ~cauchy(0,1);
    sdo ~cauchy(0,1);
    
     // Observation error model       
   for(i in 1:N){
   y[i] ~ normal(z[i], sdo); // observation model with fixed observation error
   }
    
    // Distribution for the first state
    z0 ~ normal(m0,sdp0);
    z[1] ~ normal(z0,sdp);
   
    // Distributions for states
    for(i in 2:N){
       z[i] ~ normal(z[i-1]*phi+b0, sdp);// process model with error
        }
   
 
 }

