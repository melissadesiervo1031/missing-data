
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    vector[N] light;
    real z0;
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
    phi ~ beta(1,1);
    b0 ~ normal(0, 5);
    b1 ~ normal(0,5);
    sdp ~normal(0,1);
    sdo ~normal(0.01,0.001);
    
    // Observation error model       
 
   for(i in 1:N){
   y[i] ~ normal(z[i], sdo); // observation model with fixed observation error
   }
  
    
    // Distribution for the zero and first state
    
    z[1] ~ normal(z0,sdp);
   
    // Distributions for states
    for(i in 2:N){
       z[i] ~ normal(b0+ z[i-1]*phi+light[i]*b1, sdp);// process model with error
        }
   
 
 }
  
      
