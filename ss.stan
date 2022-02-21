
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] z; // Observations
    real z0;
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    real<lower=0> sdp; // Standard deviation of the process equation
    
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
      }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
   
    // Distribution for the first state
    z[1] ~ normal(z0, sdp);
  
    // Distributions for all other states
    for(t in 2:N){
       z[t] ~ normal(phi*z[t-1], sdp);// process model with error
    }
  }
    
    
