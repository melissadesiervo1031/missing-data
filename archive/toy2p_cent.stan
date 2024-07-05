
 
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    vector[N] light;
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
  }
  transformed data {
  vector[N] light_std;
  vector[N] y_std;
  light_std = (light - mean(light)) / sd(light);
  y_std = (y_miss - mean(y_miss)) / sd(y_miss);
}
  
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    vector<lower = 0>[y_nMiss] y_imp;// Missing data
    real<lower=0> sdp_std; // Standard deviation of the process equation
    real b0_std;
    real b1_std;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    
  }
  transformed parameters { 
    vector[N] y;
    y=y_std; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
    // Prior distributions
   
    sdp_std ~ normal(0, 1);
    phi ~ beta(1,1);
    b0_std ~ normal(0,5);
    b1_std ~ normal(0,5);
    
    // Distribution for the first state
    y_std[1] ~ normal(z0, sdp_std);
    
    // Distributions for all other states
    for(t in 2:N){
      y_std[t] ~ normal(b0_std+y_std[t-1]*phi+light_std[t]*b1_std, sdp_std);
    }
    
  }
  generated quantities {
 vector[N] y_rep; // replications from posterior predictive dist
 real b0;
 real b1;
 real<lower=0> sdp;
 b0 = sd(y_miss) * (b0_std - (b1_std * mean(light)/sd(light)))+mean(y_miss);
 b1 = b1_std * sd(y_miss)/ sd(light);
 sdp= sd(y_miss)* sdp_std;
 
 y_rep[1]=normal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 y_rep[t]=normal_rng(b0+y_miss[t-1]*phi+light[t]*b1, sdp);
 
 }
 
  }

    
