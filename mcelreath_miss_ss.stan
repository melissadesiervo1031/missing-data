
 

/*----------------------- Functions --------------------------*/  
  functions{
    vector merge_missing( int[] y_index_mis, vector y_miss, vector y_imp) {
    int N = dims(y_miss)[1];
    int N_miss = dims(y_imp)[1];
    vector[N] merged;
    merged = y_miss;
    for (i in 1:N_miss)
        merged[y_index_mis[i] ] =y_imp[i];
    return merged;    
    }
  }
  
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
  }
  
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    vector<lower = 0>[y_nMiss] y_imp;// Missing data
    real<lower=0> sdp; // Standard deviation of the process equation
    //real<lower=0> sdo; // Standard deviation of the observation equation
    real b0;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
   // vector[N] z; // State time series
  }
  

  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
    //merge imputed and observed data
    vector[N] B_merge;
    B_merge= merge_missing(y_index_mis, to_vector(y_miss), y_imp);
   
    // Prior distributions
    //sdo ~ normal(0, 1);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
  
    
    // Distribution for the first state
  B_merge[1] ~ normal(z0, sdp);
   //y_miss[1]~ normal(B_merge[1],sdo);
   
    // Distributions for all other states
    for(t in 2:N){
       B_merge[t] ~ normal(b0+ B_merge[t-1]*phi, sdp);
       //y_miss[t] ~ normal(B_merge[t], sdo); // observation model with fixed observation error
    
    }
   
  }
  
  
    
