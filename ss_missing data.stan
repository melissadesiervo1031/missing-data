
    data {
    int<lower=1> N; // Number of observations
    vector[N] y_miss; // Response variable, including missing values
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
    real light[N]; //light data
    }
    
    parameters {
    vector[y_nMiss] y_imp;// Missing data
    vector[N] X; // latent state data
    real<lower = 0> sigma_proc; // process error
    real<lower = 0> sigma_obs; // observation error
    real<lower = 0, upper=1 > phi;  // auto-regressive parameter
    real b1; // light parameter 
    }
    
    transformed parameters { 
    vector[N] y;
    y=y_miss; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
    
    model {
        X[1]~normal(y[1],sigma_proc); //set initial state
    
    for(i in 1:N){
    y[i] ~ normal(X[i], sigma_obs); // observation model
    }
    
    for (t in 2:N){
    X[t] ~ normal(phi*X[t-1]+b1*light[t], sigma_proc); // process model with unknown process error //regression model with AR errors
    
    }  
    
    // error priors
    sigma_proc ~ normal(0, 1); 
    sigma_obs ~ normal(0, 2);
    
    // single parameters priors 
    b1~normal(0, 5);
   
    //Prior for AR coefficient needs to be reparameterized
    phi~beta(1,1);
    }
 
generated quantities {
 vector[N] x_rep; // replications from posterior predictive dist
 vector [N] y_rep;
 x_rep[1]=normal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 x_rep[t]=normal_rng(phi*X[t-1]+b1*light[t], sigma_proc);
 }
 for (t in 1:N) {
 y_rep[t]=normal_rng(x_rep[t], sigma_obs);
 }

  
 }
 
    
