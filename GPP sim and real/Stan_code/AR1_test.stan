//

data {
    int<lower = 1> N;    // Length of state and observation time series
    vector[N] y;         // observations
    vector[N] light;     // covariate
    vector[N] igflow; // covariate
    vector[N] turb; // covariate
    
}
  
parameters {
    vector[4] beta;                // intercept and covariate coefficients
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real<lower=0> sigma;           // standard deviation of observation error
}
    
transformed parameters {
    vector[N] mu;   // underlying process mean
    vector[N] X;
    
    mu = beta[1] + beta[2] * light + beta[3] * igflow + beta[4] * turb;
    
    X = y - mu;
    
}

model {
  
    // Prior distributions
    phi ~ beta(1,1);
    beta ~ normal(0, 5);
    sigma ~ normal(0, 1);
 
    // likelihood
    for(t in 2:N){
        X[t] ~ normal(phi * X[t-1], sigma);
    }

    
}
  
