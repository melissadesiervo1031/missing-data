//

data {
    int<lower = 1> N;    // Length of state and observation time series
    vector[N] y;         // observations
    vector[N] light;     // covariate
    vector[N] discharge; // covariate
    
}
  
parameters {
    vector[3] beta;                // intercept and covariate coefficients
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real<lower=0> sigma;           // standard deviation of observation error
}
    
transformed parameters {
    vector[N] mu;   // underlying process mean

    mu = beta[1] + beta[2] * light + beta[3] * discharge;
    
}

model {
  
    // Prior distributions
    phi ~ beta(1,1);
    beta ~ normal(0, 5);
    sigma ~ normal(0, 1);
 
    // likelihood
    for(n in 2:N){
        y[n] ~ normal(mu + phi * y[n-1], sigma);
    }
    
}
  
