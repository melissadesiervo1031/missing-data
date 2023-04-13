//

data {
    int<lower = 1> N;    // Length of state and observation time series
    int<lower = 0> K;    // Number of covariates
    vector[N] y;         // observed GPP
    matrix[N, K + 1] X;    // covariate matrix
}


parameters {
    vector[K+1] beta;                // intercept and covariate coefficients
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real<lower=0> sdo;             // standard deviation of observation error
    // real<lower=0> sdp;             // Standard deviation of the process equation
}


transformed parameters {
    vector[N] mu;   // underlying process mean
 
    mu[1] = y[1];
 
    for(i in 2:N){
        mu[i] = X[i, ] * beta + y[i-1] * phi;
    }
}

model {
    // Prior distributions
    // phi ~ beta(1,1);
    // beta ~ normal(0, 5);
    // sdp ~ normal(0, 1);
 
    // likelihood
    y[2:N] ~ normal(mu[2:N], sdo);   
 
}
 
generated quantities {
    // Posterior predictive distribution
    real y_rep[N-1] = normal_rng(mu[2:N], sdo);
}
  