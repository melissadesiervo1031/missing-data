//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int y[N];
}

parameters {
  real <lower=0> rmax;
  real <lower=0> alpha;
 }


model {
  
  vector[N] Nmod; 
  
  for ( i in 2:N)
  y[i]~ poisson(y[i-1]*exp(rmax-y[i-1]*alpha));
  
  alpha~normal(0,0.05);
  rmax~normal (0,5);
    
}

