
    data {
    int N; #Number of observations
    vector[N] temp; #Response variable, including missing values
    int Temp_nMiss;
    int Temp_index_mis[Temp_nMiss];
    }
    
    parameters {
    vector[Temp_nMiss] temp_imp;//Missing data
    real<lower = 0> sigma; //  standard deviation
     real alpha; // Constant
    real beta;  // AR
    }
    
    transformed parameters { 
    vector[N] y;
    y = temp; 
    y[Temp_index_mis] =temp_imp;
    } 
    
    model {
    for (n in 2:N){
    y[n] ~ normal(alpha+beta*y[n-1], sigma);
      }
    sigma ~ normal(1, 1);
    }
    
    
