
  
  /*----------------------- Data --------------------------*/
    data {
      int N; // Number of data
      int K; //Number of predictors
      vector[N] y; // Deer population observations
      vector[N] snow; // snow data
      vector[N] mei; // PDO data
      real y0; // Initial data prior
    }
  
  /*----------------------- Parameters --------------------------*/
    
    parameters {
      real<lower=0> sdp; // Standard deviation of the process equation
      vector[K] beta;
      real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    }
  
  
  /*----------------------- Model --------------------------*/
    model {
      // Prior distributions
      sdp ~ normal(0, 1);
      phi ~ beta(1,1);
      beta[1] ~ normal(0,5);
      beta[2] ~ normal(0,5);
      beta[2] ~ normal(0,5);
      
      // Distribution for the first data
      y[1] ~ normal(y0, sdp);
      
      // Distributions for all other data
      for(t in 2:N){
        y[t] ~ normal(beta[1]+phi*y[t-1]+ beta[2]*snow+ beta[3]*mei, sdp);// process model with error
      }
    }
   

  
