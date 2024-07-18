# File descriptions
Each script within this repository begins with a header detailing the function purpose, inputs, outputs, and examples, so the descriptions provided in this README do not attempt to duplicate those details here.

# Ricker functions
* `ricker_sim.R`: This function simulates a stochastic population abundance time series according to the Ricker model. First the expected population size ($$N'_t$$) is caluclated according to $$N'_t = N_{t-1} * exp(r-\alpha*N_{t-1})$$. Then, the actual population size ($$N_t$$) is drawn from either a Poisson distribution or a Negative Binomial. 

* `ricker_MI_function.R `: This function will fill in missing data with multiple imputation using Amelia and then fit the Ricker Model using Poisson or Negative Binomial error distribution

* `ricker_drop_function.R`: This script contains two functions. The first (`fit_ricker_cc`) uses the complete case method to analyze a Ricker time series with missing data. The second function (`fit_ricker_drop`) uses a naive data deletion approach instead.

* `ricker_count_MCMC.R`: This script contains four functions all related to fitting a Ricker time series with the Data Augmentation (DA) approach. 
  `posterior_mode`: Returns the mode of a parameter distribution
  `MH_block_sample`: A custom function to implement the Metropolis-Hastings algorithm for sampling from the posterior distribution with no missing data
  `MH_Gibbs_DA`: Modified Metropolis-Hastings algorithm within Gibbs Sampling to sample from the posterior under DA
  `fit_ricker_DA`: A single function to use a Bayesian approach to fit a Ricker time series. The script checks for missing data and will then implement the appropriate sampling algorithm. Importantly, this function relies on `ricker_drop_function.R`, `ricker_count_EM.R`,  `ricker_count_likelihood functions.R`, and the other functions defined with this file. Also, it is intended to use the `parallel` package.

* `ricker_count_likelihood functions.R`: This file contains two functions useful for the calculation of the log likelihood of a Ricker population model.
  `ricker_step`: implements a single, deterministic step according to Ricker population growth
  `ricker_count_neg_ll`: calculates the negative log likelihood of a given parameter combination for the Ricker model and assuming either a Poisson or Negative Binomial error distribution

* `ricker_count_EM.R`: This script contains two functions to fit a Ricker model to discrete time series data with missing values via Expectation Maximization (EM).
  `ricker_EM`: Function implementing the EM algorithm for fitting data to the Ricker model
  `fit_ricker_EM`: Wrapper function to implement `ricker_EM` and return a list of parameter estimates

# Arima functions
* `Arima_drop_function.R`: This script contains the `fit_arima_dropmissing` function which drops missing values before fitting the model using ARIMA

* `Arima_Kalman_function.R`: This script contains the `fit_arima_Kalman` function which uses the Kalman filter to fit the model using ARIMA

* `Arima_MI_function.R`: This script contains the function `fit_arima_MI` which uses multiple imputation (as implemented via the amelia package) to fill in missing values and then fit the model using Arima

# Other functions
* `missing_data_functions.R`:

* `model_fitting.R`:

