# File descriptions

* `AR_ARIMA_STAN_comparison.R`: This script compares model fits to simulated data using different methods. The data are simulated from the model

  $$\mu_t = \boldsymbol x_t^\top \boldsymbol \beta + \phi \mu_{t-1}$$

  $$y_t \sim \mathcal{N}(\mu_t, \sigma^2)$$

  The custom Stan code can recovers the parameters from the model, while my reading of the code and the `brms` and `forecast::Arima` methods is that they fit different models, meaning they don't recover the parameters from the simulation. These methods fit the model

  $$y_t = \boldsymbol x_t^\top \boldsymbol \beta + \eta_{t}$$

  $$\eta_t = \phi \eta_{t-1} + \epsilon_t$$

  That is, the error term is assumed to follow an AR1 process with mean 0 and variance $\sigma^2 / (1 - \phi)$.

  *Code reviewed by DG, 3 Jul 2024*


* `Bias_check_expanded.R `: This script double checks whether the bias in the EM estimator diminishes with sample size (is it asymptotically unbiased). 
*This code isn't working right for me...chat with AMY. MD, 10 Jul 2024*

* `checking_ricker_DA_function.R `: This script double checks the Ricker data augmentation function with a small simulated dataset
  *Code reviewed by MD, 10 Jul 2024*
  
* `test_AR1_model_fit.R`: This script fits the model specified by the file `/GPP\ sim\ and\ real/Stan_code/AR1_light_Q_centered.stan` to 5 simulated datasets and one real dataset to check that the code in the Stan file is working as expected.

  - **This script does not run as is** because line 23 attempts to load data that we are no longer using and because line 42 uses an old version of the `makeMissing()` function with deprecated arguments.
  
  - Otherwise *gold standard* for level of comments and code readability.
  

