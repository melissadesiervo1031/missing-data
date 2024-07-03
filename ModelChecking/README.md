# File descriptions

* `AR_ARIMA_STAN_comparison.R`: This script compares model fits to simulated data using different methods. The data are simulated from the model

  $$ \mu_t = x_t^\top \beta + \phi \mu_{t-1} $$

  $$ y_t \sim \mathcal{N}(\mu_t, \sigma^2) $$

  The custom Stan code can recover the parameters from the model, while `brms` and `forecast::Arima` fit different models, meaning they don't recover the parameters from the simulation. These methods fit the model

  $$ y_t = x_t^\top \beta + \eta_{t} $$

  $$ \eta_{t} = \phi \eta_{t-1} + \epsilon_t $$

  That is, the error term is assumed to follow an AR1 process with mean 0 and variance $\sigma^2 / (1 - \phi)$.

  *Code reviewed by DG, 3 Jul 2024*


* `file 2 `: An explanation of what this file does, and how it should be used in the analysis
