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
  
  *Code reviewed by DG, 15 July 2024*

* `fit_GPP_datasets.R`: This script fits a linear model with two covariates and an AR1 correlation structure in the errors (using `brms`) to multiple GPP datasets. Datasets from Au Sable River and Badger Mill Creek were selected, and this script cleans those data for later use.

  - Not sure this script belongs in this folder or if it could be archived.
  
  *Code reviewed by D.G. 17 July 2024*
  
* `gauss_amelia_testing.R`: This script tests multiple options available in the `Amelia` package for imputing missing data in time series. No conclusions are drawn in the comments.

  - This script uses a single dataset to test these four different options, so should maybe be expanded to test these options over many datasets? Or use author recommendations if available.
  
  - There are some hard-coded value substitutions in this script that I don't fully understand. These will likely break if the `gauss_sim_randMiss_A.rds` dataset changes.
  
  *Code reviewed by D.G. 17 July 2024*
  
* `generate_subset_NWIS_data_low_missingness.R`: This script looks like it loads some data and cleans a subset of it for potential use as an example dataset in the manuscript.

  - **This script does not run** because it uses absolute filepaths to a dataset that is not available in the project repo.
  
  - This script doesn't seem to belong in this directory.
  
  *Code reviewed by D.G. 17 July 2024*
  
* `plot_gauss_real_model_results.R`: This is a short plotting script for results.

  - **This script does not run**, likely due to reorganization of the repo and changes to file paths.
  
  *Code reviewed by D.G. 17 July 2024*
  

  

