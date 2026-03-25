# Automatically tune the Factor Slice Sampler

## Description

Runs a three-stage auto-tuning pipeline for the Factor Slice Sampler:

1. ***Width tuning (round 1):*** `[tune_widths](tune_widths)` is called
    with coordinate-aligned factor directions to find good initial interval widths.
1. ***Factor estimation:*** `[estimate_factors](estimate_factors)` is called to
    rotate the slice directions to align with the posterior covariance structure.
1. ***Width re-tuning (round 2):*** `[tune_widths](tune_widths)` is called
    again along the estimated factor directions to refine the widths.

The returned `pars_afs` object can be passed directly to
`[factor_slice_sampler](factor_slice_sampler)` for production sampling.

## Usage

```r
auto_tune_sampler(theta_init, dat, lp, control = list())
```

## Arguments

* `theta_init`: Named numeric vector of initial parameter values.
* `dat`: Data list passed to `lp`.
* `lp`: Function to evaluate the log-joint probability. Must accept a parameter
vector as its first argument and the data list as its second.
* `control`: Named list of tuning controls. Recognised keys (all optional):
* `pars_afs`: Initial `pars_afs` list; `NULL` to start
    from scratch (default `NULL`).
* `tune_w_after`: Passed to `[tune_widths](tune_widths)` (default `1`).
* `ratio_target`: Target step-out ratio; passed to
    `[tune_widths](tune_widths)` (default `0.8`).
* `stop_after`: Maximum width-tuning steps; passed to
    `[tune_widths](tune_widths)` (default `2^10 + 1`).
* `burnin`: Burn-in iterations for factor estimation; passed to
    `[estimate_factors](estimate_factors)` (default `500`).

## Value

A list with two elements (the output of the final `[tune_widths](tune_widths)`  call):

* `pars_afs`: Fully tuned parameter list with optimised `Gamma`
      (factor directions), `width` (per-factor interval widths), and
      `theta_cov` (estimated posterior covariance).
* `theta`: Named numeric vector of the last sampled parameter values,
      ready to be used as the starting point for production sampling.

# Default control parameters for

## Description

Returns a named list of default values for the factor slice sampler when used
inside a Gibbs loop. Pass the output (or a modified copy) as the
`control_sampler` argument of `[fss_in_Gibbs_sample](fss_in_Gibbs_sample)`.

## Usage

```r
control_fssig()
```

## Value

A named list with the following elements:

* `iter`: Integer. Total number of Gibbs iterations to run
      (default `1000`).
* `time_out`: Integer. Maximum contraction attempts per factor
      direction before `[slice_proposals](slice_proposals)` throws an error
      (default `10000`).
* `max_steps_out`: Integer. Maximum expansion steps per direction
      in `[step_out](step_out)` (default `10000`).
* `track_interval_steps`: Logical. Whether to record stepping-out
      and contraction counts during sampling (default `TRUE`).

# Default control parameters for

## Description

Returns a named list of default values for the tuning pipeline. Pass the
output (or a modified copy) as the `control` argument of
`[auto_tune_sampler](auto_tune_sampler)`.

## Usage

```r
control_tuning()
```

## Value

A named list with the following elements:

* `pars_afs`: `NULL`; triggers initialisation from the
      identity matrix and small random widths inside `[tune_widths](tune_widths)`.
* `tune_w_after`: Integer. Number of sampler steps between each
      width adjustment (default `1`).
* `ratio_target`: Numeric. Target ratio of stepping-out to total
      steps (default `0.8`).
* `stop_after`: Integer. Maximum width-tuning steps per round
      (default `2^10 + 1`).
* `burnin`: Integer. Burn-in iterations used by
      `[estimate_factors](estimate_factors)` to estimate the posterior covariance
      (default `500`).

# Initialize theta from a list of y vectors

## Description

Fits `fit_ricker_cc` to each series, averages estimates, and
returns a named starting vector. For neg_binom adds lpsi.

## Usage

```r
.init_theta(y_list, fam = "poisson")
```

# Build a model matrix X for a single complete y vector

## Description

Build a model matrix X for a single complete y vector

## Usage

```r
.make_X(y)
```

# Compute the Hessian-based proposal covariance for a list of series

## Description

Uses the joint NLL across all series evaluated at `theta_init`.

## Usage

```r
.proposal_sigma(theta_init, y_list, fam = "poisson")
```

# Estimate factor directions from a burn-in run

## Description

Runs the factor slice sampler for `burnin` iterations, computes the
empirical covariance matrix of the resulting samples, and sets the factor
directions (`pars_afs$Gamma`) to the eigenvectors of that covariance
matrix. This aligns the slice directions with the principal axes of the
posterior, which is the core adaptation step of the Factor Slice Sampler
(Tibbits et al. 2013).

## Usage

```r
estimate_factors(theta_init, dat, lp, pars_afs, burnin = 500)
```

## Arguments

* `theta_init`: Named numeric vector of initial parameter values.
* `dat`: Data list passed to `lp`.
* `lp`: Function to evaluate the log-joint probability. Must accept a parameter
vector as its first argument and the data list as its second.
* `pars_afs`: List of current factor slice sampler parameters. Must contain
`Gamma` and `width`. Both are passed to `[factor_slice_sampler](factor_slice_sampler)`;
only `Gamma` and `theta_cov` are updated in the returned object.
* `burnin`: Number of iterations to run for covariance estimation. Defaults
to `500`.

## Value

A list with two elements:

* `pars_afs`: Updated parameter list with `Gamma` set to the
      eigenvectors of the sample covariance and `theta_cov` set to the
      sample covariance matrix itself.
* `theta`: Named numeric vector of the last sampled parameter values
      from the burn-in run, suitable for use as `theta_init` in a subsequent
      run.

# Run the Factor Slice Sampler (Algorithm 2 of Tibbits et al. 2013)

## Description

At each iteration a random slice height is drawn beneath the log-probability of
the current point, then each factor direction is updated in turn via
`[step_out](step_out)` (to bracket the slice) and `[slice_proposals](slice_proposals)`(to draw from within the bracket). The factor directions are taken from
`pars_afs$Gamma`, so the sampler reduces to univariate slice sampling along
the coordinate axes when `Gamma` is the identity matrix.

## Usage

```r
factor_slice_sampler(
  theta_init,
  dat,
  lp,
  pars_afs,
  iter = 1000,
  control = list()
)
```

## Arguments

* `theta_init`: Named numeric vector of initial parameter values.
* `dat`: Data list passed to `lp`.
* `lp`: Function to evaluate the log-joint probability. Must accept a parameter
vector as its first argument and the data list as its second.
* `pars_afs`: List of factor slice sampler parameters. Must contain:
* `Gamma`: Matrix whose columns are the factor directions.
* `width`: Numeric vector of interval widths, one per factor.
* `iter`: Number of MCMC iterations to run. Defaults to `1000`.
* `control`: Named list of algorithm controls. Recognised keys:
* `time_out`: Max proposals per contraction step (default `10000`).
* `max_steps_out`: Max expansion steps per direction (default `10000`).
* `track_interval_steps`: Logical; whether to record stepping-out and
    contraction counts (default `TRUE`).

## Value

A list with three elements:

* `theta`: Numeric matrix with `iter` rows and one column per
      parameter; row names omitted, column names inherited from `theta_init`.
* `out_steps`: Integer matrix (`iter` x `p`) of stepping-out
      counts per iteration and factor direction; `NULL` if
      `track_interval_steps` is `FALSE`.
* `in_steps`: Integer matrix (`iter` x `p`) of contraction
      counts per iteration and factor direction; `NULL` if
      `track_interval_steps` is `FALSE`.

# Function that will drop missing values and then fit model using ARIMA

## Description

Function that will drop missing values and then fit model using ARIMA

## Usage

```r
fit_arima_dropmissing(sim_list, sim_pars)
```

## Arguments

* `sim_list`: a list of GPP data from multiple simulations
* `sim_pars`: a list of simulation parameters corresponding to the GPP simulated values

## Value

List of the ARIMA parameter estimates, errors, and true values from the simulations.

## Examples

```r
gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_drop <- fit_arima_dropmissing(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)
```

# Function that will fit the model using ARIMA while dealing with missing

## Description

Function that will fit the model using ARIMA while dealing with missing 
values using the Kalman filter ARIMA

## Usage

```r
fit_arima_Kalman(sim_list, sim_pars)
```

## Arguments

* `sim_list`: a list of GPP data from multiple simulations
* `sim_pars`: a list of simulation parameters corresponding to the GPP simulated values

## Value

List of the ARIMA parameter estimates, errors, and true values from the simulations.

## Examples

```r
gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_kalman <- fit_arima_Kalman(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)
```

# Function that will fit the model using ARIMA while dealing with missing

## Description

Function that will fit the model using ARIMA while dealing with missing 
values using imputation via the Amelia package.

## Usage

```r
fit_arima_MI(sim_list, sim_pars, imputationsnum)
```

## Arguments

* `sim_list`: a list of GPP data from multiple simulations
* `sim_pars`: a list of simulation parameters corresponding to the GPP simulated values
* `imputationsnum`: the number of imputed datasets to create. This corresponds to the
`"m"` argument to the `"amelia"` function from the `"Amelia"` package.

## Value

List of the ARIMA parameter estimates, errors, and true values from the simulations.

## Examples

```r
gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_0.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_MI <- fit_arima_MI(GPP_sim_MAR$y,GPP_sim_MAR$sim_params, imputationsnum=5)
```

# fit_brms_model: use brms to fit a missing data model to simulated data

## Description

This function fits a brms missing-data AR1 model to a simulated data set 
multiple times, once for each different level of missingness in the data.
It returns a list of parameters from each model fit as well as the parameters
used to simulate the data.

## Usage

```r
fit_brms_model(sim_list, sim_pars, iter = 4000, include_missing = FALSE)
```

## Arguments

* `sim_list`: a list of datasets simulated from the same set of parameters
each with a different level of missingness.
* `sim_pars`: a list of parameters that was used to simulate the data 
including the value for the AR1 coefficient (phi), the value of the covariate
effects (beta) and the matrix of simulated covariates (X)
* `iter`: the number of iterations of the MCMC chain to run for each model. 
Defaults to 4000 because 2000 was giving warnings for the datasets with missing.
* `include_missing`: a boolean indicating if the parameter estimates for 
the missing data should be returned or not. Defaults to FALSE to match the 
arima function outputs.

## Value

A list of parameter estimates from fitting a stan AR1 model on each
simulated dataset.

## Examples

```r
gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

brms_fit <- fit_brms_model(GPP_sim_MAR$y, GPP_sim_MAR$sim_params, include_missing = FALSE)
```

# Fit Ricker count model with complete cases only

## Description

Fit Ricker count model with complete cases only

## Usage

```r
fit_ricker_cc(y, fam = "poisson", pro_conf = "none", off_patch = F)
```

## Arguments

* `y`: Vector of population counts through time. NA should be used for missing data.
* `fam`: Error distribution to use. Options include c("poisson", "neg_binom").
* `pro_conf`: what type of projection confidence to provide, defaults to "none", "sim" returns 1000 simulated estimates, "boot" returns 1000 bootstrapped estimates

## Value

List of intrinsic growth factor and intra-specific competitive effect estimates,
standard errors, and 95

## Examples

```r
y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
fit_ricker_cc(y)
```

# Fit a Bayesian Ricker population growth model with (potentially) missing data

## Description

Estimates the posterior distribution of the Ricker model parameters $r$and $\alpha$ (and dispersion $\psi$ for `fam = "neg_binom"`) using
a data-augmentation MCMC scheme. Starting values are obtained via empirical Bayes
(complete-case MLE), the Factor Slice Sampler is auto-tuned on the complete cases,
and then `chains` independent chains are run via
`[fss_in_Gibbs_sample](fss_in_Gibbs_sample)`, which alternates between Gibbs imputation of
missing counts and factor slice updates of the parameters.

## Usage

```r
fit_ricker_DA(
  y,
  fam = "poisson",
  chains = 4,
  samples = 1000,
  burnin = 5000,
  priors_list = list(m_r = 0, sd_r = 2.5, m_lalpha = -3, sd_lalpha = 1, m_lpsi = 0,
    sd_lpsi = 1),
  nthin = 5,
  return_y = FALSE
)

fit_ricker_DA(
  y,
  fam = "poisson",
  chains = 4,
  samples = 1000,
  burnin = 5000,
  priors_list = list(m_r = 0, sd_r = 2.5, m_lalpha = -3, sd_lalpha = 1, m_lpsi = 0,
    sd_lpsi = 1),
  nthin = 5,
  return_y = FALSE
)
```

## Arguments

* `y`: Either a single numeric vector of counts, or a ***list*** of
such vectors. `NA` marks missing observations. Each series must
begin with a non-NA value.
* `fam`: Character: `"poisson"` or `"neg_binom"`.
* `chains`: Number of parallel MCMC chains.
* `samples`: Number of posterior samples to keep per chain.
* `burnin`: Number of warm-up samples per chain.
* `priors_list`: Named list of prior hyperparameters:
* m_r, sd_r: Normal prior on r
* m_lalpha, sd_lalpha: Normal prior on log(alpha)
* m_lpsi, sd_lpsi: Normal prior on log(psi); only used for neg_binom
* `nthin`: Thinning interval.
* `return_y`: Logical. If `TRUE` and there is missingness, return
posterior summaries for imputed observations.
* `off_patch`: Logical. If `TRUE`, `y` is expected to be a data
frame with columns `yt` and `ytm1` (already in lagged-pair format,
e.g. for a multi-patch study). If `FALSE` (default), `y` is a plain
vector and the lagged-pair data frame is constructed internally.
* `control_sampler`: Named list of controls passed to
`[fss_in_Gibbs_sample](fss_in_Gibbs_sample)` via `[control_fssig](control_fssig)`. Use this
to override the number of sampling iterations, time-out limits, etc.
* `control_tuner`: Named list of controls passed to
`[auto_tune_sampler](auto_tune_sampler)` via `[control_tuning](control_tuning)`. Use this
to override the width-tuning schedule, factor-estimation burn-in, etc.

## Seealso

`[auto_tune_sampler](auto_tune_sampler)`, `[fss_in_Gibbs_sample](fss_in_Gibbs_sample)`,
`[control_tuning](control_tuning)`, `[control_fssig](control_fssig)`

## Value

On success, a named list with five elements:

* `estim`: Named numeric vector of posterior means. Parameters are
      returned on the natural scale ($\alpha$, not $\log\alpha$;
      $\psi$, not $\log\psi$). If `return_y = TRUE` and data are
      missing, imputed count means are appended as `y1`, `...`, `yn`.
* `rhat`: Named numeric vector of Gelman-Rubin $\hat{R}$
      convergence diagnostics (one per parameter), computed via
      `posterior::rhat`.
* `se`: Named numeric vector of posterior standard deviations.
* `lower`: Named numeric vector of 2.5% posterior quantiles.
* `upper`: Named numeric vector of 97.5% posterior quantiles.

On failure (population extinction, insufficient non-missing pairs, sampler
  getting stuck, etc.) a two-element list `list(NA, reason = "...")` is
  returned with a descriptive warning.
A list with elements `estim`, `rhat`, `se`,
`lower`, `upper`. If `return_y = TRUE` and data have
  missing values, imputed observation summaries are appended.

# Fit Ricker count model by naively dropping NAs

## Description

This function, rather than constructing a complete cases dataframe,
uses a naive dropping approach, thus violating the assumption of equal
spacing between observations in the time series.

## Usage

```r
fit_ricker_drop(y, fam = "poisson", pro_conf = "none", off_patch = F)
```

## Arguments

* `y`: Vector of population counts through time. NA should be used for missing data.
* `fam`: Error distribution to use. Options include c("poisson", "neg_binom").
* `pro_conf`: what type of projection confidence to provide, defaults to "none", "sim" returns 1000 simulated estimates, "boot" returns 1000 bootstrapped estimates

## Value

List of intrinsic growth factor and intra-specific competitive effect estimates,
standard errors, and 95

## Examples

```r
y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
fit_ricker_cc(y)
```

# Wrapper to fit a Ricker count model to data using the EM algorithm

## Description

This function is a wrapper to fit the stochastic Ricker model with
either Poisson or Negative Binomial error distribution and missing observations
encoded as NAs.

## Usage

```r
fit_ricker_EM(y, fam = "poisson", off_patch = FALSE, ...)
```

## Arguments

* `y`: Vector of population counts, with NA in the place of missing observations. If `off_patch == TRUE`,
this should be a dataframe with `yt` and `ytm1` in separate columns, pre-filtered such
that `ytm1` for a given time series or "patch" does not start with `NA`.
* `fam`: Error family. Can be either c("poisson", "neg_binom").
* `off_patch`: Logical. Set to `TRUE` when the data come from multiple replicate time series.
* `...`: Additional arguments passed to the EM algorithm, such as initial values 
(init_theta = c()) or maximum iterations before stopping (max_iter = 50).

## Value

List of intrinsic growth factor and intra-specific competitive effect estimates,
standard errors, and 95EM algorithm, so `se = NA; lower = NA; upper = NA`.

## Examples

```r
y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
fit_ricker_EM(y)
```

# Title Fit Ricker Model to population count data with Multiple Imputation using Amelia

## Description

This function

## Usage

```r
fit_ricker_MI(
  y,
  imputationsnum = 5,
  fam = "poisson",
  method = "dual",
  p2samelia = 1,
  ameliatimeout = 60,
  pro_conf = "none",
  off_patch = F
)
```

## Arguments

* `y`: vector of population time series data with NAs in any positions
* `imputationsnum`: number of imputations to average together, more will take longer but may have greater accuracy
* `fam`: either fit the models using poisson method ("poisson") or negative binomial method ("neg_binom")
* `method`: Each missing data point will be filled in twice by Amelia, once as y_t and once as y_t-1, should Amelia only predict y_t ("forward"), only predict y_t-1 ("backward"), use both and allow discrepancies ("dual", default), or use both and average out discrepancies ("averaging")
* `p2samelia`: an integer value to pass to amelia taking either 0 for no screen output, 1 for normal screen printing of iteration numbers, and 2 for detailed screen output
* `ameliatimeout`: a number in seconds indicating when to cut off Amelia's fitting algorithm and indicate convergence issues, default is 1 minute
* `pro_conf`: what type of projection confidence to provide, defaults to "none", "sim" returns 1000 simulated estimates, "boot" returns 1000 bootstrapped estimates

## Value

List of estimates, standard errors, and confidence intervals, or NA if an error occurred

## Examples

```r
y <- readRDS("data/missingDatasets/pois_sim_randMiss_B.rds")[[1]]$y[[10]]
fit_ricker_MI(y,ameliatimeout=10)
```

# Factor Slice Sampler within a Gibbs data-augmentation loop

## Description

Alternates between two steps at each iteration:

1. ***Gibbs step:*** `fill_rng` is called to draw the missing
    observations from their full conditional distribution given the current
    parameter values $\boldsymbol\theta$.
1. ***Factor slice step:*** `[factor_slice_sampler](factor_slice_sampler)` is called
    for a single iteration to draw a new $\boldsymbol\theta$ given the
    just-imputed complete data.

This function assumes that the factor slice sampler has already been tuned
(via `[auto_tune_sampler](auto_tune_sampler)`) and that a valid `pars_afs` object
is supplied.

## Usage

```r
fss_in_Gibbs_sample(
  theta_init,
  dat,
  lp,
  fill_rng,
  pars_afs,
  control_sampler = list()
)
```

## Arguments

* `theta_init`: Named numeric vector of initial parameter values.
* `dat`: Data list. Must contain a `y` element that is a data frame with
columns `yt` and `ytm1`, with `NA` in `yt` for any missing
time steps.
* `lp`: Function to evaluate the log-joint probability. Must accept a parameter
vector as its first argument and the data list as its second.
* `fill_rng`: Function that imputes missing observations. Must accept the current
parameter vector as its first argument and the data list as its second, and return
an updated version of `dat$y` with all `NA`s in `yt` filled in.
* `pars_afs`: Tuned factor slice sampler parameter list, as returned by
`[auto_tune_sampler](auto_tune_sampler)`. Must contain `Gamma` (factor direction
matrix) and `width` (per-factor interval widths).
* `control_sampler`: Named list of sampler controls. Unrecognised keys are
ignored; missing keys fall back to `[control_fssig](control_fssig)` defaults.
Recognised keys:
* `iter`: Total number of Gibbs iterations (default `1000`).
* `time_out`: Max contraction attempts per factor direction
    (default `10000`).
* `max_steps_out`: Max expansion steps per direction
    (default `10000`).

## Value

A list with three elements:

* `theta`: Numeric matrix with `iter` rows and one column per
      parameter; column names inherited from `theta_init`.
* `y`: Numeric matrix with `iter` rows and one column per time
      step, storing the full (imputed) `yt` vector at each iteration.
      Columns corresponding to observed time steps will be constant.
* `miss_ids`: Integer vector of row indices in `dat$y` where
      `yt` was `NA`, identifying which columns of `y` are imputed.

# Function that will introduce varying types and amounts of missing data in a time series

## Usage

```r
makeMissing(
  timeSeries,
  typeMissing,
  propMiss = NULL,
  autoCorr = NULL,
  type,
  ...
)
```

## Arguments

* `timeSeries`: a time series in vector format (a single vector, not a list)
* `typeMissing`: a character string indicating how missingness is to be generated 
("random" (with varying degrees of autocorrelation) or "minMax")
* `propMiss`: an optional vector argument that gives the proportion(s) of the data 
that should be missing; if not specified, then increasingly more missingness is 
introduced from 5\itemautoCorran optional argument between 0 and 1 (non-inclusive) giving the degree of
autocorrelation in missingness that the timeseries will have. (0 is no autocorrelation, 
.99 is maximum autocorrelation). Only used if typeMissing = "random". If a value isn't 
supplied, the function assumes no autocorrelation (0)\itemtypeprovides the 'type' of time series (i.e. "Poisson" or "Gaussian"), 
which determines how the function trims missing values if the minMax option is selectedList of time series with different proportions of missing dataFunction that will introduce varying types and amounts of missing data in a time seriesricker <- readRDS("./data/ricker_0miss_datasets.rds")
makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "random", propMiss = c(.5, .05), autoCorr = .003)
makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "random", propMiss = c(.5, .05))

# Metropolis-Hastings sampling of the posterior with block updates

## Description

Metropolis-Hastings sampling of the posterior with block updates
Metropolis-Hastings sampling of the posterior with block updates

## Usage

```r
MH_block_sample(dat, lp, burnin, iter, nthin = 1)

MH_block_sample(dat, lp, burnin, iter, nthin = 1)
```

## Arguments

* `dat`: List with elements `y_list` (list of complete count vectors),
plus prior hyperparameters `m_r`, `sd_r`, `m_lalpha`,
`sd_lalpha`, and (for neg_binom) `m_lpsi`, `sd_lpsi`.
* `lp`: Log-joint probability function: `lp(theta, dat)`.
* `burnin`: Number of warm-up samples.
* `iter`: Number of samples to keep.
* `nthin`: Thinning ratio.
* `q_rng`: Proposal random number generator. This function has one argument and takes in
the current value of $\boldsymbol \theta$.
* `q_lpdf`: Log proposal density. This function takes in the proposal for $\boldsymbol \theta$
as its first argument and the current value as its second.

## Value

A list with a matrix called `theta` (one row per sample and one column per parameter).
The second element in the list is a vector of whether the proposal was accepted in each iteration.
List with `theta` matrix and `acceptance` vector.

# Block Metropolis-Hastings within Gibbs sampling of the posterior for Data Augmentation

## Description

Block Metropolis-Hastings within Gibbs sampling of the posterior for Data Augmentation
Block Metropolis-Hastings within Gibbs sampling with Data Augmentation

## Usage

```r
MH_Gibbs_DA(dat, fill_rng, lp, burnin, iter, nthin = 1)

MH_Gibbs_DA(dat, fill_rng, lp, burnin, iter, nthin = 1)
```

## Arguments

* `dat`: List with `y_list` (list of count vectors, NAs allowed),
prior hyperparameters, and `fam`.
* `fill_rng`: Function to impute missing values: `fill_rng(theta, dat)`,
returns an updated `y_list` with NAs replaced.
* `lp`: Log-joint probability: `lp(theta, dat, y_list_full)`.
* `burnin`: Number of warm-up samples.
* `iter`: Number of samples to keep.
* `nthin`: Thinning ratio.
* `q_rng`: Proposal random number generator. This function has one argument and takes in
the current value of $\boldsymbol \theta$.
* `q_lpdf`: Log proposal density. This function takes in the proposal for $\boldsymbol \theta$
as its first argument and the current value as its second.

## Value

A list with a matrix called `theta` (one row per sample and one column per parameter). 
Column names are inherited from `theta_init`. The second element in the list is a matrix `y`with one row per iteration and a column for each observation of the response. The values in these columns
will be fixed for observed values of y, but will vary for unobserved values. The final element in the 
list is a vector of whether the proposal was accepted in each iteration.
List with `theta` matrix, `y_list` (list of matrices,
  one row per kept sample), and `acceptance` vector.

# Compute the mode of an empirical distribution

## Description

Compute the mode of an empirical distribution
Compute the mode of an empirical distribution

## Usage

```r
posterior_mode(x)

posterior_mode(x)
```

## Arguments

* `x`: Vector of samples from the target distribution.

## Value

The mode of the empirical distribution.
The mode of the empirical distribution.

# Negative log-likelihood for multiple independent Ricker time series

## Description

Assumes all series share the same parameter vector `theta`. Series are
treated as conditionally independent given `theta`, so their
log-likelihoods are summed.

## Usage

```r
ricker_count_neg_ll_multi(theta, y_list, X_list, fam = "poisson")

ricker_count_neg_ll_multi(theta, y_list, X_list, fam = "poisson")
```

## Arguments

* `theta`: Parameter vector. For Poisson: c(r, alpha). For neg_binom:
c(r, alpha, psi).
* `y_list`: List of population count vectors (NAs already filled in).
* `X_list`: List of model matrices corresponding to each series.
* `fam`: Error family: "poisson" or "neg_binom".

## Value

Scalar total negative log-likelihood summed across all series.

# Negative log-likelihood function for Ricker time series with

## Description

Negative log-likelihood function for Ricker time series with
conditionally Poisson- or negative-binomial-distributed observations
Negative log-likelihood function for Ricker time series with
conditionally Poisson- or negative-binomial-distributed observations

## Usage

```r
ricker_count_neg_ll(theta, y, X = NULL, fam = "poisson")

ricker_count_neg_ll(theta, y, X = NULL, fam = "poisson")
```

## Arguments

* `theta`: Generic parameter vector. For Poisson: length 2 (r, alpha).
For neg_binom: length 3 (r, alpha, psi) where psi is the size/dispersion
parameter on the natural scale (must be positive).
* `y`: The observed count data
* `X`: ($n$ x $P$) Model matrix. The first column should be all 1's and the second
should be y. Additional covariates can be added.
* `fam`: Error family: "poisson" or "neg_binom".

## Value

Scalar value of the negative log likelihood of `theta` given the data
Scalar value of the negative log likelihood of `theta` given the data

# Fitting the Ricker population model to count data

## Description

This function uses an Expectation Maximization (EM) approach to fit a Ricker population model
with Poisson or Negative-Binomial demographic stochasticity.

## Usage

```r
ricker_EM(
  y,
  init_theta,
  fam = "poisson",
  tol = 1e-05,
  max_iter = 50,
  off_patch = FALSE
)
```

## Arguments

* `y`: Vector of population counts with NAs placed in positions in which there are missing
observations (if any).
* `init_theta`: Initial values for the parameter vector.
* `fam`: Family of the error distribution. Can be either `"poisson"` or `"neg_binom"`
* `tol`: Tolerance for convergence.
* `max_iter`: Maximum number of iterations to run the EM algorithm before stopping.

## Value

List with estimates and conditional standard errors for `theta`, estimates of the
full data vector with latent values filled in (`z_hat`), and a convergence code. 0 means 
the algorithm converged before being stopped.

# Simulate population abundance from a stochastic Ricker model

## Description

Simulate population abundance from a stochastic Ricker model

## Usage

```r
ricker_sim(n, r, alpha, N0, err_fam = "poisson", psi = NULL)
```

## Arguments

* `n`: Length of desired time series
* `r`: Log intrinsic growth factor
* `alpha`: Intra-specific density-dependence
* `N0`: Initial abundance
* `err_fam`: Error family. Either "poisson" or "neg_binom". If "neg_binom" is 
specified, then a dispersion parameter `psi` must be included.
* `psi`: Dispersion parameter for negative binomial error distribution.

# One step in the deterministic Ricker population process

## Description

One step in the deterministic Ricker population process
One step in the deterministic Ricker population process

## Usage

```r
ricker_step(theta, Nt)

ricker_step(theta, Nt)
```

## Arguments

* `theta`: Parameter vector with intrinsic growth rate $r$ as the first element
and negative density dependence $alpha$ as the second parameter.
NOTE that $alpha$ should be negative to induce negative density dependence.
* `Nt`: The population size at the current time.

## Value

Population size at the next time step.
Population size at the next time step.

# Draw a proposal from a slice interval via shrinkage (contraction)

## Description

Samples uniformly from the current interval `A_curr` along factor direction
`k` and contracts the interval toward the current point whenever the proposal
falls below the slice height `h`, following the shrinkage procedure of Neal (2003).

## Usage

```r
slice_proposals(k, A_curr, theta_curr, pars_afs, h, lp, dat, time_out = 1e+05)
```

## Arguments

* `k`: Integer index of the factor direction to sample along.
* `A_curr`: Numeric vector of length 2 giving the current lower and upper bounds
of the slice interval on the standardised factor scale.
* `theta_curr`: Numeric vector of current parameter values.
* `pars_afs`: List of factor slice sampler parameters. Must contain
`Gamma`, a matrix whose columns are the factor directions.
* `h`: Numeric scalar. Log-probability threshold defining the current slice height.
* `lp`: Function to evaluate the log-joint probability. Must accept a parameter
vector as its first argument and the data list as its second.
* `dat`: Data list passed to `lp`.
* `time_out`: Maximum number of proposals before an error is thrown. Defaults
to `100000`.

## Value

A list with two elements:

* `theta`: Numeric vector of the accepted proposal in parameter space.
* `contraction_count`: Integer count of how many times the interval
      was contracted before a valid proposal was found.

# Expand a slice interval along a factor direction using the stepping-out procedure

## Description

Starting from a randomly placed initial interval of width `pars_afs$width[k]`,
this function expands the lower and upper bounds until both lie below the
log-probability threshold `h`, implementing the stepping-out procedure of
Neal (2003) projected onto factor direction `k`.

## Usage

```r
step_out(k, theta_curr, pars_afs, h, lp, dat, max_steps_out = 10000)
```

## Arguments

* `k`: Integer index of the factor direction to step out along.
* `theta_curr`: Numeric vector giving the current parameter values.
* `pars_afs`: List of factor slice sampler parameters. Must contain:
* `width`: Numeric vector of interval widths, one per factor direction.
* `Gamma`: Matrix whose columns are the factor directions (eigenvectors).
* `h`: Numeric scalar. Log-probability threshold defining the current slice height.
* `lp`: Function to evaluate the log-joint probability. Must accept a parameter
vector as its first argument and the data list as its second.
* `dat`: Data list passed to `lp`.
* `max_steps_out`: Maximum number of expansion steps allowed in each direction
before giving up. Defaults to `10000`.

## Value

A list with two elements:

* `int`: Numeric vector of length 2 giving the lower and upper bounds
      of the expanded interval on the standardised factor scale.
* `n_expand`: Total number of expansion steps taken across both
      directions (lower and upper combined); set to 1 if no expansion was needed.

# Tune factor slice sampler interval widths

## Description

Runs the factor slice sampler and iteratively adjusts the interval width for
each factor direction so that the ratio of stepping-out steps to total steps
(stepping-out + contraction) converges toward `ratio_target`. Widths are
scaled up when the ratio exceeds the target and scaled down when it falls below.
Tuning is declared converged for a given direction once the multiplicative
adjustment falls within 10% of 1.

## Usage

```r
tune_widths(
  theta_init,
  dat,
  lp,
  pars_afs = NULL,
  tune_w_after = 1,
  ratio_target = 0.8,
  stop_after = 2^10 + 1
)
```

## Arguments

* `theta_init`: Named numeric vector of initial parameter values.
* `dat`: Data list passed to `lp`.
* `lp`: Function to evaluate the log-joint probability. Must accept a parameter
vector as its first argument and the data list as its second.
* `pars_afs`: List of factor slice sampler parameters, or `NULL` to
initialise with identity factor directions and small random widths. Must contain
`Gamma` (factor direction matrix) and `width` (per-factor widths)
when not `NULL`.
* `tune_w_after`: Number of sampler steps between width adjustments. Doubles
automatically once all step-out rates drop below 10 per adjustment window.
Defaults to `1`.
* `ratio_target`: Target ratio of stepping-out steps to total steps (stepping-out
+ contraction). Defaults to `0.8`.
* `stop_after`: Maximum total sampler steps before halting, regardless of
convergence. Defaults to `2^10 + 1`.

## Value

A list with two elements:

* `pars_afs`: Updated parameter list with tuned `width` values.
* `theta`: Named numeric vector of the last sampled parameter values,
      suitable for use as `theta_init` in a subsequent run.

