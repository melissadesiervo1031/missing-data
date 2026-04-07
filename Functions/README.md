# Table of Contents

*Documentation written with help from Claude Code*

**MCMC & Sampling**
- [`auto_tune_sampler`](#automatically-tune-the-factor-slice-sampler)
- [`control_fssig`](#default-control-parameters-for)
- [`control_tuning`](#default-control-parameters-for-1)
- [`estimate_factors`](#estimate-factor-directions-from-a-burn-in-run)
- [`factor_slice_sampler`](#run-the-factor-slice-sampler-algorithm-2-of-tibbits-et-al-2013)
- [`fss_in_Gibbs_sample`](#factor-slice-sampler-within-a-gibbs-data-augmentation-loop)
- [`MH_block_sample`](#metropolis-hastings-sampling-of-the-posterior-with-block-updates)
- [`MH_Gibbs_DA`](#block-metropolis-hastings-within-gibbs-sampling-of-the-posterior-for-data-augmentation)
- [`slice_proposals`](#draw-a-proposal-from-a-slice-interval-via-shrinkage-contraction)
- [`step_out`](#expand-a-slice-interval-along-a-factor-direction-using-the-stepping-out-procedure)
- [`tune_widths`](#tune-factor-slice-sampler-interval-widths)

**Ricker Model — Fitting**
- [`fit_ricker_cc`](#fit-ricker-count-model-with-complete-cases-only)
- [`fit_ricker_DA`](#fit-a-bayesian-ricker-population-growth-model-with-potentially-missing-data)
- [`fit_ricker_drop`](#fit-ricker-count-model-by-naively-dropping-nas)
- [`fit_ricker_EM`](#fit-a-ricker-population-model-to-count-data-using-the-em-algorithm)
- [`fit_ricker_MI`](#title-fit-ricker-model-to-population-count-data-with-multiple-imputation-using-amelia)
- [`ricker_EM`](#fitting-the-ricker-population-model-to-count-data)

**Ricker Model — Likelihood & Helpers**
- [`ricker_count_nb_fit`](#fit-ricker-count-model-with-negative-binomial-errors-via-alternating-optimization)
- [`ricker_count_neg_ll`](#negative-log-likelihood-function-for-ricker-time-series-with)
- [`ricker_count_neg_ll_cnstr`](#constrained-negative-log-likelihood-for-a-ricker-count-time-series)
- [`ricker_count_neg_ll_multi`](#negative-log-likelihood-for-multiple-independent-ricker-time-series)
- [`ricker_count_pois_fit`](#fit-ricker-count-model-with-poisson-errors)
- [`ricker_sim`](#simulate-population-abundance-from-a-stochastic-ricker-model)
- [`ricker_step`](#one-step-in-the-deterministic-ricker-population-process)
- [`zt_neg_binom_rng`](#draw-samples-from-a-zero-truncated-negative-binomial-distribution)
- [`zt_poisson_rng`](#draw-samples-from-a-zero-truncated-poisson-distribution)

**Gaussian (ARIMA) Fitting**
- [`fit_arima_dropmissing`](#function-that-will-drop-missing-values-and-then-fit-model-using-arima)
- [`fit_arima_Kalman`](#function-that-will-fit-the-model-using-arima-while-dealing-with-missing)
- [`fit_arima_MI`](#function-that-will-fit-the-model-using-arima-while-dealing-with-missing-1)
- [`fit_brms_model`](#fit_brms_model-use-brms-to-fit-a-missing-data-model-to-simulated-data)

**Utilities**
- [`makeMissing`](#function-that-will-introduce-varying-types-and-amounts-of-missing-data-in-a-time-series)
- [`posterior_mode`](#compute-the-mode-of-an-empirical-distribution)
- [`.init_theta`](#initialize-theta-from-a-list-of-y-vectors)
- [`.make_X`](#build-a-model-matrix-x-for-a-single-complete-y-vector)
- [`.proposal_sigma`](#compute-the-hessian-based-proposal-covariance-for-a-list-of-series)

---

# Automatically tune the Factor Slice Sampler

## Description

Runs a three-stage auto-tuning pipeline for the Factor Slice Sampler:

1. ***Width tuning (round 1):*** [`tune_widths`](#tune_widths) is called
    with coordinate-aligned factor directions to find good initial interval widths.
2. ***Factor estimation:*** [`estimate_factors`](#estimate_factors) is called to
    rotate the slice directions to align with the posterior covariance structure.
3. ***Width re-tuning (round 2):*** [`tune_widths`](#tune_widths) is called
    again along the estimated factor directions to refine the widths.

The returned `pars_afs` object can be passed directly to
[`factor_slice_sampler`](#factor_slice_sampler) for production sampling.

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
* `tune_w_after`: Passed to [`tune_widths`](#tune_widths) (default `1`).
* `ratio_target`: Target step-out ratio; passed to
    [`tune_widths`](#tune_widths) (default `0.8`).
* `stop_after`: Maximum width-tuning steps; passed to
    [`tune_widths`](#tune_widths) (default `2^10 + 1`).
* `burnin`: Burn-in iterations for factor estimation; passed to
    [`estimate_factors`](#estimate_factors) (default `500`).

## Value

A list with two elements (the output of the final [`tune_widths`](#tune_widths)  call):

* `pars_afs`: Fully tuned parameter list with optimised `Gamma`
      (factor directions), `width` (per-factor interval widths), and
      `theta_cov` (estimated posterior covariance).
* `theta`: Named numeric vector of the last sampled parameter values,
      ready to be used as the starting point for production sampling.

# Default control parameters for

## Description

Returns a named list of default values for the factor slice sampler when used
inside a Gibbs loop. Pass the output (or a modified copy) as the
`control_sampler` argument of [`fss_in_Gibbs_sample`](#fss_in_gibbs_sample).

## Usage

```r
control_fssig()
```

## Value

A named list with the following elements:

* `iter`: Integer. Total number of Gibbs iterations to run
      (default `1000`).
* `time_out`: Integer. Maximum contraction attempts per factor
      direction before [`slice_proposals`](#slice_proposals) throws an error
      (default `10000`).
* `max_steps_out`: Integer. Maximum expansion steps per direction
      in [`step_out`](#step_out) (default `10000`).
* `track_interval_steps`: Logical. Whether to record stepping-out
      and contraction counts during sampling (default `TRUE`).

# Default control parameters for

## Description

Returns a named list of default values for the tuning pipeline. Pass the
output (or a modified copy) as the `control` argument of
[`auto_tune_sampler`](#auto_tune_sampler).

## Usage

```r
control_tuning()
```

## Value

A named list with the following elements:

* `pars_afs`: `NULL`; triggers initialisation from the
      identity matrix and small random widths inside [`tune_widths`](#tune_widths).
* `tune_w_after`: Integer. Number of sampler steps between each
      width adjustment (default `1`).
* `ratio_target`: Numeric. Target ratio of stepping-out to total
      steps (default `0.8`).
* `stop_after`: Integer. Maximum width-tuning steps per round
      (default `2^10 + 1`).
* `burnin`: Integer. Burn-in iterations used by
      [`estimate_factors`](#estimate_factors) to estimate the posterior covariance
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
`Gamma` and `width`. Both are passed to [`factor_slice_sampler`](#factor_slice_sampler);
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
[`step_out`](#step_out) (to bracket the slice) and [`slice_proposals`](#slice_proposals)(to draw from within the bracket). The factor directions are taken from
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
[`fss_in_Gibbs_sample`](#fss_in_gibbs_sample), which alternates between Gibbs imputation of
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
[`fss_in_Gibbs_sample`](#fss_in_gibbs_sample) via [`control_fssig`](#control_fssig). Use this
to override the number of sampling iterations, time-out limits, etc.
* `control_tuner`: Named list of controls passed to
[`auto_tune_sampler`](#auto_tune_sampler) via [`control_tuning`](#control_tuning). Use this
to override the width-tuning schedule, factor-estimation burn-in, etc.

## Seealso

[`auto_tune_sampler`](#auto_tune_sampler), [`fss_in_Gibbs_sample`](#fss_in_gibbs_sample),
[`control_tuning`](#control_tuning), [`control_fssig`](#control_fssig)

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
fit_ricker_drop(
  y,
  fam = "poisson",
  pro_conf = "none",
  off_patch = F,
  patch_col = "patch"
)
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

# Fit a Ricker population model to count data using the EM algorithm

## Description

Fits the stochastic Ricker model to a single time series of population counts,
optionally with missing observations. Missing values are handled via Expectation
Maximization (EM): at each E-step, missing counts are replaced by their expected
value under the current parameter estimates; at each M-step, parameters are
re-estimated by maximum likelihood on the filled-in series. Supports Poisson
and Negative Binomial error distributions.

## Usage

```r
fit_ricker_EM(y, fam = "poisson", off_patch = FALSE, ...)
```

## Arguments

* `y`: Either (1) a numeric vector of population counts with `NA` marking
missing observations, or (2) if `off_patch = TRUE`, a data frame with columns
`yt` (count at time $t$) and `ytm1` (count at time $t-1$),
already formatted into lag pairs and filtered so that no series begins with
`ytm1 = NA`. Leading `NA`s in a vector input are stripped with a warning.
* `fam`: Error distribution family. One of `"poisson"` (default) or
`"neg_binom"`.
* `off_patch`: Logical. Set to `TRUE` when `y` has already been
formatted into lag pairs (i.e., the `yt`/`ytm1` data frame format).
Use this when fitting to data from multiple replicate time series that have
been row-bound into a single data frame prior to calling this function.
Default is `FALSE`.
* `...`: Additional arguments passed to the internal `ricker_EM` algorithm:
* `init_theta`: Named numeric vector of starting values. For Poisson,
    `c(r = 0.5, lalpha = log(0.01))`; for Negative Binomial, also includes
    `lpsi = log(5)`. Note that `alpha` is parameterized on the log scale
    internally as `lalpha = log(alpha)`, and similarly for the NB dispersion
    parameter `psi`.
* `tol`: Convergence tolerance on the change in negative log-likelihood
    between successive EM iterations. Default `1e-5`.
* `max_iter`: Maximum number of EM iterations before stopping.
    Default `50`.

## Seealso

[`ricker_EM`](#ricker_em) for the underlying EM implementation,
[`fit_ricker_DA`](#fit_ricker_da) for Bayesian estimation with posterior uncertainty,
[`fit_ricker_cc`](#fit_ricker_cc) for complete-case (listwise deletion) MLE.

## Value

A list with the following elements:

* `estim`: Named numeric vector of parameter estimates on the natural
      scale: `r` (intrinsic growth rate) and `alpha` (intraspecific
      competition coefficient, positive). For `fam = "neg_binom"`, also
      includes `psi` (NB dispersion parameter).
* `se`: Always `NA`. The EM algorithm does not directly yield
      standard errors; use `fit_ricker_DA` for posterior uncertainty
      quantification.
* `lower`: Always `NA` (see `se`).
* `upper`: Always `NA` (see `se`).
* `convergence`: Convergence code. `0` indicates the algorithm
      converged (NLL change fell below `tol`) before reaching `max_iter`;
      `1` indicates the algorithm was stopped at `max_iter` without
      converging. In the latter case a warning is issued.

Returns a list with `NA` as its first element (and a named `cause` or
`reason` element) if the series contains zeros (population extinction),
`NaN`, `Inf`, or if the M-step optimizer fails.

## Examples

```r
y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
fit_ricker_EM(y)

# Negative Binomial with custom starting values and more iterations
fit_ricker_EM(y, fam = "neg_binom", init_theta = c(r = 1, lalpha = log(0.05), lpsi = log(10)),
              max_iter = 100)
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
1. ***Factor slice step:*** [`factor_slice_sampler`](#factor_slice_sampler) is called
    for a single iteration to draw a new $\boldsymbol\theta$ given the
    just-imputed complete data.

This function assumes that the factor slice sampler has already been tuned
(via [`auto_tune_sampler`](#auto_tune_sampler)) and that a valid `pars_afs` object
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
[`auto_tune_sampler`](#auto_tune_sampler). Must contain `Gamma` (factor direction
matrix) and `width` (per-factor interval widths).
* `control_sampler`: Named list of sampler controls. Unrecognised keys are
ignored; missing keys fall back to [`control_fssig`](#control_fssig) defaults.
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

# Fit Ricker count model with negative binomial errors via alternating optimization

## Description

Estimates Ricker population model parameters under a negative binomial
observation model using a two-step alternating optimization scheme.
In step 1, the intrinsic growth rate (`r`) and log-transformed
intra-specific competition coefficient (`lalpha`) are optimized
conditional on the current overdispersion estimate. In step 2, the
log-overdispersion parameter (`lpsi`) is optimized conditional on
the current predicted means. The two steps alternate until parameter
estimates converge or the iteration limit is reached.

## Usage

```r
ricker_count_nb_fit(theta, y, tol = 1e-05, max_iter = 500)
```

## Arguments

* `theta`: Named numeric vector of initial parameter values. Must contain
elements named `"r"` (intrinsic growth rate), `"lalpha"`
(log of the intra-specific competition coefficient), and `"lpsi"`
(log of the negative binomial overdispersion parameter).
* `y`: Either a numeric vector of population counts through time, or a
data frame with columns `yt` (count at time $t$) and
`ytm1` (count at time $t-1$). If a vector is supplied it is
automatically converted to the lagged data frame format.
* `tol`: Numeric scalar. Convergence tolerance; the Euclidean distance
between successive parameter vectors must fall below this value for the
algorithm to be considered converged. Defaults to `1e-5`.
* `max_iter`: Integer. Maximum number of alternating-optimization
iterations before the algorithm stops regardless of convergence.
Defaults to `500`.

## Details

The Ricker mean function is $\mu_t = y_{t-1} \exp(r - \alpha y_{t-1})$,
where $\alpha$ is constrained to be positive via the log
parameterization `lalpha = log(alpha)`.

## Value

A named list with the following elements:

* `estim`: Named numeric vector of point estimates for
      `r`, `alpha` (back-transformed from `lalpha`), and
      `psi` (back-transformed from `lpsi`).
* `se`: Named numeric vector of standard errors for `r`
      and `alpha`. The SE for `alpha` is obtained via the delta
      method. `psi` SE is `NA` as it is estimated via
      `optimize()` without a Hessian.
* `lower`: Numeric vector of lower 95% confidence limits for
      `r` and `alpha`. `psi` lower limit is `NA`.
* `upper`: Numeric vector of upper 95% confidence limits for
      `r` and `alpha`. `psi` upper limit is `NA`.
* `convergence`: Integer flag: `0` if the algorithm
      converged within `max_iter` iterations, `1` otherwise.

## Examples

```r
y <- readRDS("data/missingDatasets/nb_sim_randMiss_A.rds")[[1]]$y[[1]]
init <- c(r = 1, lalpha = log(0.01), lpsi = log(5))
ricker_count_nb_fit(theta = init, y = y)
```

# Constrained negative log-likelihood for a Ricker count time series

## Description

Computes the negative log-likelihood of a Ricker population model with
Poisson or Negative Binomial observations, where the density-dependence
parameter $\alpha$ is log-parameterized to enforce positivity. The
log-linear mean at time $t$ is:
$$\log \mu_t = \log N_{t-1} + r - \alpha N_{t-1}$$where $\alpha = \exp(\texttt{lalpha}) > 0$, ensuring negative density
dependence. This is the likelihood used by fitting routines that optimize
over a constrained parameter space.

## Usage

```r
ricker_count_neg_ll_cnstr(theta, y, X = NULL, fam = "poisson")
```

## Arguments

* `theta`: Named numeric vector of parameters. Must contain:
* `r`: Intrinsic growth rate (unconstrained, real-valued).
* `lalpha`: Log of the intraspecific competition coefficient.
    Exponentiated internally so that $\alpha = e^{\texttt{lalpha}} > 0$.
* `lpsi`: Log of the Negative Binomial dispersion parameter
    ($\psi = e^{\texttt{lpsi}} > 0$). Required when `fam = "neg_binom"`;
    ignored otherwise.
* `y`: Observed count data. Accepted formats:
* Numeric vector: Raw time series. Lag pairs `(yt, ytm1)` are
    constructed internally.
* Data frame or matrix: Must have named columns `yt` (count at
    time $t$) and `ytm1` (count at time $t-1$), i.e., the
    offset/lag-pair format used elsewhere in this package. Rows may come
    from multiple concatenated series when using `off_patch` mode.`NA`s should be filled in before calling this function.
* `X`: Ignored. Retained for interface compatibility with
[`ricker_count_neg_ll`](#ricker_count_neg_ll). Lag pairs are always derived from `y`.
* `fam`: Error distribution family. One of `"poisson"` (default) or
`"neg_binom"`.

## Details

Unlike [`ricker_count_neg_ll`](#ricker_count_neg_ll), this function accepts `y`in several formats and constructs the lag pairs internally, so no separate
model matrix needs to be supplied.

## Seealso

[`ricker_count_neg_ll`](#ricker_count_neg_ll) for the unconstrained version,
[`ricker_step`](#ricker_step) for the deterministic Ricker step,
[`fit_ricker_cc`](#fit_ricker_cc) and [`fit_ricker_EM`](#fit_ricker_em) which use
  this function internally.

## Value

Scalar negative log-likelihood evaluated at `theta`.

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

# Fit Ricker count model with Poisson errors

## Description

Estimates Ricker population model parameters under a Poisson observation
model via a single call to `optim`. The intra-specific competition
coefficient $\alpha$ is constrained to be positive through the log
parameterization `lalpha = log(alpha)`.

## Usage

```r
ricker_count_pois_fit(theta_init, y)
```

## Arguments

* `theta_init`: Named numeric vector of initial parameter values. Must
contain elements named `"r"` (intrinsic growth rate) and
`"lalpha"` (log of the intra-specific competition coefficient).
* `y`: Either a numeric vector of population counts through time, or a
data frame with columns `yt` (count at time $t$) and
`ytm1` (count at time $t-1$). If a vector is supplied it is
converted internally to the lagged data frame format.

## Details

The Ricker mean function is $\mu_t = y_{t-1} \exp(r - \alpha y_{t-1})$,
where counts at time $t$ are assumed to follow a Poisson distribution
with mean $\mu_t$.

## Value

A named list with the following elements:

* `estim`: Named numeric vector of point estimates for
      `r` and `alpha` (back-transformed from `lalpha`).
* `se`: Named numeric vector of standard errors. The SE for
      `alpha` is obtained via the delta method.
* `lower`: Numeric vector of lower 95% confidence limits for
      `r` and `alpha`.
* `upper`: Numeric vector of upper 95% confidence limits for
      `r` and `alpha`.
* `working`: Named numeric vector of estimates on the working
      (log) scale, i.e. `r` and `lalpha`.
* `V`: Variance-covariance matrix of the working-scale
      estimates, derived from the inverse Hessian returned by
      `optim`.
* `convergence`: Integer convergence flag from `optim`:
      `0` indicates successful convergence.

## Examples

```r
y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
init <- c(r = 1, lalpha = log(0.01))
ricker_count_pois_fit(theta_init = init, y = y)
```

# Fitting the Ricker population model to count data

## Description

This function uses an Expectation Maximization (EM) approach to fit a Ricker population model
with Poisson or Negative-Binomial demographic stochasticity.

## Usage

```r
ricker_EM(y, init_theta, fam = "poisson", tol = 1e-05, max_iter = 50)
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
  init_width = 0.1,
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

# Draw samples from a zero-truncated Negative Binomial distribution

## Description

Generates random draws from a Negative Binomial distribution conditioned on
the outcome being strictly positive ($X \geq 1$). Uses the same
inverse-CDF method as [`zt_poisson_rng`](#zt_poisson_rng): a uniform variate is
drawn on $[p_0, 1]$, where $p_0 = P(X = 0)$ under the untruncated
Negative Binomial, and then mapped through the NB quantile function.
Used during data augmentation to impute missing counts under Negative
Binomial error while avoiding zeros in the likelihood.

## Usage

```r
zt_neg_binom_rng(n, mu, size)
```

## Arguments

* `n`: Number of draws to return.
* `mu`: Mean of the underlying (untruncated) Negative Binomial distribution.
Must be positive. Corresponds to the `mu` parameterization in
[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html).
* `size`: Dispersion (size) parameter of the Negative Binomial distribution
($\psi$ elsewhere in this package). Must be positive. Larger values
approach the Poisson limit.

## Seealso

[`zt_poisson_rng`](#zt_poisson_rng) for the Poisson analogue,
[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html) for the parameterization used.

## Value

Integer vector of length `n` with all elements $\geq 1$.

# Draw samples from a zero-truncated Poisson distribution

## Description

Generates random draws from a Poisson distribution conditioned on the outcome
being strictly positive ($X \geq 1$). Uses the inverse-CDF method: a
uniform variate is drawn on $[p_0, 1]$, where $p_0 = P(X = 0)$ under
the untruncated Poisson, and then mapped through the Poisson quantile function.
This guarantees all returned values are $\geq 1$ without rejection sampling.
Used during data augmentation to impute missing counts while avoiding zeros
that would cause $\log(0)$ in the likelihood.

## Usage

```r
zt_poisson_rng(n, lambda)
```

## Arguments

* `n`: Number of draws to return.
* `lambda`: Rate parameter of the underlying (untruncated) Poisson distribution.
Must be positive.

## Seealso

[`zt_neg_binom_rng`](#zt_neg_binom_rng) for the Negative Binomial analogue.

## Value

Integer vector of length `n` with all elements $\geq 1$.

