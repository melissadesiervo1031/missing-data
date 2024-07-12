# File descriptions

* `gaussian_ar1_data_sims.R`: Generates gaussian AR1 timeseries datasets with arbitrary covariate matricies and coefficients. Used to create the 1000 simulations for the remaining analyses and saves them in gauss_ar1_0miss_datasets.rds.

* `ricker_data_sims.R`: Generates integer valued time series using a ricker population model with arbitrary parameters, removing any time series where populations went extinct. Used to create the 1000 simulations for the remaining analyses, saves them in ricker_0miss_datasets.rds.

* `MakeMissingDatasets.R `: Uses the makeMissing function to generate missingness in gaussian and ricker time series. Allows for a specification of the total percent missing as well as the degree of autocorrelation in the missing time steps. This script applies makeMissing to the sets of simulated data as well as the real data examples for the productivity and population time series.