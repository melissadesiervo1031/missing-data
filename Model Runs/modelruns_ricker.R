
# This script isn't done yet!!!

library(here)
library(parallel)
library(tidyverse)


in_args <- commandArgs(trailingOnly = T)

dat <- readRDS(here(in_args[1]))

pars <- readRDS(here(in_args[2]))

n_autocorrs <- length(dat) / nrow(pars)

nmiss_props <- length(dat[[1]]$y)

pars_full <- pars[rep(1:nrow(pars), n_autocorrs), ]
pars_full <- pars_full[rep(1:nrow(pars), each = nmiss_props), ]

# flatten the list
dat_flat <- list_flatten(dat) %>% list_flatten()

# add missingness and autocorrelation specs to pars df
pars_full <- pars_full %>% mutate(
  autoCorr = as.numeric(str_extract(
    names(dat_flat),
    pattern = "autoCorr_\\d+"
  ) %>% str_extract(
    ., pattern = "\\d+"
  )),
  propMiss = rep(
    seq(0.1, 0.7, by = 0.05)
  )
)


cl <- parallel::makeCluster(as.numeric(in_args[4]))

em_fits <- parLapply(
  cl = cl,
  fun = 
)
