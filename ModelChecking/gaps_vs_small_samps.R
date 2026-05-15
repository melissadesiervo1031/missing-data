# This script explores if there are differences between 
# having gaps in the time series compared to just having
# a shorter time series.

# ---- Libraries ----

library(tidyverse)
fns <- list.files(here::here("Functions"), full.names = T, pattern = ".R$")
lapply(fns, source)

# ---- Create data ----

set.seed(260515)
reps <- 1000

## ---- real-valued data ----

phi <- 0.6
sd_err <- 0.5

p_miss <- c(0.1, 0.2, 0.4)

ns <- round(365 * (1 - p_miss))

# shorter time series, but no intermediate missingness
ys_short <- tibble(
  n = rep(ns, each = reps)
) %>%
  mutate(
    y = purrr::map(
      n,
      ~ arima.sim(model = list(ar = phi), sd = sd_err, n = .x)
    )
  )

ys_miss <- tibble(
  prop_miss = rep(p_miss, each = reps),
  y_full = lapply(1:nrow(ys_short), function(x){
    arima.sim(model = list(ar = phi), sd = sd_err, n = 365) |>
      as.double()
  })
) %>%
  mutate(
    y_miss = purrr::map2(
      y_full,
      prop_miss,
      ~ makeMissing(.x, "random", propMiss = .y, autoCorr = 0.3, type = "Gaussian")
    )
  )

# save to data dir
saveRDS(ys_short, file = here::here("data/short_ts_compare_to_miss.rds"))
saveRDS(ys_miss, file = here::here("data/miss_ts_compare_to_short.rds"))



## ---- Count data ----

r <- log(1.6)
alpha <- 0.01

ns_c <- round(60 * (1 - p_miss))

ys_short_c <- tibble(
  n = rep(ns_c, each = reps),
  y = purrr::map(
    n,
    ~ ricker_sim(.x, r, alpha, N0 = sample(10:20, size = 1), err_fam = "poisson")
  )
)

ys_miss_c <- tibble(
  
)






