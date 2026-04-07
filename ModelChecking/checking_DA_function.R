
library(tidyverse)

dat <- readRDS(here::here("data/missingDatasets/pois_sim_randMiss_B.rds"))
funlist <- list.files(
  path = here::here("Functions/"),
  pattern = "ricker",
  full.names = T
)

lapply(funlist, source)

# ---- Checks on Bursaria data ----
bursaria <- read_csv(here::here("data/Bursaria.csv")) %>%
  janitor::clean_names()

# make list of series
burs_list <- map(
  unique(bursaria$patch),
  ~ {select(bursaria, c(patch, date, number)) %>%
      filter(patch == .x)}
)

# first, just convert this into an offset dataframe without missingness
burs_dat_os <- map(
  burs_list,
  ~ data.frame(
    yt = .x$number[2:nrow(.x)],
    ytm1 = .x$number[1:(nrow(.x) - 1)]
  )
) %>% Reduce(rbind, .)












# ---- Simulated data ----

# get vector of actual proportion missing
prop_miss <- map(
  dat,
  ~ str_extract(
    names(.x$y),
    pattern = "0.\\d+"
  )
) %>% unlist() %>% as.numeric()
# flatten the list
dat_flat <- map(
  dat,
  ~ pluck(.x, "y")
) %>% list_flatten()
dat_flat[20007]

y <- dat_flat[[20007]]

fit_ricker_cc(y)
fit_ricker_EM(y)
fit_ricker_DA(y, samples = 1000, burnin = 6000)
