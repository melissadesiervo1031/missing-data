################################################################
# This script double checks the Ricker data augmentation function with a small simulated dataset
################################################################


# Test
set.seed(5526)
library(tidyverse)

flist <- list.files("Functions/", full.names = T)
lapply(flist, source)

y <- purrr::map(
  seq(50, 500, by = 50),
  ~ ricker_sim(.x, 0.8, 0.02, N0 = 40)
)

nas <- purrr::map(
  y,
  ~ {
    sample(1:length(.x), size = round(length(.x) * 0.2))
  }
)

y_miss <- map2(
  y, nas,
  ~ replace(.x, .y, NA)
)

fit_mcmc <- purrr::map(
  y,
  ~ fit_ricker_DA(
    .x,
    burnin = 5000,
    priors_list = list(
      m_r = 0,
      sd_r = 2.5,
      m_lalpha = -3,
      sd_lalpha = 1
    ),
    nthin = 5,
    return_y = T
  )
)

fit_cc <- purrr::map(
  y,
  ~ fit_ricker_cc(.x)
)


df <- tibble(
  n = rep(seq(50, 500, by = 50), each = 2) %>% rep(., 2),
  par = rep(c("r", "alpha"), 20),
  method = c(rep("DA", 20), rep("cc", 20)),
  estim = c(
    unlist(map(fit_mcmc, ~ .x$estim)),
    unlist(map(fit_cc, ~ .x$estim))
  ),
  low = c(
    unlist(map(fit_mcmc, ~ .x$lower)),
    unlist(map(fit_cc, ~ .x$lower))
  ),
  high = c(
    unlist(map(fit_mcmc, ~ .x$upper)),
    unlist(map(fit_cc, ~ .x$upper))
  )
)

ggplot(data = df, aes(x = n, y = estim, color = method)) +
  facet_wrap(vars(par), scales = "free_y") +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
  geom_point() +
  theme_classic() +
  ylab("Estimate")

ggsave(
  here("figures/DA_example.png"),
  device = "png",
  height = 3,
  width = 6,
  units = "in"
)
  