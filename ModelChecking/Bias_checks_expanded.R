################################################################
# This script double checks whether the bias in the EM estimator 
# diminishes with sample size (i.e., asymptotically unbiased)
################################################################

# load libraries
library(tidyverse)
library(here)
library(parallel)
source(here::here("Functions/missing_data_functions.R"))
ricker_funs <- list.files(here::here("Functions"), pattern = "ricker", full.names = T)
lapply(ricker_funs, source)


# establish some global params
set.seed(2533)
n_params <- 250
uns <- c(50, 100, 250, 500)
nns <- length(uns)

# create data
ns <- rep(uns, each = n_params)
r <- runif(nns * n_params, min = 0.4, max = log(2))
alpha <- runif(nns * n_params, min = 0.01, max = 0.04)  

y <- mapply(
  ricker_sim,
  n = ns, r = r, alpha = alpha, N0 = 40
)

# add these to a dataframe
dat <- tibble(
  y = y,
  r = r,
  alpha = alpha,
  n = ns
)

# find and remove simulations in which pop went extinct
dat <- dat %>% mutate(
  ext = map_dbl(
    dat$y,
    ~ as.numeric(.x[length(.x)] == 0)
  )
)
dat <- dat %>% filter(
  ext == 0
)

# create datasets with 20% missing
dat <- dat %>% mutate(
  y_miss = map(
    y,
    ~ makeMissing(.x, "random", propMiss = 0.2)
  ) %>% flatten()
)

# create vector of methods
mthds <- c("drop", "cc", "EM", "MI", "DA")

fit_all <- function(y, mthds){
  df <- data.frame(
    method = rep(mthds, each = 2),
    param = rep(c("r", "alpha"), length(mthds)),
    estims = 0
  )
  drop <- fit_ricker_drop(y)
  cc <- fit_ricker_cc(y)
  EM <- fit_ricker_EM(y)
  MI <- fit_ricker_MI(y)
  DA <- fit_ricker_DA(y)
  if(any(DA$rhat > 1.2)){
    da_estim <- c(r = NA, alpha = NA)
  } else{
    da_estim <- DA$estim
  }
  df$estims <- c(
    drop$estim,
    cc$estim,
    EM$estim,
    MI$estim,
    da_estim
  )
  return(df)
}

# fit the models using each method
# cl <- makeCluster(10)
# clusterEvalQ(
#   cl,
#   {source(here::here("Functions/missing_data_functions.R"))
#     ricker_funs <- list.files(here::here("Functions"), pattern = "ricker", full.names = T)
#     lapply(ricker_funs, source)}
# )
# 
# fits <- parLapply(
#   cl,
#   X = dat$y_miss,
#   fun = fit_all,
#   mthds = mthds
# )
# stopCluster(cl)

fits <- lapply(
  dat$y_miss,
  fit_all,
  mthds = mthds
)

dat <- dat %>% mutate(
  res = fits
)

# pivot each sub-df wider to make the next
# calculation easier
dat <- dat %>% mutate(
  res_wide = map(
    res,
    ~ pivot_wider(.x, names_from = param, values_from = estims)
  )
)

# disambiguate names
dat <- dat %>% mutate(
  res_wide = map(
    res_wide,
    ~ rename(
      .x,
      r_hat = r,
      alpha_hat = alpha
    )
  )
)

# expand this df to make a results df in long format
results <- dat %>%
  select(c(r:n, res_wide)) %>%
  unnest(cols = c(res_wide))

# Calculate the relative bias for each param
results <- results %>% mutate(
  bias_r_hat = (r_hat - r) / r,
  bias_alpha_hat = (alpha_hat - alpha) / alpha
)

# summarize
results_sum <- results %>% group_by(n, method) %>%
  summarise(
    mbias_r_hat = mean(bias_r_hat, na.rm = T),
    mbias_alpha_hat = mean(bias_alpha_hat, na.rm = T),
    se_r_hat = sd(bias_r_hat, na.rm = T),
    se_alpha_hat = sd(bias_alpha_hat, na.rm = T),
    .groups = "drop"
  )

# now pivot longer for plotting
results_long <- results_sum %>%
  select(n:mbias_alpha_hat) %>%
  pivot_longer(
    cols = c(mbias_r_hat, mbias_alpha_hat),
    names_to = "param",
    values_to = "rel_bias"
  )

results_long <- cbind(
  results_long,
  results_sum %>%
    select(c(n, method, se_r_hat, se_alpha_hat)) %>%
    pivot_longer(
      cols = c(se_r_hat, se_alpha_hat),
      names_to = "param",
      values_to = "se"
    ) %>%
    select(se)
)

# save the dataset in case we want to make changes to the figure
saveRDS(dat, file = here::here("data/bias_checks.rds"))

# create panel labels
results_long <- results_long %>% mutate(
  param = factor(
    param,
    levels = c("mbias_r_hat", "mbias_alpha_hat"),
    labels = c(expression(hat(r)), expression(hat(alpha)))
  )
)

# plot the results
ggplot(results_long, aes(x = n, y = rel_bias, color = method)) +
  facet_wrap(~param, labeller = label_parsed) +
  geom_errorbar(
    aes(ymin = rel_bias - se, ymax = rel_bias + se),
    width = 0.1
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  ylab("Relative bias") +
  xlab("Sample size (20% missing)") +
  scale_color_brewer(palette = "Dark2")

# save the plot
ggsave(
  file = here::here("figures/bias_checks_ricker.png"),
  width = 6, height = 3,
  units = "in", dpi = 300
)  




