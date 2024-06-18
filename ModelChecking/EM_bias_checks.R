################################################################
# This script double checks whether the bias in the EM estimator 
# diminishes with sample size (i.e., asymptotically unbiased)
################################################################

# load libraries
  library(tidyverse)
  library(here)
  library(patchwork)
  func_list <- list.files(here("Functions/"), pattern = ".R", full.names = T)
  lapply(func_list, source)

# establish some global params
  set.seed(9621)
  n_params <- 100
  uns <- c(50, 100, 1000, 5000, 10000)
  nns <- length(uns)
  
# create data
  ns <- rep(uns, each = n_params)
  r <- rep(runif(n_params, min = 0.5, max = log(2)), nns) 
  alpha <- rep(runif(n_params, min = 0.01, max = 0.05), nns)  

  y <- mapply(
    ricker_sim,
    n = ns, r = r, alpha = alpha, N0 = 10
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
  
# create datasets with 15% missing
  dat <- dat %>% mutate(
    y_miss = map(
      y,
      ~ makeMissing(.x, "random", propMiss = 0.2)
    ) %>% flatten()
  )
  
# fit the models
  dat <- dat %>% mutate(
    estims_full = lapply(
      y,
      FUN = ricker_EM,
      init_theta = c(0.6, -0.025)
    ),
    estims_miss = lapply(
      y_miss,
      FUN = ricker_EM,
      init_theta = c(0.6, -0.025)
    )
  )

# pull out estimates
  dat <- dat %>% mutate(
    rhat_full = map_dbl(
      estims_full,
      ~ .x$theta[1]
    ),
    ahat_full = map_dbl(
      estims_full,
      ~ .x$theta[2]
    ),
    rhat_miss = map_dbl(
      estims_miss,
      ~ .x$theta[1]
    ),
    ahat_miss = map_dbl(
      estims_miss,
      ~ .x$theta[2]
    )
  )

# results
  results <- select(
    dat,
    c(r, alpha, n, rhat_full, rhat_miss)
  ) %>% pivot_longer(
    cols = c(rhat_full, rhat_miss),
    names_to = "scenario",
    values_to = "rhat"
  ) %>% cbind(
    .,
    pivot_longer(
      data = select(
        dat,
        c(r, alpha, n, ahat_full, ahat_miss)
      ),
      cols = c(ahat_full, ahat_miss),
      names_to = "scenario",
      values_to = "ahat"
    ) %>% select(., ahat)
  ) %>% mutate(
    scenario = replace(
      scenario,
      scenario == "rhat_full",
      "full"
    ),
    scenario = replace(
      scenario,
      scenario == "rhat_miss",
      "miss"
    )
  )

# create a column with relative bias
  results <- results %>% mutate(
    rel_bias_r = (rhat - r) / r,
    rel_bias_alpha = (ahat - (-alpha)) / abs(alpha)
  )

# summarize
  results_sum <- group_by(results, n, scenario) %>%
    summarise(
      mean_bias_r = mean(rel_bias_r),
      sd_r = sd(rel_bias_r),
      mean_bias_alpha = mean(rel_bias_alpha),
      sd_a = sd(rel_bias_alpha)
    ) %>% ungroup() %>% mutate(
      n_adj = log(n) + log(0.2) * (1 - 1:nrow(.) %% 2)
    )
    

# create plots
  r <- ggplot(data = results_sum, aes(x = n_adj, y = mean_bias_r, color = scenario)) +
    geom_errorbar(
      aes(ymin = mean_bias_r - sd_r, ymax = mean_bias_r + sd_r),
      width = 0.1
    ) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_point() +
    geom_line() +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_manual(
      values = c("blue", "red"),
      labels = c("Full data", "20 percent missing")
    ) +
    xlab("ln(Number of observations)") +
    ylab("(estimate - true)/true") +
    ggtitle(expression(r))
  
  alpha <- ggplot(data = results_sum, aes(x = n_adj, y = mean_bias_alpha, color = scenario)) +
    geom_errorbar(
      aes(ymin = mean_bias_alpha - sd_a, ymax = mean_bias_alpha + sd_a),
      width = 0.1
    ) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_point() +
    geom_line() +
    theme_classic() +
    scale_color_manual(
      values = c("blue", "red"),
      labels = c("Full data", "20 percent missing")
    ) +
    xlab("ln(Number of observations)") +
    ylab("") +
    ggtitle(expression(alpha))
    

  ggsave(
    filename = here("Population sim and real/Figures/bias_checks_EM.png"),
    plot = r + alpha,
    device = "png",
    width = 8, height = 3,
    units = "in", dpi = 300
  )
  
  
  
  
  