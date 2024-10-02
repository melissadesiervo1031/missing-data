################################################################################
# This script generates datasets with missingness 
# based on multiple mechanisms and then fits the Ricker model
# using EM
################################################################################

# load libraries
  library(tidyverse)
  library(here)
  library(parallel)

# load user-defined functions
  our_functions <- list.files(here("Functions/"), pattern = ".R", full.names = T)
  lapply(our_functions, source)
  
# load data
  ricker_dat <- readRDS(here("data/ricker_0miss_datasets.rds"))
  
# set seed for reproducibility
  set.seed(85612)
  
# create tibble that stores all the data
  ricker_dat_all <- tibble(
    miss_prop = rep(seq(0, 0.6, by = 0.1), 3),
    auto_cor = rep(c(0, 0.6, NA), each = 7),
    miss_type = rep(c("random", "random", "minMax"), each = 7)
  )
  
  # create list of lists columns
  fn <- function(miss_prop, auto_cor, miss_type, dat = ricker_dat){
    if(miss_prop == 0){
      purrr::map(
        dat,
        pluck("y")
      )
    } else{
      lapply(
        dat,
        FUN = function(d, p, a, type){
          if(type == "minMax"){
            makeMissing(
              d$y,
              typeMissing = type,
              propMiss = p,
              autoCorr = NULL
            )
          } else{
            makeMissing(
              d$y,
              typeMissing = type,
              propMiss = p,
              autoCorr = a
            )
          }
        },
        p = miss_prop,
        a = auto_cor,
        type = miss_type
      )
    }
  }
  
# map the function to the rows of the tibble
  ricker_dat_all <- ricker_dat_all %>% mutate(
    dat = pmap(
      ricker_dat_all,
      fn
    )
  )
  
# fit the models using EM
#  this part takes a while, but parallelization could be
#  improved on Beartooth
  cl <- makeCluster(4)

# load libraries on the cluster
  clusterEvalQ(
    cl,
    {library(here);
      our_functions <- list.files(here("Functions/"), pattern = ".R", full.names = T);
      lapply(our_functions, source)}
  )
  
# add column to data using clusterApply
  ricker_dat_all <- ricker_dat_all %>% mutate(
    estims = lapply(
      ricker_dat_all$dat,
      FUN = function(l){
        clusterApply(
          cl,
          l,
          fun = function(x){
            ricker_EM(unlist(x), init_theta = c(0.5, -0.05))$theta
          }
        )
      }
    )
  )
  stopCluster(cl)
  
  # now make relative differences column
  ricker_dat_all <- ricker_dat_all %>% mutate(
    rel_diff = lapply(
      estims,
      function(l){
        lapply(
          1:length(l),
          function(i){
            theta_true <- c(ricker_dat[[i]]$sim_params$r, -ricker_dat[[i]]$sim_params$alpha)
            (l[[i]] - theta_true) / abs(theta_true)
          }
        )
      }
    )
  )

# now summarize for r and alpha
  ricker_dat_all <- ricker_dat_all %>% mutate(
    mean_rd_r = sapply(
      ricker_dat_all$rel_diff,
      function(l){
        mean(
          map_dbl(l, ~ .x[1])
        )
      }
    ),
    mean_rd_alpha = sapply(
      ricker_dat_all$rel_diff,
      function(l){
        mean(
          map_dbl(l, ~ .x[2])
        )
      }
    ),
    sd_rd_r = sapply(
      ricker_dat_all$rel_diff,
      function(l){
        sd(
          map_dbl(l, ~ .x[1])
        )
      }
    ),
    sd_rd_alpha = sapply(
      ricker_dat_all$rel_diff,
      function(l){
        sd(
          map_dbl(l, ~ .x[2])
        )
      }
    )
  )
  
# add more informative type column
  ricker_dat_all <- ricker_dat_all %>% mutate(
    type2 = factor(rep(
      c("rho == 0", "rho == 0.6", "MNAR"),
      each = 7
    ),
    levels = c("rho == 0", "rho == 0.6", "MNAR")
  ))

# pivot the select columns to long format
  plot_df <- pivot_longer(
    select(ricker_dat_all, miss_prop, type2, mean_rd_r, mean_rd_alpha),
    cols = mean_rd_r:mean_rd_alpha,
    values_to = "mean_rd",
    names_to = "param"
  )
  plot_df <- plot_df %>% mutate(
    se = pivot_longer(
      select(ricker_dat_all, miss_prop, type2, sd_rd_r, sd_rd_alpha),
      cols = sd_rd_r:sd_rd_alpha,
      values_to = "se",
      names_to = "param2"
    )$se
  )

# replace values in param column with better labels
  plot_df <- plot_df %>% mutate(
    param = rep(c("r", "alpha"), nrow(plot_df) / 2)
  )
  
# create a facet plot
  p <- ggplot(data = plot_df, aes(x = miss_prop, y = mean_rd)) +
    geom_hline(yintercept = 0, color = "brown", linetype = "dashed") +
    facet_grid(
      vars(param), vars(type2), scales = "free",
      labeller = label_parsed
    ) +
    geom_point() +
    geom_errorbar(aes(ymin = mean_rd - se/sqrt(865), ymax = mean_rd + se/sqrt(865)), width = 0) +
    theme_bw() +
    xlab("Proportion missing") +
    ylab("(estim - true) / true")
  
# save the plot
  ggsave(
    here("Population sim and real/Figures/mean_rel_diff_ricker_EM.png"),
    p, height = 4, width = 6, units = "in"
  )
  
  
  
  
  
  
  
  