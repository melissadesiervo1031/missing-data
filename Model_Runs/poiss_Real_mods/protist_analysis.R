
# ---- Libraries and data ----

library(tidyverse)
library(DHARMa)
library(lme4)
func_list <- list.files(here::here("Functions"), pattern = "\\.R$")
lapply(func_list, function(f){ source(here::here("Functions", f)) })

# other, one-off packages that are needed:
## - here
## - janitor


# load data
bursaria <- read_csv(here::here("data/Bursaria.csv"))


# ---- Data formatting ----

bursaria <- bursaria %>% 
  janitor::clean_names() %>%
  mutate(
    date = lubridate::mdy(date),
    patch = factor(patch)
  )

# grouping by patch, then adding a number of days column
bursaria_patchlist <- purrr::map(
  unique(bursaria$patch),
  ~ {
    filter(bursaria, patch == .x) %>%
      arrange(date) %>%
      mutate(
        days = as.numeric(date - date[1], units = "days")
      )
  }
)

# now offset each so that we have a response that is time t and explanatory variable that is
# time t - 1
bursaria_patchlist_diff <- purrr::map(
  bursaria_patchlist,
  ~ {
    df2 <- tibble(
      count_tm1 = .x$number[-nrow(.x)]
    )
    cbind(.x[-1, ], df2)
  }
)

# recombine into one dataframe
bursaria_diff <- Reduce(rbind, bursaria_patchlist_diff)
bursaria2 <- Reduce(rbind, bursaria_patchlist)


# ---- Analysis with glm ----

## ---- Poisson errors ----
fit_init <- glm(
  number ~ count_tm1,
  offset = log(count_tm1),
  family = poisson,
  data = bursaria_diff
)

resids <- simulateResiduals(fit_init)
plot(resids)

testDispersion(resids)

# data are clearly over-dispersed for Poisson model. 
# option 1. Negative binomial

## ---- Negative binomial ----
fit_nb <- MASS::glm.nb(
  number ~ count_tm1 + offset(log(count_tm1)),
  data = bursaria_diff
)

resids_nb <- simulateResiduals(fit_nb)
plot(resids_nb)

plotResiduals(resids_nb, form = bursaria_diff$count_tm1)

# NOTES: Fits better. Saturation maybe not perfectly log-linear, but
# doesn't look too bad

summary(fit_nb)


## ---- Poisson model with random effects for each rep ----
fit_pois_re <- glmer(
  number ~ count_tm1 + (1 | patch),
  offset = log(count_tm1),
  family = poisson(),
  data = bursaria_diff
)

AIC(fit_pois_re, fit_nb)

# NOTES: Negative binomial model is simpler and is a reasonable model for the data

## ---- Plotting model predictions over the data ----

estims <- coefficients(fit_nb)

# get geometric mean abundance for each time step
mpreds_df <- bursaria_diff %>% 
  group_by(days) %>%
  summarize(
    mean_tm1 = exp(mean(log(count_tm1))),
    .groups = "drop"
  )

X_pred <- cbind(1, mpreds_df$mean_tm1)
V <- vcov(fit_nb)

mpreds_df <- mpreds_df %>%
  mutate(
    estim = exp(X_pred %*% estims + log(X_pred[,2])),
    se = sqrt(diag(X_pred %*% V %*% t(X_pred))),
    lower = exp(log(estim) - 2 * se),
    upper = exp(log(estim) + 2 * se),
    low_pred = qnbinom(0.025, size = fit_nb$theta, mu = estim),
    high_pred = qnbinom(0.975, size = fit_nb$theta, mu = estim)
  )

ggplot(bursaria2, aes(x = days, y = number)) +
  geom_line(aes(group = patch), alpha = 0.9, linewidth = 0.2) +
  geom_ribbon(
    data = mpreds_df,
    aes(x = days, y = estim, ymin = low_pred, ymax = high_pred),
    fill = "blue", alpha = 0.5
  ) +
  geom_line(
    data = mpreds_df,
    aes(x = days, y = estim),
    linewidth = 2
  ) +
  theme_classic()
  
