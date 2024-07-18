# Load packages
library(tidyverse)

sim_dat <- readRDS('data/gauss_ar1_0miss_datasets.rds')
DF <- data.frame()

for(i in 1:length(sim_dat)){
  
  ss <- sim_dat[[i]]

  df <- data.frame(
    parameter = c('phi', 'beta1', 'beta2', 'beta3', 'sd'),
    value = c(ss$sim_params$phi, 
              ss$sim_params$beta, 
              1)
  )

  arima_mod <- arima(ss$y, order = c(1,0,0), xreg = ss$sim_params$X[,2:3])
  df <-df %>% 
      mutate(estimate = c(arima_mod$coef,  arima_mod$sigma2), 
             estimate_se = c(sqrt(diag(arima_mod$var.coef)), NA),
             rse = abs(value - estimate),
             model = i)
  
  DF = bind_rows(DF, df)

}


filter(DF, parameter != 'sd') %>%
  ggplot(aes(value, estimate)) +
  geom_point() + facet_wrap(.~parameter, scales = 'free') + 
  theme_bw()
 
DF %>%
  mutate(coverage = case_when(estimate_se > rse ~ 1,
                              TRUE ~ 0)) %>%
  group_by(parameter) %>%
  summarize(rmse = sqrt(mean(rse^2)),
            coverage = mean(coverage))

# note I couldn't get percent coverage for the estimate of sd because the
# arima model doesn't give a standard error on that estimate for that I can find.


# Also, we may want to convert the se into a 95% confidence interval before calculating coverage