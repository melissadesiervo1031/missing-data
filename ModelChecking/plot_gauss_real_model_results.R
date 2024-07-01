library(tidyverse)

pr <- read_csv('data/pine_river_data_prepped.csv')
fit <- readRDS('data/model_results/gauss_real_ModelResults.rds')

forecast_mar <- read_csv('data/model_results/gauss_real_MAR_brms_FORECASTpreds.csv')
forecast_mnar <- read_csv('data/model_results/gauss_real_MNAR_brms_FORECASTpreds.csv') %>%
  mutate(missingness = 'MNAR')

pars <- read_csv('data/model_results/gauss_real_MAR_brmsFORECASTS.csv') %>%
  filter(missingprop_autocor == 'y') %>%
  group_by(parameter) %>%
  summarize(mean = mean(mean))

pr$GPP_t1 = c(pr$GPP[1], pr$GPP[1:365]) 

pr$preds = pr$GPP_t1 * pars$mean[4] + pr$light.rel * pars$mean[3] + pr$Q * pars$mean[2] + pars$mean[1]

pr %>%
  select(date, GPP, preds) %>%
  pivot_longer(cols = -date, names_to = ) 
ggplot(pr, aes(date, GPP)) + geom_line() + geom_line(aes(y = preds), col = 'brown3') +
  geom_vline(xintercept = as.Date('2016-11-30'), col = 'black')+
  theme_classic()
forecast <- bind_rows(forecast_mar, forecast_mnar) %>%
  filter(missingprop_autocor != 'y') %>%
  group_by(missingness, date) %>%
  summarize(Est = mean(Estimate, na.rm = T),
            SD = sd(Estimate, na.rm = T))
  # mutate(missingprop_autocor = case_when(missingprop_autocor == 'y' ~ propMissAct_0.0_autoCorr_0.0,
  #                                        TRUE ~ missingprop_autocor),
  #        percent_missing = grep(''))

pred <- bind_rows(forecast_mar, forecast_mnar) %>%
  filter(missingprop_autocor == 'y') %>%
  group_by(date) %>%
  summarize(Est = mean(Estimate, na.rm = T))
  

ggplot(forecast, aes(date, Est, col = missingness))+
  geom_line(data = filter(pr, date >= as.Date('2016-11-01')),
            aes(date, GPP), col = 'black') +
  geom_point() + geom_line() +
  geom_ribbon(aes(ymax = Est + SD, ymin = Est - SD, fill = missingness),
              alpha = 0.2, col = NA) +
  geom_vline(xintercept = as.Date('2016-11-30'), col = 'black')+
  geom_line(data = pred, lty = 2, col = 'black', size = 0.8) +
  theme_classic() +
  ylab('GPP') + xlab('Date')


head(fit)
