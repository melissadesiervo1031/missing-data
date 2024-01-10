# Run stan model with fixed observation error on actual GPP data
# Author: AM Carter

library(tidyverse)
library(knitr)
library(rstan)
library(brms)
options(mc.cores = parallel::detectCores())

dat <- read_csv('data/NWIS_MissingTS_subset_new.csv')
# mdat <- read_csv('data/NWIS_MissingTSinfo_subset.csv')

# functions
plot_gpp <- function(df, forecast = NA, pp_fit = FALSE){
    ratio_QL <- max(df$light)/max(df$Q)
    GPP_plot <- ggplot(df, aes(date, GPP))+
      #geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="chartreuse4")+
      geom_point(size = 1, color="chartreuse4") + 
      geom_line(size = 0.8, color="chartreuse4")+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$long_name[1])+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12))
    if(pp_fit){
        GPP_plot <- df %>%
          left_join(dplyr::select(forecast, date, post_est = Estimate,
                           post_err = Est.Error), by = 'date')%>%
          mutate(GPP.lower = post_est - post_err,
                 GPP.upper = post_est + post_err) %>%
          pivot_longer(cols = c('GPP', 'post_est'), 
                       names_to = 'GPP', values_to = 'value') %>%
          ggplot(aes(date, value, col = GPP))+
          geom_errorbar(aes(ymin = GPP.lower, 
                            ymax = GPP.upper), 
                        width=0.2,color="grey")+
          geom_point(size=1) + geom_line(size=0.8)+
          scale_color_manual(values = c('chartreuse4', 'grey35'))+
          labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$long_name[1])+
          theme(legend.position = 'top', legend.title = element_blank(),
                panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
                axis.text.y = element_text(size=12),
                axis.title.y = element_text(size=12))
    }
    data_plot <- ggplot(df, aes(date, Q*ratio_QL))+
      geom_point(data=df, aes(date, light), size=1.5, color="darkgoldenrod3")+
      geom_line(size=1, color="deepskyblue4")+
      scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
      labs(y=expression('Daily PPFD'))+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text = element_text(size=12),
            axis.title.y.left = element_text(size=12, color="darkgoldenrod3"),
            axis.title.y.right = element_text(size=12, color="deepskyblue4"),
            axis.text.x = element_text(angle=25, hjust = 1),
            strip.background = element_rect(fill="white", color="black"),
            strip.text = element_text(size=15))
    
    pp <- ggpubr::ggarrange(GPP_plot, data_plot, ncol = 1, align = 'v', 
                            heights = c(1, 0.5))
    
    return(pp)
    
}

fit_brms_model <- function(df, iter = 4000,
                           forecast = TRUE, forecast_days = 365){

  if(forecast){
    df_sub <- df[1:(nrow(df)-forecast_days), ]  # Remove to save these for forecasting
  }
  
  # Make the model formula and priors
  bform <- brms::bf(GPP | mi() ~ light + Q + ar(p = 1))
  bprior <- c(prior(normal(0,1), class = 'ar'),
              prior(normal(0,5), class = 'b'))
  
  # fit model to list of datasets
  bmod <- brms::brm(bform, data = df_sub, 
                    prior = bprior, iter = iter)
  
  extract_brms_pars <- function(bfit){
    bsum <- brms::posterior_summary(bfit, probs = c(0.025, 0.5, 0.975))
    bsum <- as.data.frame(bsum) %>%
      mutate(parameter = row.names(bsum)) %>%
      filter(!(parameter %in% c('lprior', 'lp__'))) %>%
      mutate(parameter = case_when(parameter == 'ar[1]' ~ 'phi',
                                   TRUE ~ parameter)) %>%
      select(parameter, mean = Estimate, sd = Est.Error, 
             '2.5%' = Q2.5, '50%' = Q50, '97.5%' = Q97.5)
    
    row.names(bsum) <- NULL
    bsum <- bsum[grep('^Ymi', bsum$parameter, invert = TRUE),]
    
    return(bsum)
  }
  
  bpars <- extract_brms_pars(bmod)
  
  if(forecast){  
    # dat_forecast <- df %>%
    #   slice((nrow(df)-forecast_days):nrow(df)) %>%
    #   select(date, GPP, light, Q)
    
    predictions <- predict(bmod, newdata = df[,-2]) %>%
        as.data.frame() %>% mutate(date = df$date,
                                   GPP = df$GPP)
    
    return(list(brms_fit = bmod, 
                brms_forecast = predictions,
                brms_pars = bpars))
  }
  
  return(list(brms_fit = bmod,
              brms_forecast = predictions,
              brms_pars = bpars))
  
}

# site list:
sites <- dat %>% select(site_name, long_name) %>%
  unique()

# Fit site and plot 
fit_site <- function(site_id){
  ss <- filter(dat, site_name == site_id) %>%
    select(-Q) %>%
    arrange(date)
  dates <- data.frame(date = seq(ss$date[1], ss$date[nrow(ss)], by = 'day'))
  ss <- left_join(dates, ss, by = 'date') %>%
    mutate(site_name = ss$site_name[1],
           long_name = ss$long_name[1], 
           Year = year(date),
           DOY = format(date, '%j'))
  
  q <- dataRetrieval::readNWISdata(sites = substr(site_id, 6, nchar(site_id)),
                                   parameterCd = "00060",
                                   service = "dv",
                                   startDate = ss$date[1], 
                                   endDate = ss$date[nrow(ss)]) %>%
    select(date = dateTime,
           discharge_cfs = X_00060_00003) %>%
    mutate(Q = discharge_cfs / (3.28084)^3 ) # convert discharge to cms
  
  ss <- left_join(ss, select(q, date, Q), by = 'date') %>%
    mutate(Q = log(Q),
           light = zoo::na.approx(light))
  
  fit <- fit_brms_model(ss, forecast_days = 365)
  
  # p <- plot_gpp(ss, forecast = fit$brms_forecast, pp_fit = TRUE)
  # print(p)
  
  return(list(fit = fit,
              dat = ss))
  
}


# individual site fits:

pdf('figures/model_fits_to_real_gauss_datasets.pdf')
mod_list <- vector("list", nrow(sites))


for(i in 1:nrow(sites)){
 
  fit <- fit_site(sites$site_name[i])
  
  p <- plot_gpp(fit$dat, forecast = fit$fit$brms_forecast, pp_fit = TRUE)
  print(p)

  
  file_path <- "figures/model_fits_real_gauss_summary.txt"
  con <- file(file_path, "a")  # "a" stands for append mode
  
  # Write (append) the text to the file
  model_summary <- capture.output(summary(fit$fit$brms_fit))  # Capture the summary output
  writeLines(sites$long_name[i], con)
  writeLines(model_summary, con)
  
  # Close the file connection
  close(con)
  mod_list[[i]] <- fit
  

}
dev.off()


# select and prepare dataset for using in manuscript
sites
# We're using 2012-2014 in Au Sable River (2) and  2014-2016 Badger Mill Creek (4)
site_id = sites$site_name[2]
ss <- filter(dat, site_name == site_id) %>%
  select(-Q) %>%
  arrange(date)
dates <- data.frame(date = seq(as.Date('2012-01-01'), as.Date('2014-12-31'), 
                               by = 'day'))
ss <- left_join(dates, ss, by = 'date') %>%
  mutate(site_name = ss$site_name[1],
         long_name = ss$long_name[1], 
         Year = year(date),
         DOY = format(date, '%j'))

q <- dataRetrieval::readNWISdata(sites = substr(site_id, 6, nchar(site_id)),
                                 parameterCd = "00060",
                                 service = "dv",
                                 startDate = ss$date[1], 
                                 endDate = ss$date[nrow(ss)]) %>%
  select(date = dateTime,
         discharge_cfs = X_00060_00003) %>%
  mutate(Q = discharge_cfs / (3.28084)^3 ) # convert discharge to cms

ss <- left_join(ss, select(q, date, Q), by = 'date') %>%
  mutate(Q = log(Q),
         GPP = log(GPP),
         light = zoo::na.approx(light))

write_csv(ss, 'data/au_sable_river_prepped.csv')


site_id = sites$site_name[4]
ss <- filter(dat, site_name == site_id) %>%
  select(-Q) %>%
  arrange(date)
dates <- data.frame(date = seq(as.Date('2013-01-01'), as.Date('2015-12-30'), 
                               by = 'day'))
ss <- left_join(dates, ss, by = 'date') %>%
  mutate(site_name = ss$site_name[1],
         long_name = ss$long_name[1], 
         Year = year(date),
         DOY = format(date, '%j'))

q <- dataRetrieval::readNWISdata(sites = substr(site_id, 6, nchar(site_id)),
                                 parameterCd = "00060",
                                 service = "dv",
                                 startDate = ss$date[1], 
                                 endDate = ss$date[nrow(ss)]) %>%
  select(date = dateTime,
         discharge_cfs = X_00060_00003) %>%
  mutate(Q = discharge_cfs / (3.28084)^3 ) # convert discharge to cms

ss <- left_join(ss, select(q, date, Q), by = 'date') %>%
  mutate(Q = log(Q),
         GPP = log(GPP),
         light = zoo::na.approx(light))

write_csv(ss, 'data/badger_mill_creek_prepped.csv')
