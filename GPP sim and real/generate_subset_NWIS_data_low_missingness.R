library(tidyverse)

dat <- readRDS('../except_heterotrophy/data_ignored/high_quality_daily_metabolism_with_SP_covariates_complete.rds')

# look for a site that has at least 90% of the data for 3 years
sites <- dat %>% group_by(site_name, Year) %>%
  summarize(n = sum(!is.na(GPP))) %>%
  filter(n > 330) %>%
  select(site_name) 

dat %>% filter(site_name %in% sites$site) %>%
  group_by(site_name, Year) %>%
  summarize(n = sum(!is.na(GPP))/365) %>%
  mutate(above90 = case_when(n>=0.9 ~ 1,
                             TRUE ~ 0)) %>%
  ggplot(aes(Year, n, col = above90)) +
  geom_point() +
  facet_wrap(.~site_name) +
  ylim(0,1)

sites <- paste('nwis', c('02169000', '04137500', '04249000', '05435943', 
                         '08180700', '08211200', '10129900', '10133800'), 
               sep = '_')

png('figures/potential_gauss_datasets.png',
    width = 8, height = 6, res = 300, units = 'in')
dat %>% filter(site_name %in% sites) %>%
  group_by(site_name, Year) %>%
  summarize(n = sum(!is.na(GPP))/365) %>%
  mutate(per_miss = factor(case_when(n>=0.9 ~ '<10%',
                             n>=0.85 ~ '<15%',
                             TRUE ~ NA))) %>%
  ggplot(aes(Year, n, col = per_miss)) +
  geom_point() +
  facet_wrap(.~site_name) +
  ylim(0,1)
dev.off()


dat_sub <- dat %>% filter(site_name %in% sites) %>%
  group_by(site_name, Year) %>%
  summarize(n = sum(!is.na(GPP))/365) %>%
  filter(n>=0.85) %>%
  mutate(site_year = paste(site_name, Year, sep = '_')) %>%
  filter(!site_year %in% c('nwis_08180700_2013', 'nwis_02169000_2016', 
                           'nwis_04137500_2016', 'nwis_08211200_2016')) 
  
dd <- dat %>% 
  mutate(site_year = paste(site_name, Year, sep = '_')) %>%
  filter(site_year %in% dat_sub$site_year)
  
# edit dataframe to have relevant covariates:
glimpse(dd)
dd <- dd %>%
  group_by(site_name) %>%
  mutate(light = Stream_PAR / max(Stream_PAR, na.rm = T)) %>%
  select(site_name, long_name, date, Year, DOY, 
         GPP = GPP_filled, ER = ER_filled,
         light, Q = discharge) 

write_csv(dd, 'data/NWIS_MissingTS_subset_new.csv')
