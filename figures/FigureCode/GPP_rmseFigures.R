#####
# making figures for GPP forecast RMSE 
#####

# load packages
library(tidyverse)


# read in prediction data ------------------------------------------------------------
# read in real, complete dataset
realData <- read.csv("./data/pine_river_data_prepped.csv") %>% 
  mutate(date = lubridate::as_datetime(date))

# MAR arima data
MAR_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date))
# MAR brms data
MAR_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_FORECASTpreds_normPriorNB.csv") %>% 
  mutate(date = lubridate::as_datetime(date))

# MNAR arima data
MNAR_arima <- read.csv("./data/model_results/gauss_real_MNAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")
# MNAR brms data
MNAR_brms <- read.csv("./data/model_results/gauss_real_MNAR_brms_FORECASTpreds_normPriorNB.csv")%>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")

# join all data together
allDat <- rbind(MAR_arima, MAR_brms, MNAR_arima, MNAR_brms)

# remove runs w/ no missingness (missingprop_autocor = "y)
allDat <- allDat %>% 
  filter(missingprop_autocor != "y")
# reformat data -----------------------------------------------------------
# retrieve proportion missing and amount autocor from names
allDat <- allDat %>% 
  mutate(propMiss = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,2]),
         amtAutoCorr = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,4]))
allDat[allDat$missingness == "MNAR", "amtAutoCorr"] <- 0
# assign "bins" of autocorrelation
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr <=0.3 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_lowAutoCorr"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr > 0.3 & allDat$amtAutoCorr < 0.6 &!is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_medAutoCorr"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr  >= 0.6 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_highAutoCorr"

# calculate RMSE  ---------------------------------------------------------
RMSE <- allDat %>% 
  group_by(missingprop_autocor, missingness, type, propMiss, amtAutoCorr) %>% 
  summarise(RMSE = sqrt(mean((Estimate - GPP)^2, na.rm = TRUE)))

# add to all data df
allDat <- allDat %>% 
  left_join(RMSE)


# Make figure of predictions ----------------------------------------------

allDat_fig <- allDat 
# calculate low, med, and high autocorr
  
ggplot() + 
  facet_grid(.~as.factor(missingness)) +
  geom_line(data = realData[lubridate::month(realData$date) %in% c(lubridate::month(1:11)),], aes(x = date, y = GPP)) + 
  geom_line(data = allDat, aes(x = date, y = Estimate, col = type)) +
  theme_classic()


# RMSE figure -------------------------------------------------------------
ggplot(data = RMSE) +
  facet_grid(.~missingness) +
  geom_point(aes(x = propMiss, y = RMSE, col = type), alpha = .5) +
  geom_smooth(aes(x = propMiss, y = RMSE, col = type), method = "lm") +
  theme_classic()


      