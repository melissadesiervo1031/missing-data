#####
# making figures for GPP forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)


# read in prediction data for au-sable river ------------------------------------------------------------
# read in real, complete dataset
realData <- read.csv("./data/au_sable_river_prepped.csv") %>% 
  mutate(date = lubridate::as_datetime(date))

# MAR arima data
MAR_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/au_sable/gauss_real_MAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date))
# MAR brms data
MAR_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/auSable/gauss_auSable_real_MAR_brms_FORECASTpreds.csv") %>% 
  mutate(date = lubridate::as_datetime(date)) %>% 
  select(-X)

# MNAR arima data
MNAR_arima <- read.csv("./data/model_results/gauss_real_MNAR_arima_modResults/au_sable/gauss_auSable_real_MNAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")
# MNAR brms data
MNAR_brms <- read.csv("./data/model_results/gauss_real_MNAR_brms_modResults/auSable/brmspreds.csv")%>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")

# join all data together
allDat <- rbind(MAR_arima, MAR_brms, MNAR_arima, MNAR_brms)

# remove rows w/ no missingness (missingprop_autocor = "y)
allDat <- allDat %>% 
  filter(missingprop_autocor != "y") %>% 
  filter(!is.na(Estimate)) # remove NAs in Estimate column (12/31)

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

allDat_fig <- allDat %>% 
  group_by(date, missingness, type, round(propMiss,1)) %>% 
  summarize(Estimate_mean = mean(Estimate),
            Estimate_sd = sd(Estimate),
            Est.Error_mean = mean(Est.Error)) %>% 
  rename(propMiss = "round(propMiss, 1)")
  
# calculate low, med, and high autocorr
(auSableTS_fig <- 
  ggplot(realData) + 
  geom_line(aes(x = date, y = GPP)))

(forecastFig <- ggplot() + 
  facet_grid(.~as.factor(missingness) ~ as.factor(type)) +
  geom_line(data = realData[lubridate::year(realData$date) == 2014,]
            , aes(x = date, y = GPP), col = "grey60") + 
  geom_line(data = allDat_fig, aes(x = date, y = Estimate_mean, col = propMiss, group = propMiss), alpha = .8) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
# save figure
png(file = "./figures/forecastAccuracy_gaussian_auSable.png", width = 9, height = 6, units = "in", res = 700)
forecastFig
dev.off()
# RMSE figure -------------------------------------------------------------
(rmse_fig <- ggplot(data = RMSE) +
  facet_grid(.~missingness) +
  geom_point(aes(x = propMiss, y = RMSE, col = type), alpha = .5) +
  geom_smooth(aes(x = propMiss, y = RMSE, col = type), method = "lm", se = FALSE) +
  theme_classic() +
  ylim(c(0,1.25)) + 
   scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
 )
# save figure
png(file = "./figures/RMSE_gaussian_auSable.png", width = 8, height = 4, units = "in", res = 700)
rmse_fig
dev.off()

#///////////////////////////
# read in prediction data for badger creek river ------------------------------------------------------------
# read in real, complete dataset
  realData <- read.csv("./data/badger_mill_creek_prepped.csv") %>% 
  mutate(date = lubridate::as_datetime(date))

# MAR arima data
MAR_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/badger_mill/gauss_badger_real_MAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date))
# MAR brms data
MAR_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/badgerMill/gauss_badger_real_MAR_brms_FORECASTpreds.csv") %>% 
  mutate(date = lubridate::as_datetime(date)) %>% 
  select(-X)

# MNAR arima data
MNAR_arima <- read.csv("./data/model_results/gauss_real_MNAR_arima_modResults/badger_mill/gauss_badger_real_MNAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")
# MNAR brms data
MNAR_brms <- read.csv("./data/model_results/gauss_real_MNAR_brms_modResults/badgerMill/brmspreds.csv")%>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")

# join all data together
allDat <- rbind(MAR_arima, MAR_brms, MNAR_arima, MNAR_brms)

# remove rows w/ no missingness (missingprop_autocor = "y)
allDat <- allDat %>% 
  filter(missingprop_autocor != "y") %>% 
  filter(!is.na(Estimate)) # remove NAs in Estimate column (12/31)

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

allDat_fig <- allDat %>% 
  group_by(date, missingness, type, round(propMiss,1)) %>% 
  summarize(Estimate_mean = mean(Estimate),
            Estimate_sd = sd(Estimate),
            Est.Error_mean = mean(Est.Error)) %>% 
  rename(propMiss = "round(propMiss, 1)")

# calculate low, med, and high autocorr
(badgerTS_fig <- 
    ggplot(realData) + 
    geom_line(aes(x = date, y = GPP)))

(forecastFig <- ggplot() + 
    facet_grid(.~as.factor(missingness) ~ as.factor(type)) +
    geom_line(data = realData[lubridate::year(realData$date) == 2015,]
              , aes(x = date, y = GPP), col = "grey60") + 
    geom_line(data = allDat_fig, aes(x = date, y = Estimate_mean, col = propMiss, group = propMiss), alpha = .8) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
# save figure
png(file = "./figures/forecastAccuracy_gaussian_badger.png", width = 9, height = 6, units = "in", res = 700)
forecastFig
dev.off()
# RMSE figure -------------------------------------------------------------
(rmse_fig <- ggplot(data = RMSE) +
   facet_grid(.~missingness) +
   geom_point(aes(x = propMiss, y = RMSE, col = type), alpha = .5) +
   geom_smooth(aes(x = propMiss, y = RMSE, col = type), method = "lm", se = FALSE) +
   theme_classic() +
   ylim(c(0,1.25)) + 
   scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
)
# save figure
png(file = "./figures/RMSE_gaussian_badger.png", width = 8, height = 4, units = "in", res = 700)
rmse_fig
dev.off()
#///////////////////////////
# figure of mean predictions across all amounts missingness for each type  --------
second_fig <- allDat_fig %>% 
  filter(
         propMiss <=0.5)
ggplot() + 
  facet_grid(.~as.factor(missingness)~ as.factor(propMiss)) +
  #geom_ribbon(data = second_fig, aes(x = date, ymin = Estimate_mean - 1.96 * Est.Error_mean, ymax = Estimate_mean + 1.96 * Est.Error_mean, fill = type, group =type), alpha = .3) +
  geom_line(data = realData[lubridate::month(realData$date) %in% c(lubridate::month(11:12)),], aes(x = date, y = GPP)) + 
  geom_line(data = second_fig, aes(x = date, y = Estimate_mean, col = type, group = type)) +
  theme_classic() 

# figure of mean predictions across all amounts of missingness for each type  --------
third_fig <- allDat_fig %>% 
  filter(propMiss <=0.5) %>% 
  group_by(date, missingness, type) %>% 
  summarize(Estimate_mean = mean(Estimate_mean),
            Est.Error_mean = mean(Est.Error_mean))
  
ggplot() + 
  facet_grid(.~as.factor(missingness)~ as.factor(type)) +
  geom_ribbon(data = third_fig, aes(x = date, ymin = Estimate_mean - 1.96 * Est.Error_mean, ymax = Estimate_mean + 1.96 * Est.Error_mean, fill = type, group =type), alpha = .3) +
  geom_line(data = realData[lubridate::month(realData$date) %in% c(lubridate::month(11:12)),], aes(x = date, y = GPP)) + 
  geom_line(data = third_fig, aes(x = date, y = Estimate_mean, col = type, group = type)) +
  theme_classic() 


# figure of preds vs. full real time series -------------------------------
allDat_all <- allDat %>% 
  mutate(ID = paste0(type, missingness, propMiss, amtAutoCorr))

# calculate low, med, and high autocorr

( ggplot() + 
    geom_line(data = realData, aes(x = date, y = GPP)) + 
    geom_line(data = allDat_all, aes(x = date, y = Estimate, col = type, group = ID), alpha = .3) +
    theme_classic() + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
  
)



      