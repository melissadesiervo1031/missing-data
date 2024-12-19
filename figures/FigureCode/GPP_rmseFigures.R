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
MAR_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/au_sable/gauss_auSable_real_MAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X, -...1, -...11, -...12, -...13) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date))
# MAR brms data
MAR_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/auSable/gauss_auSable_real_MAR_brms_FORECASTpreds.csv") %>% 
  mutate(date = lubridate::as_datetime(date)) %>% 
  select(-...1, -X, -...12)

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
allDat_temp <- rbind(MAR_arima, MAR_brms, MNAR_arima, MNAR_brms)

# remove rows w/ no missingness (missingprop_autocor = "y)
allDat_temp <- allDat_temp %>% 
  filter(missingprop_autocor != "y") %>% 
  filter(!is.na(Estimate)) # remove NAs in Estimate column (12/31)

# reformat data -----------------------------------------------------------
# retrieve proportion missing and amount autocor from names
# fix issue with autocorrelation calculation
allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))] <- 
  str_sub(allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))], 
          start = 1, end = 30)
allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))] <- 
  str_sub(allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))], 
          start = 1, end = 16)

# bin autocorrelation
allDat <- allDat_temp %>% 
  mutate(propMiss = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,2]),
         amtAutoCorr = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,4]))
allDat[allDat$missingness == "MNAR", "amtAutoCorr"] <- 0
# assign "bins" of autocorrelation
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr <=0.3 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_lowAutoCorr"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr > 0.3 & allDat$amtAutoCorr < 0.6 &!is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_medAutoCorr"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr  >= 0.6 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_highAutoCorr"

# bin autocorr and propMiss 
allDat$propMiss_bin <- NA
allDat$amtAutoCorr_bin  <- NA

allDat[allDat$propMiss == 0 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0
allDat[allDat$propMiss > 0 & allDat$propMiss < 0.3 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.2
allDat[allDat$propMiss >= 0.3 & allDat$propMiss < 0.5 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.4
allDat[allDat$propMiss >= 0.5 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0.6

allDat[allDat$amtAutoCorr == 0 &  !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0
allDat[allDat$amtAutoCorr > 0 & allDat$amtAutoCorr < 0.375 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.25
allDat[allDat$amtAutoCorr >= 0.375 & allDat$amtAutoCorr < 0.625 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.5
allDat[allDat$amtAutoCorr >= 0.625 &  !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.75
  
# calculate RMSE  ---------------------------------------------------------
RMSE <- allDat %>% 
  group_by(missingprop_autocor, missingness, type, propMiss_bin, amtAutoCorr_bin) %>% 
  summarise(RMSE = sqrt(mean((Estimate - GPP)^2, na.rm = TRUE))) 

# add to all data df
allDat <- allDat %>% 
  left_join(RMSE)

# Make figure of predictions ----------------------------------------------
allDat_fig <- allDat %>%
  group_by(date, missingness, type, propMiss_bin) %>%
  summarize(Estimate_mean = mean(Estimate),
            Estimate_sd = sd(Estimate),
            Est.Error_mean = mean(Est.Error)) %>%
  rename(propMiss = "propMiss_bin") %>% 
  mutate(propMiss_fac = factor(propMiss))

# calculate low, med, and high autocorr
(auSableTS_fig <-
  ggplot(realData) +
  geom_line(aes(x = date, y = GPP)))

# ggplot(allDat_fig[allDat_fig$missingness == "MAR_highAutoCorr" & allDat_fig$type == "Kalman Filter",]) + 
#   geom_line(aes(x = date, y = Estimate_mean, col = propMiss_fac))

(forecastFig <- ggplot() +
  facet_grid(. ~ as.factor(type) ~ as.factor(missingness)) +
  geom_line(data = realData[lubridate::year(realData$date) == 2014,]
            , aes(x = date, y = GPP), col = "grey60") +
  geom_line(data = allDat_fig, aes(x = date, y = Estimate_mean, col = propMiss_fac, group = propMiss_fac), alpha = .8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
# save figure
png(file = "./figures/forecastAccuracy_gaussian_auSable.png", width = 9, height = 6, units = "in", res = 700)
forecastFig
dev.off()

# Make figure of predictions aucorr x missingness ----------------------------------------------
(forecastFig_2 <- ggplot() +
    facet_grid(. ~ as.factor(propMiss_fac) ~ as.factor(missingness)) +
    geom_line(data = realData[lubridate::year(realData$date) == 2014,]
              , aes(x = date, y = GPP), col = "grey60") +
    geom_line(data = allDat_fig, aes(x = date, y = Estimate_mean, col = type, group = type), alpha = .8) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
 
)
# save figure
png(file = "./figures/forecastAccuracy_v2_gaussian_auSable.png", width = 9, height = 6, units = "in", res = 700)
forecastFig_2
dev.off()

# RMSE figure -------------------------------------------------------------
(rmse_fig <- ggplot(data = RMSE) +
  facet_grid(.~missingness) +
  geom_point(aes(x = jitter(propMiss_bin), y = RMSE, col = type), alpha = .5) +
  geom_smooth(aes(x = jitter(propMiss_bin), y = RMSE, col = type), method = "lm", se = FALSE) +
  theme_classic() +
  #ylim(c(0,1.25)) + 
   scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
 )
# save figure
png(file = "./figures/RMSE_gaussian_auSable.png", width = 8, height = 4, units = "in", res = 700)
rmse_fig
dev.off()

# RMSE figure - boxplot style -------------------------------------------------------------
(rmse_boxFig <- ggplot(data = RMSE) +
   facet_grid(.~missingness) +
   geom_boxplot(aes(x = as.factor(propMiss_bin), y = RMSE, fill = type), alpha = .7) +
   #geom_smooth(aes(x = propMiss, y = RMSE, col = type), method = "lm", se = FALSE) +
   theme_classic() +
   #ylim(c(0,1.25)) + 
   scale_fill_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
)
# save figure
png(file = "./figures/RMSE_boxplot_gaussian_auSable.png", width = 8, height = 4, units = "in", res = 700)
rmse_boxFig
dev.off()

# # read in prediction data for badger mill creek 
# # read in real, complete dataset
# realData <- read.csv("./data/badger_mill_creek_prepped.csv") %>% 
#   mutate(date = lubridate::as_datetime(date))
# 
# # MAR arima data
# MAR_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/badger_mill/gauss_badger_real_MAR_arima_FORECASTpreds.csv") %>% 
#   rename(Estimate = pred, Est.Error = se) %>% 
#   mutate(Q2.5 = NA, Q97.5 = NA) %>% 
#   select(-CurSim, -X, -...1, -...11, -...12) %>% 
#   relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
#   mutate(date = lubridate::as_datetime(date))
# # MAR brms data
# MAR_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/badgerMill/gauss_badger_real_MAR_brms_FORECASTpreds.csv") %>% 
#   mutate(date = lubridate::as_datetime(date)) %>% 
#   select(-...1, -X, -...12)
# 
# # MNAR arima data
# MNAR_arima <- read.csv("./data/model_results/gauss_real_MNAR_arima_modResults/badger_mill/gauss_badger_real_MNAR_arima_FORECASTpreds.csv") %>% 
#   rename(Estimate = pred, Est.Error = se) %>% 
#   mutate(Q2.5 = NA, Q97.5 = NA) %>% 
#   select(-CurSim, -X) %>% 
#   relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
#   mutate(date = lubridate::as_datetime(date),
#          missingness = "MNAR")
# # MNAR brms data
# MNAR_brms <- read.csv("./data/model_results/gauss_real_MNAR_brms_modResults/badgerMill/brmspreds.csv")%>% 
#   mutate(date = lubridate::as_datetime(date),
#          missingness = "MNAR")
# 
# # join all data together
# allDat_temp <- rbind(MAR_arima, MAR_brms, MNAR_arima, MNAR_brms) %>% 
#   filter(lubridate::year(date) == 2015)
# 
# # remove rows w/ no missingness (missingprop_autocor = "y)
# allDat_temp <- allDat_temp %>% 
#   filter(missingprop_autocor != "y") %>% 
#   filter(!is.na(Estimate)) # remove NAs in Estimate column (12/31)
# 
# # reformat data --
# # retrieve proportion missing and amount autocor from names
# # fix issue with autocorrelation calculation
# allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))] <- 
#   str_sub(allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))], 
#           start = 1, end = 30)
# allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))] <- 
#   str_sub(allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))], 
#           start = 1, end = 16)
# 
# # bin autocorrelation
# allDat <- allDat_temp %>% 
#   mutate(propMiss = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,2]),
#          amtAutoCorr = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,4]))
# allDat[allDat$missingness == "MNAR", "amtAutoCorr"] <- 0
# # assign "bins" of autocorrelation
# allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr <=0.3 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_lowAutoCorr"
# allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr > 0.3 & allDat$amtAutoCorr < 0.6 &!is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_medAutoCorr"
# allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr  >= 0.6 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR_highAutoCorr"
# 
# # bin autocorr and propMiss 
# allDat$propMiss_bin <- NA
# allDat$amtAutoCorr_bin  <- NA
# 
# allDat[allDat$propMiss == 0 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0
# allDat[allDat$propMiss > 0 & allDat$propMiss < 0.3 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.2
# allDat[allDat$propMiss >= 0.3 & allDat$propMiss < 0.5 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.4
# allDat[allDat$propMiss >= 0.5 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0.6
# 
# allDat[allDat$amtAutoCorr == 0 &  !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0
# allDat[allDat$amtAutoCorr > 0 & allDat$amtAutoCorr < 0.375 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.25
# allDat[allDat$amtAutoCorr >= 0.375 & allDat$amtAutoCorr < 0.625 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.5
# allDat[allDat$amtAutoCorr >= 0.625 &  !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.75
# 
# # calculate RMSE  
# RMSE <- allDat %>% 
#   group_by(missingprop_autocor, missingness, type, propMiss_bin, amtAutoCorr_bin) %>% 
#   summarise(RMSE = sqrt(mean((Estimate - GPP)^2, na.rm = TRUE))) 
# 
# # add to all data df
# allDat <- allDat %>% 
#   left_join(RMSE)
# 
# # Make figure of predictions 
# allDat_fig <- allDat %>%
#   group_by(date, missingness, type, propMiss_bin) %>%
#   summarize(Estimate_mean = mean(Estimate),
#             Estimate_sd = sd(Estimate),
#             Est.Error_mean = mean(Est.Error)) %>%
#   rename(propMiss = "propMiss_bin") %>% 
#   mutate(propMiss_fac = factor(propMiss))
# 
# (forecastFig <- ggplot() +
#     facet_grid(. ~ as.factor(type) ~ as.factor(missingness)) +
#     geom_line(data = realData[lubridate::year(realData$date) == 2015,]
#               , aes(x = date, y = GPP), col = "grey60") +
#     geom_line(data = allDat_fig[lubridate::year(allDat_fig$date) == 2015,], aes(x = date, y = Estimate_mean, col = propMiss_fac, group = propMiss_fac), alpha = .8) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# )
# # save figure
# png(file = "./figures/forecastAccuracy_gaussian_badger.png", width = 9, height = 6, units = "in", res = 700)
# forecastFig
# dev.off()
# 
# # Make figure of predictions aucorr x missingness 
# (forecastFig_2 <- ggplot() +
#    facet_grid(. ~ as.factor(propMiss_fac) ~ as.factor(missingness)) +
#    geom_line(data = realData[lubridate::year(realData$date) == 2015,]
#              , aes(x = date, y = GPP), col = "grey60") +
#    geom_line(data = allDat_fig[lubridate::year(allDat_fig$date) == 2015,], aes(x = date, y = Estimate_mean, col = type, group = type), alpha = .8) +
#    theme_classic() +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
#                         labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
#  
# )
# # save figure
# png(file = "./figures/forecastAccuracy_v2_gaussian_badger.png", width = 9, height = 6, units = "in", res = 700)
# forecastFig_2
# dev.off()
# 
# # RMSE figure -
# (rmse_fig <- ggplot(data = RMSE) +
#    facet_grid(.~missingness) +
#    geom_point(aes(x = jitter(propMiss_bin), y = RMSE, col = type), alpha = .5) +
#    geom_smooth(aes(x = jitter(propMiss_bin), y = RMSE, col = type), method = "lm", se = FALSE) +
#    theme_classic() +
#    #ylim(c(0,1.25)) + 
#    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
#                         labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
# )
# # save figure
# png(file = "./figures/RMSE_gaussian_badger.png", width = 8, height = 4, units = "in", res = 700)
# rmse_fig
# dev.off()
# 
# # RMSE figure - boxplot style -
# (rmse_boxFig <- ggplot(data = RMSE) +
#    facet_grid(.~missingness) +
#    geom_boxplot(aes(x = as.factor(propMiss_bin), y = RMSE, fill = type), alpha = .7) +
#    #geom_smooth(aes(x = propMiss, y = RMSE, col = type), method = "lm", se = FALSE) +
#    theme_classic() +
#    #ylim(c(0,1.25)) + 
#    scale_fill_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
#                        labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
# )
# # save figure
# png(file = "./figures/RMSE_boxplot_gaussian_badger.png", width = 8, height = 4, units = "in", res = 700)
# rmse_boxFig
# dev.off()
# 
# # 
# # #///////////////////////////
# # # figure of mean predictions across all amounts missingness for each type  
# # second_fig <- allDat_fig %>% 
# #   filter(
# #          propMiss <=0.5)
# # ggplot() + 
# #   facet_grid(.~as.factor(missingness)~ as.factor(propMiss)) +
# #   #geom_ribbon(data = second_fig, aes(x = date, ymin = Estimate_mean - 1.96 * Est.Error_mean, ymax = Estimate_mean + 1.96 * Est.Error_mean, fill = type, group =type), alpha = .3) +
# #   geom_line(data = realData[lubridate::month(realData$date) %in% c(lubridate::month(11:12)),], aes(x = date, y = GPP)) + 
# #   geom_line(data = second_fig, aes(x = date, y = Estimate_mean, col = type, group = type)) +
# #   theme_classic() 
# # 
# # # figure of mean predictions across all amounts of missingness for each type  -
# # third_fig <- allDat_fig %>% 
# #   filter(propMiss <=0.5) %>% 
# #   group_by(date, missingness, type) %>% 
# #   summarize(Estimate_mean = mean(Estimate_mean),
# #             Est.Error_mean = mean(Est.Error_mean))
# #   
# # ggplot() + 
# #   facet_grid(.~as.factor(missingness)~ as.factor(type)) +
# #   geom_ribbon(data = third_fig, aes(x = date, ymin = Estimate_mean - 1.96 * Est.Error_mean, ymax = Estimate_mean + 1.96 * Est.Error_mean, fill = type, group =type), alpha = .3) +
# #   geom_line(data = realData[lubridate::month(realData$date) %in% c(lubridate::month(11:12)),], aes(x = date, y = GPP)) + 
# #   geom_line(data = third_fig, aes(x = date, y = Estimate_mean, col = type, group = type)) +
# #   theme_classic() 
# # 
# # 
# # # figure of preds vs. full real time series -
# # allDat_all <- allDat %>% 
# #   mutate(ID = paste0(type, missingness, propMiss, amtAutoCorr))
# # 
# # # calculate low, med, and high autocorr
# # 
# # ( ggplot() + 
# #     geom_line(data = realData, aes(x = date, y = GPP)) + 
# #     geom_line(data = allDat_all, aes(x = date, y = Estimate, col = type, group = ID), alpha = .3) +
# #     theme_classic() + 
# #     scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
# #                          labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
# #   
# # )
# # 
# # 
# # 
# #       