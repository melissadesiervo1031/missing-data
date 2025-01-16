#####
# making figures for GPP forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)
library(magick)
library("grid")
library("ggplotify")

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
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr <=0.3 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR: low autocorrelation"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr > 0.3 & allDat$amtAutoCorr < 0.6 &!is.na(allDat$amtAutoCorr), "missingness"] <- "MAR: medium autocorrelation"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr  >= 0.6 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MAR: high autocorrelation"
allDat[allDat$missingness == "MNAR", "missingness"] <- "Missing NOT at Random"
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

## fix issue w/ some runs not having real GPP values... can get from other runs (every date should have the same "real" values)
allDat$GPP <- round(allDat$GPP, 9)
realGPP <- unique(allDat[,c("date", "GPP")])

# add back real GPP values to the main data.frame 

test <- allDat %>% 
  select(-GPP) %>% 
  left_join(realGPP, by = "date")


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

# # calculate low, med, and high autocorr
# (auSableTS_fig <-
#   ggplot(realData) +
#   geom_line(aes(x = date, y = GPP)))
# 
# # ggplot(allDat_fig[allDat_fig$missingness == "MAR_highAutoCorr" & allDat_fig$type == "Kalman Filter",]) + 
# #   geom_line(aes(x = date, y = Estimate_mean, col = propMiss_fac))
# 
# (forecastFig <- ggplot() +
#   facet_grid(. ~ as.factor(type) ~ as.factor(missingness)) +
#   geom_line(data = realData[lubridate::year(realData$date) == 2014,]
#             , aes(x = date, y = GPP), col = "grey60") +
#   geom_line(data = allDat_fig, aes(x = date, y = Estimate_mean, col = propMiss_fac, group = propMiss_fac), alpha = .8) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   )
# # save figure
# png(file = "./figures/forecastAccuracy_gaussian_auSable.png", width = 9, height = 6, units = "in", res = 700)
# forecastFig
# dev.off()
# 
# # Make figure of predictions aucorr x missingness ----------------------------------------------
# (forecastFig_2 <- ggplot() +
#     facet_grid(. ~ as.factor(propMiss_fac) ~ as.factor(missingness)) +
#     geom_line(data = realData[lubridate::year(realData$date) == 2014,]
#               , aes(x = date, y = GPP), col = "grey60") +
#     geom_line(data = allDat_fig, aes(x = date, y = Estimate_mean, col = type, group = type), alpha = .8) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
#                         labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
#  
# )
# # save figure
# png(file = "./figures/forecastAccuracy_v2_gaussian_auSable.png", width = 9, height = 6, units = "in", res = 700)
# forecastFig_2
# dev.off()

# RMSE figure -------------------------------------------------------------
(rmse_fig <- ggplot(data = RMSE) +
  facet_grid(.~missingness) +
  geom_point(aes(x = jitter(propMiss_bin, factor = 1.5), y = RMSE, col = type), alpha = .3) +
  geom_smooth(aes(x = propMiss_bin, y = RMSE, col = type), method = "lm", se = FALSE) +
  theme_classic() +
  #ylim(c(0,1.25)) + 
   ggplot2::labs(x = "Proportion of Missing Data", 
                 y = "Root Mean Square Error (RMSE)")+
    scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imp."), 
                          )  +
   guides(col = guide_legend(title = "Model Type"))
 )

# Data deletion simple: #E69F00 (orange)
#   Data deletion complete:  #D55E00  (red)
#   Multiple imputation: #009E73 (green)
#   Kalman filter: #0072B2 (blue)
#   Data augmentation: #CC79A7 (pink)
#  Expectation maximization: #BBBBBB (gray, rather than black)

# save figure
png(file = "./figures/RMSE_gaussian_auSable.png", width = 11, height = 5, units = "in", res = 700)
rmse_fig
dev.off()

# RMSE figure only med. autocorr and MNAR-------------------------------------------------------------
(rmse_fig2 <- RMSE %>% 
   filter(missingness %in% c("Missing NOT at Random", "MAR: medium autocorrelation")) %>%
   ungroup() %>% 
   ggplot() +
   facet_grid(.~missingness) +
   geom_point(aes(x = jitter(propMiss_bin, factor = 1.5), y = RMSE, col = type), alpha = .3) +
   geom_smooth(aes(x = propMiss_bin, y = RMSE, col = type), method = "lm", se = FALSE) +
   theme_classic() +
   #ylim(c(0,1.25)) + 
   ggplot2::labs(x = "Proportion of Missing Data", 
                 y = "Root Mean Square Error (RMSE)")+
   scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
                        labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imp."), 
   )  +
   guides(col = guide_legend(title = "Model Type"))
)

# Data deletion simple: #E69F00 (orange)
#   Data deletion complete:  #D55E00  (red)
#   Multiple imputation: #009E73 (green)
#   Kalman filter: #0072B2 (blue)
#   Data augmentation: #CC79A7 (pink)
#  Expectation maximization: #BBBBBB (gray, rather than black)

# save figure
png(file = "./figures/RMSE_medAutoCorrAndMNAR_gaussian_auSable.png", width = 7, height = 5, units = "in", res = 700)
rmse_fig2
dev.off()


# # RMSE figure - boxplot style -------------------------------------------------------------
# (rmse_boxFig <- ggplot(data = RMSE) +
#    facet_grid(.~missingness) +
#    geom_boxplot(aes(x = as.factor(propMiss_bin), y = RMSE, fill = type), alpha = .7) +
#    #geom_smooth(aes(x = propMiss, y = RMSE, col = type), method = "lm", se = FALSE) +
#    theme_classic() +
#    #ylim(c(0,1.25)) + 
#    scale_fill_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
#                         labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
# )
# # save figure
# png(file = "./figures/RMSE_boxplot_gaussian_auSable.png", width = 8, height = 4, units = "in", res = 700)
# rmse_boxFig
# dev.off()

# RMSE figure - lines with error style -------------------------------------------------------------
# reformat data
# get rid of RMSE low and high autocorrelation
RMSE <- RMSE %>% 
  filter(missingness %in% c("Missing NOT at Random", "MAR: medium autocorrelation")) %>% 
  mutate(missingness = str_replace(missingness, pattern = "MAR: medium autocorrelation", replacement = "Missing at Random"))

RMSE_errorBar <- RMSE %>% 
  # filter(missingness %in% c("Missing NOT at Random", "MAR: medium autocorrelation"))%>% 
  # mutate(missingness = str_replace(missingness, pattern = "MAR: medium autocorrelation", replacement = "Missing at Random")) %>%
  group_by(missingness, type, propMiss_bin) %>% 
  summarize(RMSE_mean = mean(RMSE), 
            RMSE_sd = sd(RMSE)) %>% 
  rowwise() %>% 
  mutate(low_95CI = (RMSE_mean - RMSE_sd*1.96), 
         high_95CI = (RMSE_mean + RMSE_sd*1.96))


## order the missingness type factor
RMSE_errorBar <- RMSE_errorBar %>% 
  mutate(type = factor(type, levels =c("dropNA_simple", "dropNA_complete", "Multiple Imputations",
                       "Kalman Filter", "brms"), ordered = TRUE))
RMSE <- RMSE %>% 
  mutate(type = factor(type, levels =c("dropNA_simple", "dropNA_complete", "Multiple Imputations",
                                       "Kalman Filter", "brms"), ordered = TRUE))

(rmse_fig2 <- RMSE %>% 
    ggplot() +
    facet_grid(.~missingness) +
    geom_point(aes(x = jitter(propMiss_bin, factor = 1.5), y = RMSE, col = type), alpha = .3) +
    geom_smooth(aes(x = propMiss_bin, y = RMSE, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    #ylim(c(0,1.25)) + 
    ggplot2::labs(x = "Proportion of Missing Data", 
                  y = "Root Mean Square Error (RMSE)") +
    scale_color_discrete(
                         type = c( "#E69F00","#D55E00","#009E73", "#0072B2", "#CC79A7"),
                         labels = c("Data Deletion-Simple", "Data Deletion-Complete", "Multiple Imp.", "Kalman Filter", "Data Augmentation")) + 
    guides(col = guide_legend(title = "Model Type"))
)


(rmse_lineErrorBar <- ggplot(data = RMSE_errorBar) +
   facet_grid(.~missingness) +
    geom_linerange(aes(x = propMiss_bin, ymin = low_95CI, ymax = high_95CI, color = type), alpha = .7, position = position_dodge(width = .1)) +
   geom_point(aes(x = propMiss_bin, y = RMSE_mean, color = type), alpha = .7, position = position_dodge(width = .1)) +
    geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
   theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
   #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
    labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imp."), 
)  +
  guides(col = guide_legend(title = "Model Type"))
)

(rmse_NoLineErrorBar <- ggplot(data = RMSE_errorBar) +
    facet_grid(.~missingness) +
    geom_linerange(aes(x = propMiss_bin, ymin = low_95CI, ymax = high_95CI, color = type), alpha = .7, position = position_dodge(width = .1)) +
    geom_point(aes(x = propMiss_bin, y = RMSE_mean, color = type), alpha = .7, position = position_dodge(width = .1)) +
    #geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imp."), 
    )  +
    guides(col = guide_legend(title = "Model Type", position = "top", direction = "vertical", nrow = 2))
)
## get complete time series figure 
aus <- read.csv("./data/au_sable_river_prepped.csv", header=TRUE)

ausNew <- aus %>% 
  mutate(date = as.POSIXct(date)) 

(tsFigGrob <- 
ggplot() + 
  geom_rect(aes(xmin = as.POSIXct("2014-01-01T00:00:00Z"), xmax = as.POSIXct("2014-12-31T00:00:00Z"), 
                ymin = -6, ymax = 3.5), fill = "grey80") +
  geom_line(data = ausNew, aes(x = date, y = GPP)) + 
  geom_rug(data = ausNew[is.na(ausNew$GPP),], aes(date), col = "red") +
  theme_classic() + 
  labs(x = "Year", y = "GPP (scaled)"))

tsPlusRmse_withLine <- ggpubr::ggarrange(tsFigGrob, rmse_lineErrorBar, ncol = 1, nrow = 2, 
                  heights = c(.5, 1), legend = "bottom", labels = c("A", "B"))

tsPlusRmse_NoLine <- ggpubr::ggarrange(tsFigGrob, rmse_NoLineErrorBar, ncol = 1, nrow = 2, 
                                         heights = c(.5, 1), legend = "bottom", labels = c("A", "B"))

# save figure
png(file = "./figures/RMSE_FullFigure_lineWithErrorBar_gaussian_auSable.png", width = 9, height = 8, units = "in", res = 700)
tsPlusRmse_withLine
dev.off()

png(file = "./figures/RMSE_FullFigure_NoLineWithErrorBar_gaussian_auSable.png", width = 6.5, height = 8, units = "in", res = 700)
tsPlusRmse_NoLine
dev.off()

# coverage figure  --------------------------------------------------------
# calculate coverage

## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
figDat_covTemp <- allDat %>% 
  mutate(autoCor = round(amtAutoCorr_bin, 1), 
         amtMiss = round(propMiss_bin, 1),
         # calculate coverage
         CI95_lower = Estimate - 1.96 * Est.Error,
         CI95_upper = Estimate + 1.96 * Est.Error
         ) %>% 
  filter(!is.na(Est.Error))


figDat_covTemp$coverage <- c(figDat_covTemp$GPP >= figDat_covTemp$CI95_lower & 
                               figDat_covTemp$GPP <= figDat_covTemp$CI95_upper)

fitDat_cov <- figDat_covTemp %>% 
  group_by(missingprop_autocor, missingness, type, run_no, autoCor, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage, na.rm = TRUE), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)

# make figure
(coverage_fig <- ggplot(data = fitDat_cov) +
    facet_grid(.~missingness, scales = "free") +
    geom_boxplot(aes(x = as.factor(amtMiss), y = coveragePerc, col = type), alpha = .5) +
    #geom_smooth(aes(x = jitter(propMiss_bin), y = RMSE, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    #ylim(c(0,1.25)) + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
)


(rmse_boxFig <- ggplot(data = RMSE) +
    facet_grid(.~missingness) +
    geom_boxplot(aes(x = as.factor(propMiss_bin), y = RMSE, fill = type), alpha = .7) +
    #geom_smooth(aes(x = propMiss, y = RMSE, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    #ylim(c(0,1.25)) + 
    scale_fill_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Del.-Complete", "Data Aug.", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) 
)
