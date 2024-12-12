#####
# making figures for Population forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)
library(RColorBrewer)


# read in prediction data ------------------------------------------------------------
#allDat <- readRDS("./data/model_results/RickerForecast_resultTableAll.rds")
allDat <- readRDS("./data/model_results/RickerForecast_resultTableAll_10.rds")
#allDat <- readRDS("./data/model_results/RickerForecast_resultTableAll_20.rds")
allDat[is.nan(allDat$actAutoCorr_trim), "actAutoCorr_trim"] <- 0
# read in real data
realDat <- read.csv("./data/Wytham_tits.csv")
realDat$timeStep <- seq(1:59)
train_length=49

# remove runs where the truncated "missing data" ended in an NA (because we
# can't accurately intiate the model predicting the last five years in taht
# case)
allDat_new <- allDat[sapply(allDat$trimmed_ts, FUN = function(x) 
  length(x) ==train_length),]

# put predictions into a 'long' format
forecasts_long <- list_rbind(apply(allDat_new, MARGIN = 1, FUN = function(x) {
  startInd <- which(!is.na(x$forecasts$timeStep))[2]
  data.frame(
    "newName" = x$newName, "autoCorr" = x$actAutoCorr_trim, "propMiss" = x$actPropMiss_trim,
    "forecast_dropNA" = x$forecasts$dropNA_est[startInd:59], 
    "forecast_dropCC" = x$forecasts$dropCC_est[startInd:59], 
    "forecast_MI" = x$forecasts$MI_est[startInd:59], 
    "forecast_EM" = x$forecasts$EM_est[startInd:59], 
    "forecast_DA" = x$forecasts$DA_est[startInd:59],
    "timeStep" = x$forecasts$timeStep[startInd:59]
  )
}))
# bin according to autocorrelation
forecasts_long$autocorr_binned <- NA
forecasts_long[forecasts_long$autoCorr < 0.3, "autocorr_binned"] <- "low_autocorr"
forecasts_long[forecasts_long$autoCorr >= 0.3 & forecasts_long$autoCorr < 0.6, "autocorr_binned"] <- "med_autocorr"
forecasts_long[forecasts_long$autoCorr >=  0.6, "autocorr_binned"] <- "high_autocorr"
forecasts_long$autocorr_binned <- factor(forecasts_long$autocorr_binned, 
                                         levels = c("low_autocorr",  "med_autocorr",  "high_autocorr"))

# add back in year data
forecasts_long <- forecasts_long %>% 
  left_join(realDat[,c("Year", "timeStep")])

# make it even longer (so that there is only one column for predicted value,
# with the type of missingness in "missingness" column)
forecasts_long <- forecasts_long %>% 
  rename("dropNA_simple"="forecast_dropNA", "dropNA_CC"="forecast_dropCC", 
         "MI"="forecast_MI", "EM"="forecast_EM", "DA" ="forecast_DA") %>% 
  pivot_longer(cols = c("dropNA_simple", "dropNA_CC", "MI", "EM", "DA"), 
               values_to = "Estimate", names_to = "missingness")

# put RMSE data into it's own data frame
RMSE_df <- list_rbind(apply(allDat_new, MARGIN = 1, FUN = function(x) {
  data.frame(
    "newName" = x$newName, "autoCorr" = x$actAutoCorr_trim, "propMiss" = x$actPropMiss_trim,
    "RMSE" = x$RMSE, 
    "modelType" = names(x$RMSE)
  )
}))
rownames(RMSE_df) <- NULL
# bin according to autocorrelation
RMSE_df$autocorr_binned <- NA
RMSE_df[RMSE_df$autoCorr < 0.3, "autocorr_binned"] <- "low_autocorr"
RMSE_df[RMSE_df$autoCorr >= 0.3 & RMSE_df$autoCorr < 0.6, "autocorr_binned"] <- "med_autocorr"
RMSE_df[RMSE_df$autoCorr >=  0.6, "autocorr_binned"] <- "high_autocorr"
RMSE_df$autocorr_binned <- factor( RMSE_df$autocorr_binned, 
                                   levels = c("low_autocorr",  "med_autocorr",  "high_autocorr"))

RMSE_df$propMissCat=NA
RMSE_df$propMissCat[which(RMSE_df$propMiss==0)]=0
RMSE_df$propMissCat[which(RMSE_df$propMiss>=0.15&RMSE_df$propMiss<=0.25)]=0.2
RMSE_df$propMissCat[which(RMSE_df$propMiss>=0.35&RMSE_df$propMiss<=0.45)]=0.4
RMSE_df$propMissCat[which(RMSE_df$propMiss>=0.55&RMSE_df$propMiss<=0.65)]=0.6

# Plot RMSE against missingness -------------------------------------------
rmse_missingness_p <- ggplot(RMSE_df) +
  geom_point(aes(x = propMiss, y = RMSE, col = modelType), alpha = .5) +
  geom_smooth(aes(x = propMiss, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
  facet_wrap(~autocorr_binned, ) +
  scale_color_discrete(type = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
                       labels = c("Data Aug.","Data Del.-Complete",  "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
  theme_classic()

# save figure
png(file = "./figures/RMSE_poisson.png", width = 8, height = 4, units = "in", res = 700)
rmse_missingness_p
dev.off()


# Plot RMSE against missingness new -------------------------------------------
rmse_missingness_box <- ggplot(RMSE_df) +
  geom_boxplot(aes(x = factor(propMissCat), y = RMSE, fill = modelType), alpha = 0.7) +
  facet_wrap(~autocorr_binned) +
  scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
                    labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
  labs(
    x = "Proportion Missingness",
    y = "Root Mean Squared Error (RMSE)",
    fill = "Model Type"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(size = 12),  # Customize facet labels
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
  )

rmse_missingness_box


# Make figure of predictions ----------------------------------------------
forecasts_long$propMiss_binned <- NA
forecasts_long[forecasts_long$propMiss <0.2, "propMiss_binned"] <- "0-0.2"
forecasts_long[c(forecasts_long$propMiss >=0.2 & forecasts_long$propMiss < 0.4), "propMiss_binned"] <- "0.2-0.4"
forecasts_long[forecasts_long$propMiss >=0.4 & forecasts_long$propMiss < 0.6, "propMiss_binned"] <- "0.4-0.6"
forecasts_long[forecasts_long$propMiss >=0.6 & forecasts_long$propMiss < 0.8, "propMiss_binned"] <- "0.6-0.8"
forecasts_long[forecasts_long$propMiss >=0.8, "propMiss_binned"] <- "0.8-1"

forecasts_long$propMissCat=NA
forecasts_long$propMissCat[which(forecasts_long$propMiss==0)]=0
forecasts_long$propMissCat[which(forecasts_long$propMiss>=0.15&forecasts_long$propMiss<=0.25)]=0.2
forecasts_long$propMissCat[which(forecasts_long$propMiss>=0.35&forecasts_long$propMiss<=0.45)]=0.4
forecasts_long$propMissCat[which(forecasts_long$propMiss>=0.55&forecasts_long$propMiss<=0.65)]=0.6

forecasts_avg <- forecasts_long %>% 
  #filter(Year > 2013) %>% 
  group_by(Year, missingness, autocorr_binned, propMissCat) %>% 
  summarize(Estimate_mean = mean(Estimate),
            Estimate_sd = sd(Estimate)) #%>% 
#,Est.Error_mean = mean(Est.Error)) %>% 
#rename(propMiss = propMiss_binned)

# calculate low, med, and high autocorr
### set color palettes ------------------------------------------------------
palCont <- rev(RColorBrewer::brewer.pal(n = 4, name = "Reds"))

(forecastFig_p <- ggplot() + 
    facet_grid(.~as.factor(missingness) ~ as.factor(autocorr_binned)) +
    geom_line(data = realDat, aes(x = Year, y = Broods)) + 
    geom_line(data = forecasts_avg, aes(x = Year, y = Estimate_mean, col = as.factor(propMissCat), group = propMissCat), alpha = .8) +
    theme_classic() +
    scale_color_viridis_d(option = "plasma", name = "prop. missing") +
    xlim(2005,2019)
)
#
# save figure
png(file = "./figures/forecastAccuracy_poisson.png", width = 9, height = 6, units = "in", res = 700)
forecastFig_p
dev.off()


# figure of mean predictions across all amounts missingness for each type  --------
forecasts_fig <- forecasts_avg
#forecasts_fig <- forecasts_avg %>% 
#  filter(propMiss <=0.5)

ggplot() + 
  facet_grid(.~as.factor(autocorr_binned)~ as.factor(propMissCat)) +
  #geom_ribbon(data = second_fig, aes(x = date, ymin = Estimate_mean - 1.96 * Est.Error_mean, ymax = Estimate_mean + 1.96 * Est.Error_mean, fill = type, group =type), alpha = .3) +
  geom_line(data = realDat, aes(x = Year, y = Broods)) + 
  geom_line(data = forecasts_fig, aes(x = Year, y = Estimate_mean, col = missingness, group = missingness)) +
  theme_classic() +
  xlim(2009,2019)

# figure of mean predictions across all amounts of missingness for each type  --------
forecasts_fig3 <- forecasts_long %>% 
  group_by(Year, missingness, autocorr_binned) %>% 
  summarize(Estimate_mean = mean(Estimate),
            Estimate_sd = sd(Estimate)) 

ggplot() + 
  facet_grid(.~as.factor(missingness)~ as.factor(autocorr_binned)) +
  geom_ribbon(data = forecasts_fig3, aes(x = Year, ymin = Estimate_mean - 1.96 * Estimate_sd, ymax = Estimate_mean + 1.96 * Estimate_sd, 
                                         fill = missingness, group = missingness), alpha = .3) +
  geom_line(data = realDat, aes(x = Year, y = Broods)) + 
  geom_line(data = forecasts_fig3, aes(x = Year, y= Estimate_mean, group = missingness), alpha = .3) +
  theme_classic() +
  xlim(2009,2019)





