#####
# making figures for GPP forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)

# read in prediction data ------------------------------------------------------------
allDat <- readRDS("./data/model_results/RickerForecast_resultTableAll.rds")
allDat[is.nan(allDat$actAutoCorr_trim), "actAutoCorr_trim"] <- 0

# put predictions into a 'long' format
forecasts_long <- list_rbind(apply(allDat, MARGIN = 1, FUN = function(x) {
  startInd <- which(!is.na(x$forecasts$timeStep))[2]
 data.frame(
 "newName" = x$newName, "autoCorr" = x$actAutoCorr_trim, "propMiss" = x$actPropMiss_trim,
 "forecast_dropNA" = x$forecasts$dropNA_est[startInd:59], 
 "forecast_dropCC" = x$forecasts$dropCC_est[startInd:59], 
 "forecast_MI" = x$forecasts$MI_est[startInd:59], 
 "forecast_EM" = x$forecasts$EM_est[startInd:59], 
 "forecast_DA" = x$forecasts$DA_est[startInd:59]
 )
}))
# bin according to autocorrelation
forecasts_long$autocorr_binned <- NA
forecasts_long[forecasts_long$autoCorr < 0.3, "autocorr_binned"] <- "low_autocorr"
forecasts_long[forecasts_long$autoCorr >= 0.3 & forecasts_long$autoCorr < 0.6, "autocorr_binned"] <- "med_autocorr"
forecasts_long[forecasts_long$autoCorr >=  0.6, "autocorr_binned"] <- "high_autocorr"
forecasts_long$autocorr_binned <- factor(forecasts_long$autocorr_binned, 
                                   levels = c("low_autocorr",  "med_autocorr",  "high_autocorr"))


# put RMSE data into it's own data frame
RMSE_df <- list_rbind(apply(allDat, MARGIN = 1, FUN = function(x) {
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

# Plot RMSE against missingness -------------------------------------------
ggplot(RMSE_df) +
  geom_point(aes(x = propMiss, y = RMSE, col = modelType)) +
  geom_smooth(aes(x = propMiss, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
  facet_wrap(~autocorr_binned, ) +
  theme_classic()

