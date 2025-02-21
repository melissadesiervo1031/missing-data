#####
# making figures for Population forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)
library(RColorBrewer)
library(ggpubr)
library(here)


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
# can't accurately intiate the model predicting the last ten years in that
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

# remove the zero missing category
RMSE_df=RMSE_df[-which(RMSE_df$propMissCat==0),]

# Plot RMSE against missingness line plot-------------------------------------------
# facet_labels <- c(
#   "low_autocorr" = "Low Autocorrelation",
#   "med_autocorr" = "Medium Autocorrelation",
#   "high_autocorr" = "High Autocorrelation"
# )
# 
# rmse_missingness_p <- ggplot(RMSE_df) +
#   geom_point(aes(x = propMiss, y = RMSE, col = modelType), alpha = .5) +
#   geom_smooth(aes(x = propMiss, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
#   facet_wrap(~autocorr_binned, labeller = labeller(autocorr_binned = facet_labels)) +
#   scale_color_discrete(type = c("#CC79A7","#D55E00", "#E69F00", "#BBBBBB","#009E73"),
#                        labels = c("Data Augmentation","Data Deletion-Complete",  "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputations")) +
#   xlab("Proportion Missing")+
#   ylab("Root Mean Square Error (RMSE)")+
#   labs(color="Model Type")+
#   theme_classic()
# 
# 
# # save figure
# png(file = "./figures/RMSE_poisson.png", width = 8, height = 4, units = "in", res = 700)
# rmse_missingness_p
# dev.off()


# plot just the medium autocorrelation
# RMSE_df_med=RMSE_df[which(RMSE_df$autocorr_binned=="med_autocorr"),]
# rmse_missingness_med <- ggplot(RMSE_df_med) +
#   geom_point(aes(x = propMiss, y = RMSE, col = modelType), alpha = .5) +
#   geom_smooth(aes(x = propMiss, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
#   scale_color_discrete(type = c("#CC79A7","#D55E00", "#E69F00", "#BBBBBB","#009E73"),
#                        labels = c("Data Augmentation","Data Deletion-Complete",  "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputations")) +
#   xlab("Proportion Missing")+
#   ylab("Root Mean Square Error (RMSE)")+
#   labs(color="Model Type")+
#   theme_classic()
# 
# rmse_missingness_med



# Plot RMSE against missingness boxplots -------------------------------------------
# deprecated
# rmse_missingness_box <- ggplot(RMSE_df) +
#   geom_boxplot(aes(x = factor(propMissCat), y = RMSE, fill = modelType), alpha = 0.7) +
#   facet_wrap(~autocorr_binned) +
#   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   labs(
#     x = "Proportion Missingness",
#     y = "Root Mean Squared Error (RMSE)",
#     fill = "Model Type"
#   ) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 12),  # Customize facet labels
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
#   )
# 
# rmse_missingness_box
# 
# Plot RMSE against missingness intervals plus lines -------------------------------------------
# for this we may have to custom create segments to go with each method

group_means <- tapply(RMSE_df$RMSE, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), mean, na.rm = TRUE)
# start using standard deviation for error bars
# group_SD <- tapply(RMSE_df$RMSE, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), sd, na.rm = TRUE)
# group_LQ<- group_means-1.96*group_SD
# group_HQ <- group_means+1.96*group_SD

#stop using quantiles
group_LQ<- tapply(RMSE_df$RMSE, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.25)
group_HQ <- tapply(RMSE_df$RMSE, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.75)
custom_seg=data.frame(
  x=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0.0,0.012,0.024),times=9),
  xend=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0.0,0.012,0.024),times=9),
  y=c(as.vector(group_LQ[,1,]),as.vector(group_LQ[,2,]),as.vector(group_LQ[,3,])),
  yend=c(as.vector(group_HQ[,1,]),as.vector(group_HQ[,2,]),as.vector(group_HQ[,3,])),
  autocorr_binned=rep(c("low_autocorr","med_autocorr","high_autocorr"),each=nrow(group_LQ[,1,])*ncol(group_LQ[,1,])),
  modelType=rep(rownames(group_LQ[,1,]),9),
  means1=c(as.vector(group_means[,1,]),as.vector(group_means[,2,]),as.vector(group_means[,3,]))
)


RMSE_df$autocorr_binned <- factor(RMSE_df$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
custom_seg$autocorr_binned <- factor(custom_seg$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))


rmse_missingness_int <- ggplot(RMSE_df) +
  geom_smooth(aes(x = propMiss, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
  geom_segment(data=custom_seg,aes(x=x,y=y,xend=xend,yend=yend,col=modelType),size=0.6)+
  geom_point(data=custom_seg,aes(x=x,y=means1,col=modelType),size=1)+
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

rmse_missingness_int


# get only medium autocorrelation
RMSE_df_med=RMSE_df[which(RMSE_df$autocorr_binned=="med_autocorr"),]
custom_seg=custom_seg[which(custom_seg$autocorr_binned=="med_autocorr"),]

# right factor order
RMSE_df$modelType=factor(RMSE_df$modelType,levels=c("dropNA","dropCC","MI","EM","DA"))
custom_seg$modelType=factor(custom_seg$modelType,levels=c("dropNA","dropCC","MI","EM","DA"))

rmse_missingness_int_med <- ggplot(RMSE_df_med) +
  #geom_smooth(aes(x = propMiss, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
  geom_segment(data=custom_seg,aes(x=x,y=y,xend=xend,yend=yend,col=modelType),size=0.6)+
  geom_point(data=custom_seg,aes(x=x,y=means1,col=modelType),size=1)+
  scale_x_continuous(breaks = c(0.2,0.4,0.6)) +
  scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73", "#BBBBBB","#CC79A7"),
                       labels = c("Data Deletion-Simple", "Data Deletion-Complete", "Multiple Imputations", "Expectation Maximization", "Data Augmentation")) +
  labs(
    x = "Proportion Missing",
    y = "Root Mean Square Error (RMSE)",
    color = "Model Type"
  ) +
  theme_classic() +
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  guides(color = guide_legend(nrow = 2))+
  theme(
    strip.text = element_text(size = 12),  # Customize facet labels
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
  )

rmse_missingness_int_med

# add the wytham tits to the plot
wytham <- read.csv(here("data/Wytham_tits.csv"), header=TRUE)

wytham_tits = ggplot() + 
  geom_rect(aes(xmin = 2008.5, xmax = 2018.5, 
                ymin = 100, ymax = 500), fill = "grey80") +
  geom_line(data = wytham, aes(x = Year, y = Broods)) + 
  geom_point(data = wytham, aes(x = Year, y = Broods),fill="gray",color="black",size=3,shape=21) + 
  theme_classic() + 
  labs(x = "Year", y = "Number of Broods")
wytham_tits

png(file = "./figures/RMSE_pois_combined.png", width = 6, height = 7, units = "in", res = 700)
ggarrange(wytham_tits,rmse_missingness_int_med,
          labels = c("A", "B"),
          ncol=1,nrow=2,
          heights = c(1, 2))
dev.off()

# Make figure of CI coverage -------------------------------------------
# deprecated for now...
# Error/confidence intervals should expand over time like for simulations with temporal autocorrelation
# This is because we are basing each next data point prediction on the last data point prediction
# We will therefore need to re-calculate the CI- I think we would do it in the modelRuns_rickerForecasts.R script
# Use the lower and upper estimates to expand a lower and upper CI over time
# We didn't need to do this with the simulated data because it was just measuring coverage of a simple CI
# Will need to add extra prediction columns for each

# make coverage into proportion
# RMSE_df$coverage=RMSE_df$coverage/10
# 
# # simple scatterplot- kind of ugly honestly ----------------------------------------------
# coverage_plot <- ggplot(RMSE_df) +
#   geom_point(aes(x = propMiss, y = coverage, col = modelType), alpha = .5) +
#   geom_smooth(aes(x = propMiss, y = coverage, col = modelType), method = "lm", se = FALSE) + 
#   facet_wrap(~autocorr_binned, ) +
#   scale_color_discrete(type = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                        labels = c("Data Aug.","Data Del.-Complete",  "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   theme_classic()
# 
# coverage_plot

# simple boxplot- even more ugly honestly ----------------------------------------------
# coverage_box <- ggplot(RMSE_df) +
#   geom_boxplot(aes(x = factor(propMissCat), y = coverage, fill = modelType), alpha = 0.7) +
#   facet_wrap(~autocorr_binned) +
#   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   labs(
#     x = "Proportion Missingness",
#     y = "Coverage",
#     fill = "Model Type"
#   ) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 12),  # Customize facet labels
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
#   )
# 
# coverage_box

# Plot coverage against missingness intervals plus lines -------------------------------------------
# for this we may have to custom create segments to go with each method
# group_means <- tapply(RMSE_df$coverage, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), mean, na.rm = TRUE)
# group_LQ<- tapply(RMSE_df$coverage, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.25)
# group_HQ <- tapply(RMSE_df$coverage, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.75)
# custom_seg=data.frame(
#   x=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
#   xend=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
#   y=c(as.vector(group_LQ[,1,]),as.vector(group_LQ[,2,]),as.vector(group_LQ[,3,])),
#   yend=c(as.vector(group_HQ[,1,]),as.vector(group_HQ[,2,]),as.vector(group_HQ[,3,])),
#   autocorr_binned=rep(c("low_autocorr","med_autocorr","high_autocorr"),each=nrow(group_LQ[,1,])*ncol(group_LQ[,1,])),
#   modelType=rep(rownames(group_LQ[,1,]),12),
#   means1=c(as.vector(group_means[,1,]),as.vector(group_means[,2,]),as.vector(group_means[,3,]))
# )
# 
# 
# RMSE_df$autocorr_binned <- factor(RMSE_df$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# custom_seg$autocorr_binned <- factor(custom_seg$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# 
# 
# coverage_int <- ggplot(RMSE_df) +
#   geom_smooth(aes(x = propMiss, y = coverage, col = modelType), method = "lm", se = FALSE) + 
#   geom_segment(data=custom_seg,aes(x=x,y=y,xend=xend,yend=yend,col=modelType),size=0.6)+
#   geom_point(data=custom_seg,aes(x=x,y=means1,col=modelType),size=1)+
#   facet_wrap(~autocorr_binned) +
#   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   labs(
#     x = "Proportion Missingness",
#     y = "Coverage",
#     fill = "Model Type"
#   ) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 12),  # Customize facet labels
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
#   )
# 
# coverage_int


# Make figure of ALTERNATE CI coverage -------------------------------------------
# deprecated for now
# Error/confidence intervals should expand over time like for simulations with temporal autocorrelation
# This is because we are basing each next data point prediction on the last data point prediction
# We will therefore need to re-calculate the CI- I think we would do it in the modelRuns_rickerForecasts.R script
# Use the lower and upper estimates to expand a lower and upper CI over time
# We didn't need to do this with the simulated data because it was just measuring coverage of a simple CI
# Will need to add extra prediction columns for each

# make coverage into proportion
# RMSE_df$coverage_opt2=RMSE_df$coverage_opt2/10
# 
# # simple scatterplot- kind of ugly honestly ----------------------------------------------
# coverage_opt2_plot <- ggplot(RMSE_df) +
#   geom_point(aes(x = propMiss, y = coverage_opt2, col = modelType), alpha = .5) +
#   geom_smooth(aes(x = propMiss, y = coverage_opt2, col = modelType), method = "lm", se = FALSE) + 
#   facet_wrap(~autocorr_binned, ) +
#   scale_color_discrete(type = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                        labels = c("Data Aug.","Data Del.-Complete",  "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   theme_classic()
# 
# coverage_opt2_plot
# 
# # simple boxplot- even more ugly honestly ----------------------------------------------
# coverage_opt2_box <- ggplot(RMSE_df) +
#   geom_boxplot(aes(x = factor(propMissCat), y = coverage_opt2, fill = modelType), alpha = 0.7) +
#   facet_wrap(~autocorr_binned) +
#   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   labs(
#     x = "Proportion Missingness",
#     y = "Coverage",
#     fill = "Model Type"
#   ) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 12),  # Customize facet labels
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
#   )
# 
# coverage_opt2_box
# 
# # Plot coverage_opt2 against missingness intervals plus lines -------------------------------------------
# # for this we may have to custom create segments to go with each method
# group_means <- tapply(RMSE_df$coverage_opt2, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), mean, na.rm = TRUE)
# group_LQ<- tapply(RMSE_df$coverage_opt2, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.25)
# group_HQ <- tapply(RMSE_df$coverage_opt2, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.75)
# custom_seg=data.frame(
#   x=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
#   xend=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
#   y=c(as.vector(group_LQ[,1,]),as.vector(group_LQ[,2,]),as.vector(group_LQ[,3,])),
#   yend=c(as.vector(group_HQ[,1,]),as.vector(group_HQ[,2,]),as.vector(group_HQ[,3,])),
#   autocorr_binned=rep(c("low_autocorr","med_autocorr","high_autocorr"),each=nrow(group_LQ[,1,])*ncol(group_LQ[,1,])),
#   modelType=rep(rownames(group_LQ[,1,]),12),
#   means1=c(as.vector(group_means[,1,]),as.vector(group_means[,2,]),as.vector(group_means[,3,]))
# )
# 
# 
# RMSE_df$autocorr_binned <- factor(RMSE_df$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# custom_seg$autocorr_binned <- factor(custom_seg$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# 
# 
# coverage_opt2_int <- ggplot(RMSE_df) +
#   geom_smooth(aes(x = propMiss, y = coverage_opt2, col = modelType), method = "lm", se = FALSE) + 
#   geom_segment(data=custom_seg,aes(x=x,y=y,xend=xend,yend=yend,col=modelType),size=0.6)+
#   geom_point(data=custom_seg,aes(x=x,y=means1,col=modelType),size=1)+
#   facet_wrap(~autocorr_binned) +
#   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   labs(
#     x = "Proportion Missingness",
#     y = "Coverage",
#     fill = "Model Type"
#   ) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 12),  # Customize facet labels
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
#   )
# 
# coverage_opt2_int



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





