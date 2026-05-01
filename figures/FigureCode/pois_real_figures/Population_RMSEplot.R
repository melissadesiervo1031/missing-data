
#####
# making figures for Population forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)
library(RColorBrewer)
library(ggpubr)
library(here)
library(janitor)
library(paletteer)


# read in prediction data ------------------------------------------------------------
# MCAR
# all methods except MI
# MCAR
ricker_MCAR <- readRDS("./data/model_results/pois_real_randMissRev1.rds")
ricker_MCAR_MI <- readRDS("./data/model_results/RickerRealMAR_MIRev1.rds")
ricker_MCAR <- ricker_MCAR %>% 
  left_join(ricker_MCAR_MI) %>% 
  rename(prop_miss_act = missing_result, 
         autocor_act = autocor_result2, 
         ) %>% 
  mutate(index = NA, 
         autocor = NA, 
         missingnessType = "MCAR") %>% 
  select(prop_miss_act, index, LOO_rep, drop_fits, cc_fits, DA_fits, EM_fits, MI_fits, noMiss_fits, 
         rmse_drop, rmse_cc, rmse_DA, rmse_EM, rmse_MI, rmse_noMiss, autocor, autocor_act, missingnessType)

# MNAR
# first half
ricker_MNAR <- readRDS("./data/model_results/pois_real_minMaxMiss_dropccEMDA.rds") %>% 
  mutate(autcor = NA, 
         autocor_act = NA,
         missingnessType = "MNAR") %>% 
  rename(prop_miss_act = missing_result, 
         index = missing_i)
# remove model runs for 40% missingness, which are duplicated in the next dataset
ricker_MNAR <- ricker_MNAR %>% 
  filter(prop_miss_act<0.4)

# second half 
ricker_MNAR_2 <- readRDS("./data/model_results/pois_real_minMaxMiss_dropccEMDA_secondhalf.rds") %>% 
  mutate(autcor = NA, 
         autocor_act = NA,
         missingnessType = "MNAR") %>% 
  rename(prop_miss_act = missing_result, 
         index = missing_i)

ricker_MNAR <- ricker_MNAR %>% 
  rbind(ricker_MNAR_2)

# get MI data for the second half of data 
ricker_MNAR_MI <- readRDS("./data/model_results/RickerRealMNAR_MIRev1.rds")
ricker_MNAR_MI <- ricker_MNAR_MI %>% 
  mutate(autocor = NA, 
         autocor_act = NA,
         missingnessType = "MNAR") %>% 
  rename(prop_miss_act = missing_result, 
         index = missing_i)
ricker_MNAR_3 <- ricker_MNAR %>% left_join(ricker_MNAR_MI) %>% 
  select(names(ricker_MCAR))

allDat <- rbind(ricker_MCAR, ricker_MNAR_3) #allDat <- readRDS("./data/model_results/RickerForecast_resultTableAll.rds")

# read in real data
bursaria <- read.csv("./data/Bursaria.csv")
# tidy 
bursaria <- bursaria %>% 
  janitor::clean_names() %>%
  mutate(
    date = lubridate::mdy(date),
    patch = factor(patch)
  )

# grouping by patch, then adding a number of days column
bursaria_patchlist <- purrr::map(
  unique(bursaria$patch),
  ~ {
    filter(bursaria, patch == .x) %>%
      arrange(date) %>%
      mutate(
        days = as.numeric(date - date[1], units = "days")
      )
  }
)

# now offset each so that we have a response that is time t and explanatory variable that is
# time t - 1
bursaria_patchlist_diff <- purrr::map(
  bursaria_patchlist,
  ~ {
    df2 <- tibble(
      count_tm1 = .x$number[-nrow(.x)]
    )
    cbind(.x[-1, ], df2)
  }
)

# recombine into one dataframe
bursaria_diff <- Reduce(rbind, bursaria_patchlist_diff)
# realDat$timeStep <- seq(1:59)
# train_length=49
# rename the patches into something straitforeward
patch_lu <- data.frame("patch" = unique(bursaria$patch),
                       "PatchName" = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j") %>% stringr::str_to_upper()
)
bursaria <- bursaria %>% 
  left_join(patch_lu, by = "patch")

# prep model result data 
forecasts_long <- allDat %>% 
  mutate(datasetID = (1:nrow(allDat))) %>% 
  pivot_longer(cols = c(rmse_drop:rmse_noMiss,rmse_MI, rmse_EM), names_to = "modType", values_to = "RMSE")

# allDat %>% 
#   filter(missingnessType == "MCAR" & 
#            autocor_act>=0.3 & autocor_act < 0.6 & 
#            prop_miss_act > 0.3 & prop_miss_act <=0.5 &
#          !is.na(rmse_EM)) %>% nrow()
# bin according to autocorrelation
forecasts_long$autocorr_binned <- NA
forecasts_long[forecasts_long$missingnessType == "MCAR" & forecasts_long$autocor_act <= 0.25, "autocorr_binned"] <- "low_autocorr"
forecasts_long[forecasts_long$missingnessType == "MCAR" & forecasts_long$autocor_act > 0.25 & forecasts_long$autocor_act < 0.65, "autocorr_binned"] <- "med_autocorr"
forecasts_long[forecasts_long$missingnessType == "MCAR" & forecasts_long$autocor_act >=  0.65, "autocorr_binned"] <- "high_autocorr"
forecasts_long$autocorr_binned <- factor(forecasts_long$autocorr_binned, 
                                         levels = c("low_autocorr",  "med_autocorr",  "high_autocorr"))

# bin according to missingness
forecasts_long$prop_miss_act <- as.numeric(forecasts_long$prop_miss_act)
forecasts_long$propMiss_binned <- NA
forecasts_long[!is.na(forecasts_long$prop_miss_act) & forecasts_long$prop_miss_act <= 0.3 , "propMiss_binned"] <- 0.2
forecasts_long[!is.na(forecasts_long$prop_miss_act) & forecasts_long$prop_miss_act > 0.3 & forecasts_long$prop_miss_act <= 0.5 , "propMiss_binned"] <- 0.4
forecasts_long[!is.na(forecasts_long$prop_miss_act) & forecasts_long$prop_miss_act > 0.5 , "propMiss_binned"] <- 0.6


# randomly sample the model runs so the n in the figure bins is consistent--------
## for everything aside from DA and EM, want 870 runs (for MCAR)
forecasts_long$Status <- "undefined"
#forecasts_long[forecasts_long$modType %in% c("rmse_DA", "rmse_MI"), ]$Status <- "Keep"
forecasts_long[!is.na(forecasts_long$autocorr_binned) &  forecasts_long$autocorr_binned != "med_autocorr", ]$Status <- "Discard"

for (i in seq_along(c("rmse_cc", "rmse_drop", "rmse_EM", "rmse_DA", "rmse_MI"))) {
  type_i <- c("rmse_cc", "rmse_drop", "rmse_EM", "rmse_DA", "rmse_MI")[i]
  print(type_i)
  for (j in 1:3) {
    miss_j <- c(0.2, 0.4, 0.6)[j]
    print(miss_j)
    # get 870 model runs per missingness bins for MCAR data
    # if (type_i %in% c( "rmse_DA", "rmse_MI")) {
    #   forecasts_long[forecasts_long$modType ==type_i & 
    #                    forecasts_long$propMiss_binned == miss_j & 
    #                    forecasts_long$Status == "undefined" & 
    #                    forecasts_long$missingnessType == "MCAR", "Status"] <- "Keep"
    # } else {
      temp_MCAR <- forecasts_long[forecasts_long$modType ==type_i & 
                                    forecasts_long$propMiss_binned == miss_j & 
                                    forecasts_long$Status == "undefined" & 
                                    forecasts_long$missingnessType == "MCAR",] 
      
      temp_MCAR[sample(x = 1:nrow(temp_MCAR), size = 840, replace = FALSE), "Status"] <- "Keep"
      forecasts_long[forecasts_long$modType ==type_i & 
                       forecasts_long$propMiss_binned == miss_j & 
                       forecasts_long$Status == "undefined" & 
                       forecasts_long$missingnessType == "MCAR", "Status"]  <- temp_MCAR$Status
    #}
    
    # get 96 model runs per missingness bins for MNAR data
    # if (type_i == "rmse_MI" & 
    #     miss_j == 0.6) {
      # forecasts_long[forecasts_long$modType ==type_i & 
      #                  forecasts_long$propMiss_binned == miss_j & 
      #                  forecasts_long$Status == "undefined" & 
      #                  forecasts_long$missingnessType == "MNAR", "Status"]  <-  "Keep"
   # } else {
      temp_MNAR <- forecasts_long[forecasts_long$modType ==type_i & 
                                    forecasts_long$propMiss_binned == miss_j & 
                                    forecasts_long$Status == "undefined" & 
                                    forecasts_long$missingnessType == "MNAR",] 
      
      temp_MNAR[sample(x = 1:nrow(temp_MNAR), size = 96, replace = FALSE), "Status"] <- "Keep"
      forecasts_long[forecasts_long$modType ==type_i & 
                       forecasts_long$propMiss_binned == miss_j & 
                       forecasts_long$Status == "undefined" & 
                       forecasts_long$missingnessType == "MNAR", "Status"]  <- temp_MNAR$Status
   # }
    
  }
}

forecasts_long <- forecasts_long %>% 
  filter(Status == "Keep")

# now, group rmse by bins
group_means_MCAR <- forecasts_long %>% 
  filter(missingnessType == "MCAR") %>% 
  #filter(!is.infinite(RMSE)) %>% 
  mutate("RMSE_2" = RMSE,
         "RMSE_3" = RMSE) %>% 
  group_by(modType, autocorr_binned, 
           propMiss_binned) %>% 
  summarize(RMSE = median(RMSE, na.rm = TRUE),
            quantile_low =  quantile(RMSE_2, na.rm = TRUE, probs = 0.25),
            quantile_high =  quantile(RMSE_3, na.rm = TRUE, probs = 0.75),
            RMSE_n = sum(!is.na(RMSE_2))) %>% 
  mutate(missingnessType = "MCAR") 

## for some reason there are two models w/ complete case missingness that have an RMSE of infinity? not sure why this is? (should ask Amy)
group_means_MNAR <- forecasts_long %>% 
  filter(missingnessType == "MNAR") %>% 
  #filter(is.finite(RMSE)) %>% 
  mutate("RMSE_2" = RMSE,
         "RMSE_3" = RMSE) %>% 
  group_by(modType,  propMiss_binned) %>% 
  summarize(RMSE = median(RMSE, na.rm = TRUE),
            quantile_low =  quantile(RMSE_2, na.rm = TRUE, probs = 0.25),
            quantile_high =  quantile(RMSE_3, na.rm = TRUE, probs = 0.75), 
            RMSE_n = sum(!is.na(RMSE_2))) %>% 
  mutate(missingnessType = "MNAR", 
         autocorr_binned = NA)  %>% 
  select(names(group_means_MCAR))

## add together
RMSE_df <- group_means_MCAR %>% 
  rbind(group_means_MNAR)
# get only medium autocorrelation
RMSE_df_med <- RMSE_df[RMSE_df$autocorr_binned=="med_autocorr" | RMSE_df$missingnessType == "MNAR",]
# remove RMSE for no missingness
RMSE_df_med <- RMSE_df_med %>% 
  filter(modType != "rmse_noMiss")
# add info about symbols
RMSE_df_med <- RMSE_df_med %>% 
  mutate(symbolID = "n = 840")
RMSE_df_med[ RMSE_df_med$missingnessType == "MNAR", ]$symbolID <- "n = 96"
RMSE_df_med[RMSE_df_med$RMSE_n < 840 & RMSE_df_med$missingnessType == "MCAR", ]$symbolID <- "400 < n < 840"
RMSE_df_med[RMSE_df_med$RMSE_n < 96 & RMSE_df_med$missingnessType == "MNAR", ]$symbolID <- "5 < n < 96"

RMSE_df_med$symbolID <- factor(RMSE_df_med$symbolID, levels = c("n = 840", "400 < n < 840", "n = 96", "5 < n < 96"))
  
RMSE_df_med[RMSE_df_med$missingnessType == "MCAR", "missingnessType"] <- "Missing Completely at Random"
RMSE_df_med[RMSE_df_med$missingnessType == "MNAR", "missingnessType"] <- "Missing Not at Random"

(rmse_NoLineErrorBar <- ggplot(data = RMSE_df_med) +
    facet_grid(.~missingnessType, scales = "free_y") +
    geom_linerange(aes(x = propMiss_binned, ymin = quantile_low, ymax = quantile_high, color = modType), alpha = 1, position = position_dodge(width = .1)) +
    geom_point(aes(x = propMiss_binned, y = RMSE, color = modType, shape = as.factor(symbolID)#, size =RMSE_n 
                   ), alpha = 1, position = position_dodge(width = .1)) +
   #geom_text(aes(x = propMiss_binned, y = 11.5, label = paste0("N=",RMSE_n), group = modType),  hjust = 0.2, angle =90, position = position_dodge(width = .1) ) + 
    #geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    #scale_shape_manual(values = ) 
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    scale_shape_manual(values = c(19, 1, 17, 2)) +
    scale_color_discrete(type = c(
      "#D55E00", "#CC79A7", "#E69F00",  "#8c8c8c", "#009E73"),
                         labels = c(
                           "Data Deletion-Complete","Data Augmentation", "Data Deletion-Simple", "Expectation Maximization" ,  "Multiple Imputation")
    )  +
    
    # scale_color_discrete(type = c("#CC79A7", "#E69F00", "#D55E00", "#8c8c8c", "#009E73"),
    #                      labels = c("Data Augmentation", "Data Deletion-Simple", "Data Deletion-Complete", "Expectation Maximization", "Multiple Imputation")
    # )  +
    guides(col = guide_legend(title = NULL, position = "top", direction = "vertical", nrow = 2), 
           shape = guide_legend(title = "Sample \nSize"))
)

# make plot of bursaria data time series 
(bursaria_plot <- ggplot(data = bursaria, aes(x = date, y = number, col = PatchName)) + 
  geom_line() + 
  geom_point() + 
  theme_classic() + 
  scale_color_discrete(palette = c(paletteer::paletteer_d(palette = "ggthemes::Classic_Cyclic", n = 10))) +
  ylab("Population size") + 
  xlab("Date") +
  guides(col = guide_legend(title = "Patch Name", position = "right", direction = "vertical", nrow = 5))
)
png(file = "./figures/RMSE_FullFigure_NoLineWithErrorBar_poisson_EMPIRICAL.png", width = 8.5, height = 9, units = "in", res = 700)
ggarrange(bursaria_plot,
rmse_NoLineErrorBar, ncol = 1, heights = c(.5, 1))
dev.off()

## get plot of simulated poisson RMSE
rmse_NoLineErrorBar_sim <- readRDS("./figures/RMSEfig_poissSim.rds")

ggarrange(bursaria_plot,
          rmse_NoLineErrorBar, rmse_NoLineErrorBar_sim, ncol = 1, heights = c(.5, 1, 1),
          widths = c(1,1,.8),
          align = "h",
          labels = c("A", "B", "C"))
library(patchwork)
bigPlot <- bursaria_plot/ rmse_NoLineErrorBar/ rmse_NoLineErrorBar_sim
bigPlot <- bigPlot + plot_annotation(tag_levels = c('A', '1'))

png(file = "./figures/RMSE_forAllPoissData.png", width = 8, height = 12, units = "in", res = 700)
bigPlot
dev.off()
