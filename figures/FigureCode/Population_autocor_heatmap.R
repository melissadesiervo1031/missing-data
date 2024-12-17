## draft of the Population heat map figure ##

# Load packages
library(tidyverse)
library(ggpubr)

## read in data 
# read in data and prepare for figures------------------------------------------

ricDat_tempA <- readRDS("./data/model_results/RickerA_resultTableAll.rds")
ricDat_tempB <- readRDS("./data/model_results/RickerB_resultTableAll.rds")
ricDat_temp <- rbind(ricDat_tempA, ricDat_tempB)



## remove DA and MI runs from previous data, since these have been rerun
ricDat_tempNew <- ricDat_temp %>% 
  select(-DA_fits, -MI_fits)

## load re-runs of DA and MI models 
# load DA re-runs
fileNames_DA <- c(paste0("./data/model_results/ricker_Sim_reruns/DA_reruns/RickerA_resultTable",1:30,".rds"),
                  paste0("./data/model_results/ricker_Sim_reruns/DA_reruns/RickerB_resultTable",1:30,".rds") 
)


DA_temp <- map_df(fileNames_DA, function(x)
  readRDS(x)
)

# add DA to other ricker result data
ricDat_tempNew <- ricDat_tempNew %>% 
  left_join(DA_temp)



# add DA to other ricker result data
ricDat_tempNew <- ricDat_tempNew %>% 
  left_join(DA_temp)

## will add MI later (Since it's already in the format we ultimately want

# for now, remove the simulation 176 (really tiny simulation parameters)
ricDat_tempNew <- ricDat_tempNew[ricDat_tempNew$SimNumber != 176,]
## for now, remove rows where simulated population went extinct (was 13% of data!)
# badRows <- as.vector(sapply(ricDat_tempNew$drop_fits, function(x)
#   ifelse(sum(str_detect(names(x), pattern = "cause"))>=1, 
#              yes = "bad", 
#              no = "good")
# ))
# ricDat_tempNew <- ricDat_tempNew[rownames(ricDat_tempNew) %in% which(badRows == "good"),]
# rename columns 
ricDat_tempNew <- ricDat_tempNew %>% 
  rename(r_sim = r, alpha_sim = alpha, N0_sim = N0)
# put input autoCor and propMiss data into one column
ricDat_tempNew <- ricDat_tempNew %>% 
  mutate(input_args = paste0("a=", actAutoCorr, "_", "p=", actPropMiss)) %>% 
  select(-autoCorr, -propMiss) %>% 
  rename("autoCorr" = "actAutoCorr", "propMiss" = "actPropMiss")

# extract parameter information from the list columns
ricDat <- rbind(
  #drop na fits
  cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
                          "N0_sim","input_args", "autoCorr", "propMiss")], 
        map_df(ricDat_tempNew$drop_fits, function(x) 
          ##
          if (length(names(x)) < 3) {
            data.frame(
              "type" = "dropNA",
              "r_est" =NA, 
              "alpha_est" = NA,
              "r_se" = NA,
              "alpha_se" = NA,
              "status" = "missingnessLimitReached")
          } else {
            data.frame(
              "type" = "dropNA",
              "r_est" = x$estim["r"], 
              "alpha_est" = x$estim["alpha"],
              "r_se" = x$se["r"],
              "alpha_se" = x$se["alpha"],
              "status" = "good")
          }
        ),
        data.frame("listName" = names(ricDat_tempNew$drop_fits))
  ),
  #dropNA complete case fits
  cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
                          "N0_sim","input_args", "autoCorr", "propMiss")], 
        map_df(ricDat_tempNew$cc_fits, function(x) 
          if (length(names(x)) < 3) {
            data.frame(
              "type" = NA,
              "r_est" =NA, 
              "alpha_est" = NA,
              "r_se" = NA,
              "alpha_se" = NA,
              "status" = "missingnessLimitReached")
          } else {
            data.frame(
              "type" = "CompleteCaseDropNA",
              "r_est" = x$estim["r"], 
              "alpha_est" = x$estim["alpha"],
              "r_se" = x$se["r"],
              "alpha_se" = x$se["alpha"],
              "status" = "good")
          }
        ),
        data.frame("listName" = names(ricDat_tempNew$drop_fits))
  ),
  #EM fits
  cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
                          "N0_sim","input_args", "autoCorr", "propMiss")], 
        map_df(ricDat_tempNew$EM_fits, function(x) 
          if (length(names(x)) < 3) {
            data.frame(
              "type" = NA,
              "r_est" =NA, 
              "alpha_est" = NA,
              "r_se" = NA,
              "alpha_se" = NA,
              "status" = "missingnessLimitReached")
          } else {
            data.frame(
              "type" = "ExpectationMaximization",
              "r_est" = x$estim["r"], 
              "alpha_est" = x$estim["alpha"],
              "r_se" = x$se["r"],
              "alpha_se" = x$se["alpha"],
              "status" = "good")
          }
        ),
        data.frame("listName" = names(ricDat_tempNew$drop_fits))
        
  ),
  #DA fits
  cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
                          "N0_sim","input_args", "autoCorr", "propMiss")], 
        map_df(ricDat_tempNew$DA_fits, function(x) 
          if (length(names(x)) < 3) {
            data.frame(
              "type" = NA,
              "r_est" =NA, 
              "alpha_est" = NA,
              "r_se" = NA,
              "alpha_se" = NA,
              "status" = "missingnessLimitReached")
          } else {
            data.frame(
              "type" = "DataAugmentation",
              "r_est" = x$estim["r"], 
              "alpha_est" = x$estim["alpha"],
              "r_se" = x$se["r"],
              "alpha_se" = x$se["alpha"],
              "status" = "good")
          }
        ),
        data.frame("listName" = names(ricDat_tempNew$drop_fits))
        
  )
  #,
  # #MI fits
  # cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
  #                      "N0_sim","input_args", "autoCorr", "propMiss")], 
  #       map_df(ricDat_tempNew$MI_fits, function(x) 
  #         if (length(names(x)) < 3) {
  #           data.frame(
  #             "type" = NA,
  #             "r_est" =NA, 
  #             "alpha_est" = NA,
  #             "r_se" = NA,
  #             "alpha_se" = NA,
  #             "status" = "missingnessLimitReached")
  #         } else {
  #           data.frame(
  #             "type" = "MultipleImputations",
  #             "r_est" = x$estim["r"], 
  #             "alpha_est" = x$estim["alpha"],
  #             "r_se" = x$se["r"],
  #             "alpha_se" = x$se["alpha"],
  #             "status" = "good")
  #         }
  #       ),
  #       data.frame("listName" = names(ricDat_tempNew$drop_fits))
  #       
  # )
)

# remove NAs (from missingness Limit Reached issue)
ricDat_new <- ricDat[ricDat$status != "missingnessLimitReached",]
ricDat_new <- unique(ricDat_new)
## note, "autoCorr" and "propMiss" in ricDat are actual values, not input values
# get autocorrelation amount 
#ricDat_new$actAutoCorr <- as.numeric(str_split(ricDat_new$listName, pattern = "_", simplify = TRUE)[,9])

## add MI data back in 
# load MI re-runs
fileNames_MI <- c(paste0("./data/model_results/ricker_Sim_reruns/MI_reruns/results_rickerA_",1:10,".csv"),
                  paste0("./data/model_results/ricker_Sim_reruns/MI_reruns/results_rickerB_",1:10,".csv") 
)

MI_temp <- map_df(fileNames_MI, function(x) {
  temp <- read.csv(x) %>% 
    select(-X.1, -X)
  return(temp)
}
)
# ## load MI extinction data (remove these values from the model results)
# extinctMI <- read.csv("./data/model_results/ricker_Sim_reruns/MI_reruns/extinctMI_ALL.csv") %>% 
#   select(-X.1, -X) %>% 
#   filter(!is.na(simNumber)) %>% 
#   rename(SimNumber = simNumber,
#          actAutoCorr = autocorr_act, 
#          autoCorr = autocorrs, 
#          actPropMiss = prop_miss) %>% 
#   mutate(extinct = TRUE)
# 
# MI_tempTemp <- MI_temp %>% 
#   left_join(extinctMI %>% 
#               select(SimNumber, actPropMiss, actAutoCorr, r, alpha, N0, extinct))


# remove runs where population went extinct 
good <- !str_locate(MI_temp$estim_r, pattern = " ")[,1]
good[is.na(good)] <- TRUE

MI_temp <- MI_temp[good,]

test <- MI_temp %>% 
  mutate(input_args = paste0("a=",autoCorr,"_p=",propMiss),
         type = "MultipleImputations",
         status = NA, 
         listName = NA,
         estim_r = as.numeric(estim_r)
  ) %>% 
  select(-id, -autoCorr, -propMiss) %>% 
  rename(r_sim = r, 
         alpha_sim = alpha,
         N0_sim = N0, 
         autoCorr = actAutoCorr,
         propMiss = actPropMiss,
         r_est = estim_r, 
         alpha_est = estim_alpha,
         r_se = se_r,
         alpha_se = se_alpha
  ) %>% 
  select(SimNumber, r_sim, alpha_sim, N0_sim, input_args, autoCorr, 
         propMiss, type, r_est, alpha_est, r_se, alpha_se, status, listName) %>% 
  # remove runs that have an NA for parameter estimates
  filter(!is.na(SimNumber)) %>% 
  unique()


ricDat_new <- ricDat_new %>% 
  mutate(r_est = as.numeric(r_est)) %>% 
  rbind(test)

# actual zeros in autocorrelation are listed as NAs--fix this 
ricDat_new[is.na(ricDat_new$autoCorr),"autoCorr"] <- 0

#make into long data.frame 
paramEstLong <- ricDat_new %>% 
  #select(SimNumber, r_sim, actAutoCorr, actPropMiss, type, r_est, r_se, status, r_paramDiff) %>% 
  pivot_longer(cols = c(r_est, alpha_est), 
               values_to = "paramEst", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  select(-r_sim, -alpha_sim, -N0_sim, -r_se, -alpha_se)

paramSimLong <- ricDat_new %>% 
  pivot_longer(cols = c(r_sim, alpha_sim), 
               values_to = "paramSim", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  select(-r_est, -alpha_est, -N0_sim, -r_se, -alpha_se) %>% 
  unique()

paramSELong <- ricDat_new %>% 
  pivot_longer(cols = c(r_se, alpha_se), 
               values_to = "paramSE", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  select(-r_sim, -alpha_sim, -N0_sim, -r_est, -alpha_est) %>% 
  # remove the values for Expectation Maximization, since we don't have SE for that method
  filter(type != "ExpectationMaximization")


ricDat_new_long <- left_join(paramEstLong, paramSimLong) %>% 
  left_join(paramSELong)

#calculate standardized parameter estimates
ricDat_new_long <- ricDat_new_long %>%
  mutate("paramDiff" = (paramEst - paramSim)/abs(paramSim),
         "paramDiff_abs" = abs(paramEst - paramSim)/abs(paramSim))

# filter for low and high autocor
ricDat_new_long[ricDat_new_long$autoCorr <=0.3 & !is.na(ricDat_new_long$autoCorr), "missingness"] <- "MAR: Low AC"
ricDat_new_long[ricDat_new_long$autoCorr >0.3 & ricDat_new_long$autoCorr <0.6 & !is.na(ricDat_new_long$autoCorr), "missingness"] <- "MAR: Med. AC"
ricDat_new_long[ricDat_new_long$autoCorr  >= 0.6 & !is.na(ricDat_new_long$autoCorr), "missingness"] <- "MAR: High AC"


ricDat_new_long$type <- factor(ricDat_new_long$type, levels = c("DataAugmentation", "CompleteCaseDropNA", 
                                                                "dropNA", "ExpectationMaximization", "MultipleImputations"))

ric_sim_figDat <- ricDat_new_long
##make heat map! ##
# Make heatmaps for Gaussian MAR data -------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- ric_sim_figDat %>%
  mutate(autoCorr = round(autoCorr, 1), 
         propMiss = round(propMiss, 1)) %>% 
  group_by(missingness, type, param, autoCorr, propMiss) %>% 
  summarize(paramDiffAbsDiff_mean = mean(paramDiff_abs, na.rm = TRUE),
            paramDiffAbsDiff_med = median(paramDiff_abs, na.rm = TRUE),
            paramDiffAbsDiff_SD = sd(paramDiff_abs, na.rm = TRUE),
            n_paramDiffAbsDiff = length(paramDiff_abs),
            paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  #filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
  filter(propMiss <=.5)
# only consider missingness of 50% or less

## make heatmap for mean of parameter recovery
(heatMap_median_MAR <-ggplot(data = figDat, aes(x=propMiss, y=autoCorr)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    scale_fill_viridis_c(name = "value" , option = "A") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("alpha", "r")) ~ type) +
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    guides(fill = guide_colorbar("Bias")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Parameter bias across simulations")) 

## make heatmap for mean of *absolute value* of parameter recovery
(heatMap_SE_MAR <-ggplot(data = figDat, aes(x=propMiss, y=autoCorr)) + 
    geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
    scale_fill_viridis_c(name = "value" , option = "D") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("alpha", "r")) ~ type) +
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    guides(fill = guide_colorbar("Stand. Error")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Absolute Standard Error of Parameter Recover")) 

## make heatmap for coverage for parameter recovery
figDat_cov_temp <- ricDat_new_long %>% 
  mutate(CI95_lower = paramEst - 1.96*paramSE,
         CI95_upper = paramEst + 1.96*paramSE) %>% 
  filter(param != "sigma") %>% 
  filter(!is.na(paramSE))# randomly there are some model runs that don't have SE? 

# is the true parameter within the 95% CI? 
figDat_cov_temp$coverage <- c(figDat_cov_temp$paramSim >= figDat_cov_temp$CI95_lower & 
                                figDat_cov_temp$paramSim <= figDat_cov_temp$CI95_upper)

## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
figDat_cov <- figDat_cov_temp %>% 
  filter(param != "sigma",
         param != "intercept",
         propMiss <=.5,) %>% 
  mutate(autoCorr = round(autoCorr, 1), 
         propMiss = round(propMiss, 1)
  ) %>% 
  group_by(missingness, type, param, autoCorr, propMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)

(heatMap_cov_MAR <- ggplot(data = figDat_cov, aes(x=propMiss, y=autoCorr)) + 
    geom_tile(aes(fill=coveragePerc*100), size=5) + 
    facet_grid(~factor(figDat_cov$param, levels = c("alpha", "r")) ~ type) +
    #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +
    scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    guides(fill = guide_colorbar("% Coverage")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("% of model runs where the 95% CI \n includes the simulation parameter"))

## save figures
png(file = "./figures/heatmap_PoissonMAR_median.png", width = 7, height = 6, units = "in", res = 700)
heatMap_median_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMAR_SE.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SE_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMAR_coverage.png", width = 7, height = 6, units = "in", res = 700)
heatMap_cov_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_median_MAR, heatMap_SE_MAR, heatMap_cov_MAR)
dev.off()
