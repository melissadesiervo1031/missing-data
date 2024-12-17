#/////////////////
# This script makes figures for Population simulated data 
# 4 December 2023
#/////////////////

# Load packages
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

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

# reformat data 
figDat_lines <- ricDat_new_long %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  
  summarize(paramDiffAbsDiff_mean = mean(paramDiff_abs, na.rm = TRUE),
            paramDiffAbsDiff_med = median(paramDiff_abs, na.rm = TRUE),
            paramDiffAbsDiff_SD = sd(paramDiff_abs, na.rm = TRUE),
            n_paramDiffAbsDiff = length(paramDiff_abs),
            paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n_paramDiff = length(paramDiff),
            SE_mean = mean(paramSE, na.rm = TRUE) # the mean of the parameter standard error (not standardized, but maybe should be?)
  ) %>% 
  
  #filter(n  > 50)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.5)
  # make types in the 'correct' order
 

# parameter recovery bias -------------------------------------------------
(poiss_paramRecovery_bias_MAR <- ggplot(data = figDat_lines, aes(x = amtMiss, y = paramDiff_med)) +
   facet_grid(~factor(param, levels = c("r","alpha")) 
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")), 
              scales = "free_y") + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
   #size=0.3, width=0, position = position_dodge(width=0.03))+
   #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
   geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top")+
   theme(legend.title=element_blank())+
   ylab("Median of parameter bias across sims.")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.55)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A", "#E6AB02","#7570B3"),
                       labels = c("Data Augmentation", "Data Deletion-Complete", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))#+
)

# parameter recovery SE ---------------------------------------------------
(poiss_paramRecovery_SE_MAR <- ggplot(data = figDat_lines, aes(x = amtMiss, y = paramDiffAbsDiff_med)) +
   facet_grid(~factor(param, levels = c( "r","alpha")) 
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")), 
              scales = "free_y") + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
   #size=0.3, width=0, position = position_dodge(width=0.03))+
   #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
   geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top")+
   theme(legend.title=element_blank())+
   ylab("Absolute Standard Error of Parameter Recovery")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.55)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Augmentation", "Data Deletion-Complete", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))#+
)
# parameter recovery coverage ---------------------------------------------
# does the confidence interval contain the true parameter? 
# calculate 95% CI for each param
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
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)
  ) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)


(poiss_paramRecovery_coverage_MAR <- ggplot(data = figDat_cov, aes(x = amtMiss, y = coveragePerc)) +
    facet_grid(~factor(param, levels = c( "r","alpha")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR")), 
               scales = "free_y") + 
    #geom_col(aes(x = amtMiss, y = coveragePerc, color = as.factor(type), fill = as.factor(type)), 
    # position = "dodge", alpha = .5) + 
    geom_hline(aes(yintercept = .95), colour = "grey") +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("% of model runs where the 95% CI \n includes the simulation parameter")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    xlim(c(-0.03,0.55)) + 
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#7570B3"),
                       labels = c("Data Augmentation", "Data Deletion-Complete", "Data Deletion-Simple", "Multiple Imputation"))#+
  
    #ylim(c(0,.3)) +
   # scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A","#7570B3"),
 #                        labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple",  "Multiple Imputation"))#+
)

# put figures together ----------------------------------------------------
(poiss_paramRecovMAR <- ggarrange(poiss_paramRecovery_bias_MAR, poiss_paramRecovery_SE_MAR, 
                                  poiss_paramRecovery_coverage_MAR, common.legend = TRUE, nrow = 1))

## save results
png(file = "./figures/parameterRecoveryPoisson_MAR.png", width = 12, height = 4, units = "in", res = 700)
poiss_paramRecovMAR
dev.off()

# ricDat_new_long <- ricDat_new_long %>% 
#   select(-actAutoCorr)
# # save data to file for use later...
# write_rds(ricDat_new_long, file = "./data/model_results/ricker_sim_ModelResultsLong.rds")

# # add in extinct ts data
# ricDat_new_ext <- readRDS("./data/model_results/RickerExtinct_resultTableAll.rds")
# ricDat_new_ext <- ricDat_new_ext %>% 
#   rename(r_sim = r, alpha_sim = alpha, N0_sim = N0)
# # put input autoCor and propMiss data into one column
# ricDat_new_ext <- ricDat_new_ext %>% 
#   mutate(input_args = paste0("a=", actAutoCorr, "_", "p=", actPropMiss))
# 
# # extract parameter information from the list columns
# ricDat_new_extinctAll <- rbind(
#   #drop na fits
#   cbind(ricDat_new_ext[,c("simNumber", "r_sim", "alpha_sim", 
#                        "N0_sim","input_args", "actAutoCorr", "actPropMiss")], 
#         map_df(ricDat_new_ext$drop_fits, function(x) 
#           ##
#           if (length(names(x)) < 3) {
#             data.frame(
#               "type" = "dropNA",
#               "r_est" =NA, 
#               "alpha_est" = NA,
#               "r_se" = NA,
#               "alpha_se" = NA,
#               "status" = "missingnessLimitReached")
#           } else {
#             data.frame(
#               "type" = "dropNA",
#               "r_est" = x$estim["r"], 
#               "alpha_est" = x$estim["alpha"],
#               "r_se" = x$se["r"],
#               "alpha_se" = x$se["alpha"],
#               "status" = "good")
#           }
#         ),
#         data.frame("listName" = names(ricDat_new_ext$drop_fits))
#   ),
#   #dropNA complete case fits
#   cbind(ricDat_new_ext[,c("simNumber", "r_sim", "alpha_sim", 
#                        "N0_sim","input_args", "actAutoCorr", "actPropMiss")], 
#         map_df(ricDat_new_ext$cc_fits, function(x) 
#           if (length(names(x)) < 3) {
#             data.frame(
#               "type" = NA,
#               "r_est" =NA, 
#               "alpha_est" = NA,
#               "r_se" = NA,
#               "alpha_se" = NA,
#               "status" = "missingnessLimitReached")
#           } else {
#             data.frame(
#               "type" = "CompleteCaseDropNA",
#               "r_est" = x$estim["r"], 
#               "alpha_est" = x$estim["alpha"],
#               "r_se" = x$se["r"],
#               "alpha_se" = x$se["alpha"],
#               "status" = "good")
#           }
#         ),
#         data.frame("listName" = names(ricDat_new_ext$drop_fits))
#   ),
#   #EM fits
#   cbind(ricDat_new_ext[,c("simNumber", "r_sim", "alpha_sim", 
#                        "N0_sim","input_args", "actAutoCorr", "actPropMiss")], 
#         map_df(ricDat_new_ext$EM_fits, function(x) 
#           if (length(names(x)) < 3) {
#             data.frame(
#               "type" = NA,
#               "r_est" =NA, 
#               "alpha_est" = NA,
#               "r_se" = NA,
#               "alpha_se" = NA,
#               "status" = "missingnessLimitReached")
#           } else {
#             data.frame(
#               "type" = "ExpectationMaximization",
#               "r_est" = x$estim["r"], 
#               "alpha_est" = x$estim["alpha"],
#               "r_se" = x$se["r"],
#               "alpha_se" = x$se["alpha"],
#               "status" = "good")
#           }
#         ),
#         data.frame("listName" = names(ricDat_new_ext$drop_fits))
#         
#   ),
#   #DA fits
#   # cbind(ricDat_new_ext[,c("simNumber", "r_sim", "alpha_sim", 
#   #                      "N0_sim","input_args", "actAutoCorr", "actPropMiss")], 
#   #       map_df(ricDat_new_ext$DA_fits, function(x) 
#   #         if (length(names(x)) < 3) {
#   #           data.frame(
#   #             "type" = NA,
#   #             "r_est" =NA, 
#   #             "alpha_est" = NA,
#   #             "r_se" = NA,
#   #             "alpha_se" = NA,
#   #             "status" = "missingnessLimitReached")
#   #         } else {
#   #           data.frame(
#   #             "type" = "DataAugmentation",
#   #             "r_est" = x$estim["r"], 
#   #             "alpha_est" = x$estim["alpha"],
#   #             "r_se" = x$se["r"],
#   #             "alpha_se" = x$se["alpha"],
#   #             "status" = "good")
#   #         }
#   #       ),
#   #       data.frame("listName" = names(ricDat_new_ext$drop_fits))
#   #       
#   # ),
#   #MI fits
#   cbind(ricDat_new_ext[,c("simNumber", "r_sim", "alpha_sim", 
#                        "N0_sim","input_args", "actAutoCorr", "actPropMiss")], 
#         map_df(ricDat_new_ext$MI_fits, function(x) 
#           if (length(names(x)) < 3) {
#             data.frame(
#               "type" = NA,
#               "r_est" =NA, 
#               "alpha_est" = NA,
#               "r_se" = NA,
#               "alpha_se" = NA,
#               "status" = "missingnessLimitReached")
#           } else {
#             data.frame(
#               "type" = "MultipleImputations",
#               "r_est" = x$estim["r"], 
#               "alpha_est" = x$estim["alpha"],
#               "r_se" = x$se["r"],
#               "alpha_se" = x$se["alpha"],
#               "status" = "good")
#           }
#         ),
#         data.frame("listName" = names(ricDat_new_ext$drop_fits))
#         
#   )
# )
# 
# # remove NAs (from missingness Limit Reached issue)
# ricDat_new_extinctAll <- ricDat_new_extinctAll[ricDat_new_extinctAll$status != "missingnessLimitReached",]
# 
# #make into long data.frame 
# paramEstLong_ext <- ricDat_new_extinctAll %>% 
#   #select(SimNumber, r_sim, actAutoCorr, actPropMiss, type, r_est, r_se, status, r_paramDiff) %>% 
#   pivot_longer(cols = c(r_est, alpha_est), 
#                values_to = "paramEst", 
#                names_to = "param", 
#                names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
#   select(-r_sim, -alpha_sim, -N0_sim, -r_se, -alpha_se)
# paramSimLong_ext <- ricDat_new_extinctAll %>% 
#   pivot_longer(cols = c(r_sim, alpha_sim), 
#                values_to = "paramSim", 
#                names_to = "param", 
#                names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
#   select(-r_est, -alpha_est, -N0_sim, -r_se, -alpha_se)
# 
# paramSELong_ext <- ricDat_new_extinctAll %>% 
#   pivot_longer(cols = c(r_se, alpha_se), 
#                values_to = "paramSE", 
#                names_to = "param", 
#                names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
#   select(-r_sim, -alpha_sim, -N0_sim, -r_est, -alpha_est)
# 
# ricDat_new_extinctAll_long <- left_join(paramEstLong_ext, paramSimLong_ext) %>% 
#   left_join(paramSELong_ext)
# 
# #calculate standardized parameter estimates
# ricDat_new_extinctAll_long <- ricDat_new_extinctAll_long %>%
#   mutate("paramDiff" = (paramEst - paramSim)/abs(paramSim))
# 
# # filter for low and high autocor
# ricDat_new_extinctAll_long[ricDat_new_extinctAll_long$actAutoCorr <=0.3 & !is.na(ricDat_new_extinctAll_long$actAutoCorr), "missingness"] <- "MAR: Low AC"
# ricDat_new_extinctAll_long[ricDat_new_extinctAll_long$actAutoCorr >0.3 & ricDat_new_extinctAll_long$actAutoCorr <0.6 & !is.na(ricDat_new_extinctAll_long$actAutoCorr), "missingness"] <- "MAR: Med. AC"
# ricDat_new_extinctAll_long[ricDat_new_extinctAll_long$actAutoCorr  >= 0.6 & !is.na(ricDat_new_extinctAll_long$actAutoCorr), "missingness"] <- "MAR: High AC"

## add the regular data to the extinct data data frame
ricDat_new_long_all <- ricDat_new_long #ricDat_new_extinctAll_long %>% 
#rename(SimNumber = simNumber, autoCorr = actAutoCorr, propMiss = actPropMiss) %>% 
#bind_rows(ricDat_new_long)

## save data for later (regular runs + runs for extinct data)
# # save data to file for use later...
#write_rds(ricDat_new_long_all, file = "./data/model_results/rickerRegAndExtinct_sim_ModelResultsLong.rds")

# Make Figures ------------------------------------------------------------
##compare the parameter estimates to the true values 
alpha_plot <- ggplot(data = ricDat_new_long_all[ricDat_new_long_all$param=="alpha",]) +
  geom_point(aes(x =  paramSim, y = paramEst)) +
  geom_abline(aes(slope = 1, intercept = 0), col = "grey") +
  facet_grid(~type, scales = "free") +
  ggtitle("alpha")
r_plot <- ggplot(data = ricDat_new_long_all[ricDat_new_long_all$param=="r",]) +
  geom_point(aes(x =  paramSim, y = paramEst)) +
  geom_abline(aes(slope = 1, intercept = 0), col = "grey") +
  facet_grid(~type, scales = "free")  +
  ggtitle("r")
ggpubr::ggarrange(alpha_plot, r_plot)

# generate summary statistics for each group of missingness and autocorrelation
ricDat_new_lines <- ricDat_new_long_all %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Med. AC", "MAR: Low AC")) %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff_abs, na.rm = TRUE),
            paramDiff_med = median(paramDiff_abs, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff_abs, na.rm = TRUE),
            n = length(paramDiff_abs)) %>% 
  filter(n  > 50)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.5)

# Figure of parameter recovery (mean and sd in separate panels) 
# figure of means for each model type and level of missingness (with shortened x-axis)
(pois_sim_MedsFig_trimmed <- ggplot(data = ricDat_new_lines, aes(x = amtMiss, y = paramDiff_med)) +
   facet_grid(~factor(param, levels = c( "alpha", "r")) 
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")),
              scales = "free_y") + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
   #geom_line(aes(x = amtMiss, y = paramDiff_med, color = as.factor(type)), position = position_dodge(width=0.03), lty = 2) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top", legend.title=element_blank())+
   ylab("Median of parameter bias across sims.")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
   scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))#+
 #ylim(c(-0.85,0.55))
)
# figure of SD for each model type and level of missingness
(pois_sim_SDFig_trimmed <- ggplot(data = ricDat_new_lines, aes(x = amtMiss, y = paramDiff_SD)) +
    facet_grid(~factor(param, levels = c( "alpha", "r")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")))+# ,
    #scales = "free_y") + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    geom_hline(aes(yintercept = 0), colour = "grey") + theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))    
) 

# put into one figure
pois_paramRecov_trimmed <- ggarrange(pois_sim_MedsFig_trimmed, pois_sim_SDFig_trimmed, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Poisson_medsSD_trimmed.png", width = 9, height = 4, units = "in", res = 700)
pois_paramRecov_trimmed
dev.off()



#
ricDat_new_verylong <- ricDat_new_long_all %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Med. AC", "MAR: Low AC")) %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 100)  

# make a figure like the one above, but without a trimmed x axis
# figure of means for each model type and level of missingness
(poiss_sim_MedsFig_reg<- ggplot(data = ricDat_new_verylong, aes(x = amtMiss, y = paramDiff_med)) +
    facet_grid(~factor(param, levels = c("alpha", "r")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC"))) + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
    #size=0.3, width=0, position = position_dodge(width=0.03))+
    #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Median of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))
  
)
# figure of SD for each model type and level of missingness
(poiss_sim_SDFig_reg<- ggplot(data = ricDat_new_verylong, aes(x = amtMiss, y = paramDiff_SD)) +
    facet_grid(~factor(param, levels = c("alpha", "r")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC","MAR: High AC"))) + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    #geom_ribbon(aes(ymin = 0, ymax = paramDiff_SD, fill = as.factor(type), color = as.factor(type)), alpha = .2, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))
) 
# put into one figure
poiss_paramRecov <- ggarrange(poiss_sim_MedsFig_reg, poiss_sim_SDFig_reg, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Poisson_medsSD.png", width = 9, height = 4, units = "in", res = 700)
poiss_paramRecov
dev.off()


ricDat_new_violin <- ricDat_new_long_all %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Med. AC", "MAR: Low AC")) %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)) %>% 
  filter(amtMiss <=0.5)


(pois_sim_violin <- ggplot(data = ricDat_new_violin) +
    facet_grid(~factor(param, levels = c( "alpha", "r")) 
               #~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")),
               ~ factor(type),
               scales = "free_y") + 
    geom_point(aes(x = as.factor(amtMiss), y = paramDiff, color = type), alpha = .8, position = position_dodge(width=0.07)) +
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    geom_violin(aes(x = as.factor(amtMiss), y = paramDiff, color = type), draw_quantiles = c(.50)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Mean standardized parameter estimate")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    ylim(c(-8,8))
)

