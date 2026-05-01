#/////////////////
# This script makes figures for Population simulated data 
# 4 December 2023
#/////////////////

# Load packages
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# # read in MCAR data and prepare for figures------------------------------------------
# ricDat_tempA <- readRDS("./data/model_results/RickerA_resultTableRev1.rds")
# ricDat_tempB <- readRDS("./data/model_results/RickerB_resultTableRev1.rds")
# ricDat_temp <- rbind(ricDat_tempA, ricDat_tempB)
# 
# ## for now, remove rows where simulated population went extinct (was 13% of data!)
# # badRows <- as.vector(sapply(ricDat_tempNew$drop_fits, function(x)
# #   ifelse(sum(str_detect(names(x), pattern = "cause"))>=1, 
# #              yes = "bad", 
# #              no = "good")
# # ))
# # ricDat_tempNew <- ricDat_tempNew[rownames(ricDat_tempNew) %in% which(badRows == "good"),]
# # rename columns 
# ricDat_tempNew <- ricDat_temp %>% 
#   rename(r_sim = r, alpha_sim = alpha, N0_sim = N0)
# # put input autoCor and propMiss data into one column
# ricDat_tempNew <- ricDat_tempNew %>% 
#   mutate(input_args = paste0("a=", actAutoCorr, "_", "p=", actPropMiss)) %>% 
#   select(-autoCorr, -propMiss) %>% 
#   rename("autoCorr" = "actAutoCorr", "propMiss" = "actPropMiss")
# 
# # extract parameter information from the list columns
# ricDat <- rbind(
#   #drop na fits
#   cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_tempNew$drop_fits, function(x) 
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
#         data.frame("listName" = names(ricDat_tempNew$drop_fits))
#   ),
#   #dropNA complete case fits
#   cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_tempNew$cc_fits, function(x) 
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
#         data.frame("listName" = names(ricDat_tempNew$drop_fits))
#   ),
#   #EM fits
#   cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_tempNew$EM_fits, function(x) 
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
#         data.frame("listName" = names(ricDat_tempNew$drop_fits))
#         
#   ),
#   #DA fits
#   cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_tempNew$DA_fits, function(x) 
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
#               "type" = "DataAugmentation",
#               "r_est" = x$estim["r"], 
#               "alpha_est" = x$estim["alpha"],
#               "r_se" = x$se["r"],
#               "alpha_se" = x$se["alpha"],
#               "status" = "good")
#           }
#         ),
#         data.frame("listName" = names(ricDat_tempNew$drop_fits))
#         
#   )
#   #,
#   # #MI fits
#   # cbind(ricDat_tempNew[,c("SimNumber", "r_sim", "alpha_sim", 
#   #                      "N0_sim","input_args", "autoCorr", "propMiss")], 
#   #       map_df(ricDat_tempNew$MI_fits, function(x) 
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
#   #             "type" = "MultipleImputations",
#   #             "r_est" = x$estim["r"], 
#   #             "alpha_est" = x$estim["alpha"],
#   #             "r_se" = x$se["r"],
#   #             "alpha_se" = x$se["alpha"],
#   #             "status" = "good")
#   #         }
#   #       ),
#   #       data.frame("listName" = names(ricDat_tempNew$drop_fits))
#   #       
#   # )
# )
# 
# # remove NAs (from missingness Limit Reached issue)
# ricDat_new <- ricDat[ricDat$status != "missingnessLimitReached",]
# ricDat_new <- unique(ricDat_new)
# ## note, "autoCorr" and "propMiss" in ricDat are actual values, not input values
# # get autocorrelation amount 
# #ricDat_new$actAutoCorr <- as.numeric(str_split(ricDat_new$listName, pattern = "_", simplify = TRUE)[,9])
# 
# ## add MI data in (was run independently)
# # load MI model run data
# ricDat_tempB_MI <- readRDS("./data/model_results/RickerB_resultTableMIRev1.rds")
# ricDat_tempA_MI <- readRDS("./data/model_results/RickerA_resultTableMIRev1.rds")
# names(ricDat_tempA_MI) <- names(ricDat_tempB_MI)
# ricDat_temp_MI <- rbind(ricDat_tempA_MI, ricDat_tempB_MI)
# MI_temp <- ricDat_temp_MI
# # ## load MI extinction data (remove these values from the model results)
# # extinctMI <- read.csv("./data/model_results/ricker_Sim_reruns/MI_reruns/extinctMI_ALL.csv") %>% 
# #   select(-X.1, -X) %>% 
# #   filter(!is.na(simNumber)) %>% 
# #   rename(SimNumber = simNumber,
# #          actAutoCorr = autocorr_act, 
# #          autoCorr = autocorrs, 
# #          actPropMiss = prop_miss) %>% 
# #   mutate(extinct = TRUE)
# # 
# # MI_tempTemp <- MI_temp %>% 
# #   left_join(extinctMI %>% 
# #               select(SimNumber, actPropMiss, actAutoCorr, r, alpha, N0, extinct))
# 
# 
# # remove runs where population went extinct 
# good <- !str_locate(MI_temp$estim_r, pattern = " ")[,1]
# good[is.na(good)] <- TRUE
# 
# MI_temp <- MI_temp[good,]
# 
# test <- MI_temp %>% 
#   mutate(input_args = paste0("a=",autoCorr,"_p=",propMiss),
#          type = "MultipleImputations",
#          status = NA, 
#          listName = NA,
#          estim_r = as.numeric(estim_r)
#   ) %>% 
#   select(-id, -autoCorr, -propMiss) %>% 
#   rename(r_sim = r, 
#          alpha_sim = alpha,
#          N0_sim = N0, 
#          autoCorr = actAutoCorr,
#          propMiss = actPropMiss,
#          r_est = estim_r, 
#          alpha_est = estim_alpha,
#          r_se = se_r,
#          alpha_se = se_alpha
#   ) %>% 
#   select(SimNumber, r_sim, alpha_sim, N0_sim, input_args, autoCorr, 
#          propMiss, type, r_est, alpha_est, r_se, alpha_se, status, listName) %>% 
#   # remove runs that have an NA for parameter estimates
#   filter(!is.na(SimNumber)) %>% 
#   unique()
# 
# 
# ricDat_new <- ricDat_new %>% 
#   mutate(r_est = as.numeric(r_est)) %>% 
#   rbind(test)
# 
# ricDat_new <- ricDat_new %>% 
#   mutate(missingnessType = "MCAR")
# # actual zeros in autocorrelation are listed as NAs--fix this 
# ricDat_new[is.na(ricDat_new$autoCorr),"autoCorr"] <- 0
# 
# 
# # read in MNAR results ----------------------------------------------------
# ricDat_MNAR <- readRDS("./data/model_results/RickerMinMaxMissRev1.rds")
# 
# ricDat_MNAR <- ricDat_MNAR %>% 
#   rename(r_sim = r, alpha_sim = alpha, N0_sim = N0) %>% 
#   mutate(actAutoCorr = NA, autoCorr = NA)
# # put input autoCor and propMiss data into one column
# ricDat_MNARNew <- ricDat_MNAR %>% 
#   mutate(input_args = paste0("a=", actAutoCorr, "_", "p=", actPropMiss)) %>% 
#   select(-autoCorr, -propMiss) %>% 
#   rename("autoCorr" = "actAutoCorr", "propMiss" = "actPropMiss")
# 
# # extract parameter information from the list columns
# ricDatMNAR <- rbind(
#   #drop na fits
#   cbind(ricDat_MNARNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_MNARNew$drop_fits, function(x) 
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
#         data.frame("listName" = names(ricDat_MNARNew$drop_fits))
#   ),
#   #dropNA complete case fits
#   cbind(ricDat_MNARNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_MNARNew$cc_fits, function(x) 
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
#         data.frame("listName" = names(ricDat_MNARNew$drop_fits))
#   ),
#   #EM fits
#   cbind(ricDat_MNARNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_MNARNew$EM_fits, function(x) 
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
#         data.frame("listName" = names(ricDat_MNARNew$drop_fits))
#         
#   ),
#   #DA fits
#   cbind(ricDat_MNARNew[,c("SimNumber", "r_sim", "alpha_sim", 
#                           "N0_sim","input_args", "autoCorr", "propMiss")], 
#         map_df(ricDat_MNARNew$DA_fits, function(x) 
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
#               "type" = "DataAugmentation",
#               "r_est" = x$estim["r"], 
#               "alpha_est" = x$estim["alpha"],
#               "r_se" = x$se["r"],
#               "alpha_se" = x$se["alpha"],
#               "status" = "good")
#           }
#         ),
#         data.frame("listName" = names(ricDat_MNARNew$drop_fits))
#         
#   )
# )
# 
# # remove NAs (from missingness Limit Reached issue)
# ricDatMNAR_new <- ricDatMNAR[ricDatMNAR$status != "missingnessLimitReached",]
# ricDatMNAR_new <- unique(ricDatMNAR_new)
# ## note, "autoCorr" and "propMiss" in ricDat are actual values, not input values
# # get autocorrelation amount 
# #ricDatMNAR_new$actAutoCorr <- as.numeric(str_split(ricDatMNAR_new$listName, pattern = "_", simplify = TRUE)[,9])
# 
# ## add MI data in (was run independently)
# # load MI model run data
# ricDat_MNAR_MI <- readRDS("./data/model_results/RickerMinMaxMiss_MIRev1.rds")
# ricDat_MNAR_MI <- ricDat_MNAR_MI[,2:17] %>% 
#   mutate(actAutoCorr = NA, autoCorr = NA)
# # ## load MI extinction data (remove these values from the model results)
# # extinctMI <- read.csv("./data/model_results/ricker_Sim_reruns/MI_reruns/extinctMI_ALL.csv") %>% 
# #   select(-X.1, -X) %>% 
# #   filter(!is.na(simNumber)) %>% 
# #   rename(SimNumber = simNumber,
# #          actAutoCorr = autocorr_act, 
# #          autoCorr = autocorrs, 
# #          actPropMiss = prop_miss) %>% 
# #   mutate(extinct = TRUE)
# # 
# # MI_tempTemp <- MI_temp %>% 
# #   left_join(extinctMI %>% 
# #               select(SimNumber, actPropMiss, actAutoCorr, r, alpha, N0, extinct))
# 
# 
# # remove runs where population went extinct 
# good <- !str_locate(ricDat_MNAR_MI$estim_r, pattern = " ")[,1]
# good[is.na(good)] <- TRUE
# 
# ricDat_MNAR_MI <- ricDat_MNAR_MI[good,]
# 
# test_MNAR <- ricDat_MNAR_MI %>% 
#   mutate(input_args = paste0("a=",autoCorr,"_p=",propMiss),
#          type = "MultipleImputations",
#          status = NA, 
#          listName = NA,
#          estim_r = as.numeric(estim_r)
#   ) %>% 
#   select(-id, -autoCorr, -propMiss) %>% 
#   rename(r_sim = r, 
#          alpha_sim = alpha,
#          N0_sim = N0, 
#          autoCorr = actAutoCorr,
#          propMiss = actPropMiss,
#          r_est = estim_r, 
#          alpha_est = estim_alpha,
#          r_se = se_r,
#          alpha_se = se_alpha
#   ) %>% 
#   select(SimNumber, r_sim, alpha_sim, N0_sim, input_args, autoCorr, 
#          propMiss, type, r_est, alpha_est, r_se, alpha_se, status, listName) %>% 
#   # remove runs that have an NA for parameter estimates
#   filter(!is.na(SimNumber)) %>% 
#   unique()
# 
# 
# ricDatMNAR_new <- ricDatMNAR_new %>% 
#   mutate(r_est = as.numeric(r_est)) %>% 
#   rbind(test_MNAR)
# 
# 
# ricDatMNAR_new <- ricDatMNAR_new %>% 
#   mutate(missingnessType = "MNAR")
# # actual zeros in autocorrelation are listed as NAs--fix this 
# # ricDatMNAR_new[is.na(ricDatMNAR_new$autoCorr),"autoCorr"] <- 0
# 
# 
# # prepare for figures -----------------------------------------------------
# # add MCAR and MNAR data together
# ricDat_new <- ricDat_new %>% 
#   rbind(ricDatMNAR_new)
# 
# ## find those runs that have the same amount of missingness and autocorr, but aren't differentiated by the "listName"
# ricDat_new <- ricDat_new %>% 
#   rename(listName_Old = listName) %>% 
#   mutate(listName = paste0("pois_sim", SimNumber, "_", missingnessType, "_autoCorr_", round(autoCorr,2), "_propMissAct_", round(propMiss,2)))
# ricDat_new[is.na(ricDat_new$propMiss), "listName"] <- "y_noMiss"
# # find those 'listName"s that have >5 occurrences
# tooManyNames <- table(ricDat_new[ricDat_new$listName != "y_noMiss", "listName"])
# tooManyNames <- data.frame(tooManyNames[which(tooManyNames>5)])
# RicDat_new_normalNames <- ricDat_new[!(ricDat_new$listName %in% tooManyNames$Var1),]
# RicDat_new_tooManyNames <- apply(X = tooManyNames, MARGIN = 1, FUN = function(x) {
#   repeatVals <-  rep(c(1:(as.numeric(x["Freq"])/5)), length.out =  x["Freq"], times = 4) 
#   tempOut <- ricDat_new[ricDat_new$listName == x["Var1"],] 
#   tempOut$rep <- repeatVals
#   return(tempOut)
# }) %>% 
#   purrr::list_rbind()
# 
# RicDat_new_tooManyNames_2 <- RicDat_new_tooManyNames %>% 
#   mutate(listName = paste0(listName,"_", rep))
# 
# ricDatNew_new <- RicDat_new_normalNames %>% 
#   rbind(RicDat_new_tooManyNames_2 %>% select(-rep))
# ricDat_new <- ricDatNew_new
# ## save this output 
# saveRDS(ricDat_new, "./data/model_results/ricker_sim_ParamsForFigures.rds")

ricDat_new <- readRDS("./data/model_results/ricker_sim_ParamsForFigures.rds")
#make into long data.frame 
paramEstLong <- ricDat_new %>% 
  #select(SimNumber, r_sim, actAutoCorr, actPropMiss, type, r_est, r_se, status, r_paramDiff) %>% 
  pivot_longer(cols = c(r_est, alpha_est), 
               values_to = "paramEst", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_est", "")) %>% 
  select(-r_sim, -alpha_sim, -N0_sim, -r_se, -alpha_se)

paramSimLong <- ricDat_new %>% 
  pivot_longer(cols = c(r_sim, alpha_sim), 
               values_to = "paramSim", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_sim", "")) %>% 
  select(-r_est, -alpha_est, -N0_sim, -r_se, -alpha_se) %>% 
  unique()

paramSELong <- ricDat_new %>% 
  pivot_longer(cols = c(r_se, alpha_se), 
               values_to = "paramSE", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_se", "")) %>% 
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
ricDat_new_long[ricDat_new_long$autoCorr <=0.25 & !is.na(ricDat_new_long$autoCorr), "missingnessType"] <- "MCAR: Low AC"
ricDat_new_long[ricDat_new_long$autoCorr >0.25 & ricDat_new_long$autoCorr <0.65 & !is.na(ricDat_new_long$autoCorr), "missingnessType"] <- "MCAR: Med. AC"
ricDat_new_long[ricDat_new_long$autoCorr  >= 0.65 & !is.na(ricDat_new_long$autoCorr), "missingnessType"] <- "MCAR: High AC"


ricDat_new_long$type <- factor(ricDat_new_long$type, levels = c("DataAugmentation", "CompleteCaseDropNA", 
                                                                "dropNA", "ExpectationMaximization", "MultipleImputations"))

# for data w/ no missingness, replace NA in propMiss with 0
ricDat_new_long[ricDat_new_long$listName == "y_noMiss", "propMiss"] <- 0
ricDat_new_long[ricDat_new_long$listName == "y_noMiss", "missingnessType"] <- "MCAR: Med. AC"
ricDat_new_long_NOMISS <- ricDat_new_long[ricDat_new_long$listName == "y_noMiss", ]
ricDat_new_long_NOMISS$missingnessType <- "MNAR"
ricDat_new_long <- ricDat_new_long %>% 
  rbind(ricDat_new_long_NOMISS)
ricDat_new_long_NOMISS_ForMI <-  ricDat_new_long[ricDat_new_long$listName == "y_noMiss", ]
ricDat_new_long_NOMISS_ForMI <- ricDat_new_long_NOMISS_ForMI %>% 
  filter(type == "dropNA") %>% 
  mutate(type = "MultipleImputations")
ricDat_new_long <- ricDat_new_long %>% 
  rbind(ricDat_new_long_NOMISS_ForMI)
# reformat data 
figDat_lines <- ricDat_new_long %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = propMiss,
         amtMiss = replace(amtMiss, propMiss <=0.3 & propMiss > 0, 0.2),
         amtMiss = replace(amtMiss, propMiss > 0.3 & propMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, propMiss > 0.5, 0.6)#,
         #amtMiss = round(propMiss, 1)
         ) %>% 
  group_by(missingnessType, type, param, amtMiss) %>% 
  
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
  filter(amtMiss <=.65)
# make types in the 'correct' order


# reorder factor levels for plotting ##
figDat_lines <- figDat_lines %>% 
  mutate(type=fct_relevel(type,c("dropNA", "CompleteCaseDropNA" ,"MultipleImputations","ExpectationMaximization","DataAugmentation")))


figDat_lines2<-figDat_lines%>% filter(missingnessType %in% c("MCAR: Med. AC", "MNAR"))


# parameter recovery bias -------------------------------------------------
param.labs <- c("r", "alpha")
names(param.labs) <- c("r", expression(alpha))
figDat_lines2 <- figDat_lines2 %>% 
  mutate(figDat_lines2 = replace(param, param == "r", "'r'")) %>% 
  mutate(figDat_lines2 = replace(param, param == "alpha", 'alpha'))

(poiss_paramRecovery_bias_MAR <- ggplot(data = figDat_lines2, aes(x = amtMiss, y = paramDiff_med)) +
   ggh4x::facet_grid2(factor(param, levels = c( "r", "alpha"), labels = c("r", 'alpha'))
                      ~ factor(missingnessType, levels = c("MCAR: Med. AC", "MNAR"), labels =c("'Missing Completely at Random'", "'Missing Not at Random'")),
                      labeller =  label_parsed,
                      scales = "free_y")+
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
   ylab("Median Error")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.65)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73","#8c8c8c", "#CC79A7"),
                        labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputation","Expectation Maximization", "Data Augmentation"))
)

# parameter recovery SE ---------------------------------------------------
(poiss_paramRecovery_SE_MAR <- ggplot(data = figDat_lines2, aes(x = amtMiss, y = paramDiffAbsDiff_med)) +
    geom_hline(aes(yintercept = 0), colour = "grey") +
   ggh4x::facet_grid2(factor(param, levels = c( "r", "alpha"), labels = c("r", 'alpha'))
                      ~ factor(missingnessType, levels = c("MCAR: Med. AC", "MNAR"), labels =c("'Missing Completely at Random'", "'Missing Not at Random'")),
                      labeller =  label_parsed,
                      scales = "free_y")+
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
   ylab("Median Absolute Error")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.65)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73","8c8c8c", "#CC79A7"),
                        labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputation","Expectation Maximization", "Data Augmentation"))
 
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
  filter(#param != "sigma",
         #param != "intercept",
         propMiss <=.65,) %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = propMiss,
         amtMiss = replace(amtMiss, propMiss <=0.3 & propMiss > 0, 0.2),
         amtMiss = replace(amtMiss, propMiss > 0.3 & propMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, propMiss > 0.5, 0.6)
  ) %>% 
  group_by(missingnessType, type, param, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage)),# the total number of models 
            coverageTest = mean(coverage)
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN) %>% 
  ungroup()



figDat_cov <- figDat_cov %>% 
  filter(type != "ExpectationMaximization") %>% 
  mutate(type=fct_relevel(type,c("dropNA", "CompleteCaseDropNA" ,"MultipleImputations","DataAugmentation"))) ## No expectation maximization for coverage###


figDat_cov2<-figDat_cov%>% filter(missingnessType %in% c("MCAR: Med. AC", "MNAR"))
# 
# figDat_cov_temp %>% 
#   filter(#param != "sigma",
#       #param != "intercept",
#       propMiss <=.65,) %>% 
#       mutate(autoCor = round(autoCorr, 1), 
#              amtMiss =   propMiss,
#              amtMiss = replace(amtMiss, propMiss <=0.3 & propMiss > 0, 0.2),
#              amtMiss = replace(amtMiss, propMiss > 0.3 & propMiss <=0.5, 0.4),
#              amtMiss = replace(amtMiss, propMiss > 0.5, 0.6)
#   ) %>% 
#   filter(missingnessType %in% c("MNAR", "MCAR: Med. AC")) %>% 
# ggplot() + 
#   facet_grid(rows = vars(param), cols = vars(missingnessType)) + 
#   geom_point(aes(x = amtMiss, y = as.numeric(coverage), col = type)) +
#   #geom_boxplot(aes(x = amtMiss, y = as.numeric(coverage), col = type))#+ 
#   geom_smooth(aes(y = as.numeric(coverage), x = amtMiss, col = type), method = "glm") + 
#   geom_point(data = figDat_cov2, aes(x = amtMiss, y = coveragePerc, col = type), pch = 2)

(poiss_paramRecovery_coverage_MAR <- ggplot(data = figDat_cov2, aes(x = amtMiss, y = coveragePerc)) +
     ggh4x::facet_grid2(factor(param, levels = c( "r", "alpha"), labels = c("r", 'alpha'))
                      ~ factor(missingnessType, levels = c("MCAR: Med. AC", "MNAR"), labels =c("'Missing Completely at Random'", "'Missing Not at Random'")),
                      labeller =  label_parsed,
                      scales = "free_y") + #geom_col(aes(x = amtMiss, y = coveragePerc, color = as.factor(type), fill = as.factor(type)), 
    # position = "dodge", alpha = .5) + 
    geom_hline(aes(yintercept = .95), colour = "grey") +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Coverage")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    xlim(c(-0.03,0.65)) + 
    #scale_colour_brewer(palette = "Dark2",labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputations","Data Augmentation" ))
    
    #ylim(c(0,.3)) +
    scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73", "#CC79A7"),
                         labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputation", "Data Augmentation"))
  #                        
)




# put figures together long version----------------------------------------------------


poiss_paramRecovery_bias_MAR2<-poiss_paramRecovery_bias_MAR+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.text = element_text(size=7))+guides(color = guide_legend(nrow = 2))

poiss_paramRecovery_SE_MAR2<-poiss_paramRecovery_SE_MAR+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.text = element_text(size=7))+guides(color = guide_legend(nrow = 2))

poiss_paramRecovery_coverage_MAR2<-poiss_paramRecovery_coverage_MAR+theme(legend.text = element_text(size=7))+guides(color = guide_legend(nrow = 2))

(poiss_paramRecovMAR <- ggarrange(poiss_paramRecovery_bias_MAR2, poiss_paramRecovery_SE_MAR2, 
                                  poiss_paramRecovery_coverage_MAR2,common.legend = TRUE, ncol = 1))

## save results
png(file = "./figures/parameterRecoveryPoisson_MCARlong.png", width = 5, height = 8, units = "in", res = 700)
poiss_paramRecovMAR
dev.off()


# (poiss_paramRecovMAR <- ggarrange(poiss_paramRecovery_bias_MAR2, poiss_paramRecovery_SE_MAR2, 
#                                   poiss_paramRecovery_coverage_MAR2, legend = FALSE, common.legend = TRUE, ncol = 1))
# # save the figure object itself for subsequent plotting
# saveRDS(poiss_paramRecovMAR, "./figures/parameterRecoveryPoiss_MCARlong_FIGUREOBJECT.rds")
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
#ricDat_new_long_all <- ricDat_new_long #ricDat_new_extinctAll_long %>% 
#rename(SimNumber = simNumber, autoCorr = actAutoCorr, propMiss = actPropMiss) %>% 
#bind_rows(ricDat_new_long)

## save data for later (regular runs + runs for extinct data)
# # save data to file for use later...
#write_rds(ricDat_new_long_all, file = "./data/model_results/rickerRegAndExtinct_sim_ModelResultsLong.rds")

