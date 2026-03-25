#/////////////////
# This script makes figures for Population simulated data 
# 4 December 2023
#/////////////////

# Load packages
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# read in MCAR data and prepare for figures------------------------------------------
# MCAR
ricker_MCAR <- readRDS("./data/model_results/pois_real_randMiss_drop_cc_MI_EM.rds")
ricker_MCAR$missingnessType <- "MCAR"
# MNAR
ricker_MNAR <- readRDS("./data/model_results/pois_real_minMaxMiss_drop_cc_MI_EM.rds")
ricker_MNAR$missingnessType <- "MNAR" 
ricker_MNAR$autocor = NA
ricker_MNAR$autocor_act = NA
ricker_MNAR <- ricker_MNAR %>% 
  select(names(ricker_MCAR))


allDat <- rbind(ricker_MCAR, ricker_MNAR) #allDat <- readRDS("./data/model_results/RickerForecast_resultTableAll.rds")

# extract parameter information from the list columns
ricDat <- rbind(
  #drop na fits
  cbind(allDat[,c("autocor", "index" , "prop_miss_act" , "autocor_act","LOO_rep", "missingnessType")], 
        map_df(allDat$drop_fits, function(x) {
          data.frame(
            "type" = "dropNA",
            "r_est" = x$estim["r"], 
            "alpha_est" = x$estim["alpha"],
            "r_se" = x$se["r"],
            "alpha_se" = x$se["alpha"],
            "status" = "good")
        }
    
        )#,
       # data.frame("listName" = names(allDat$drop_fits))
  ),
  #dropNA complete case fits
  cbind(allDat[,c("autocor", "index" , "prop_miss_act" , "autocor_act","LOO_rep", "missingnessType")], 
        map_df(allDat$cc_fits, function(x) {
          data.frame(
            "type" = "dropNA_cc",
            "r_est" = x$estim["r"], 
            "alpha_est" = x$estim["alpha"],
            "r_se" = x$se["r"],
            "alpha_se" = x$se["alpha"],
            "status" = "good")
        }
        
        )#,
        # data.frame("listName" = names(allDat$drop_fits))
  ),
  #EM fits
  cbind(allDat[,c("autocor", "index" , "prop_miss_act" , "autocor_act","LOO_rep", "missingnessType")], 
        map_df(allDat$EM_fits, function(x) {
          data.frame(
            "type" = "EM",
            "r_est" = x$estim["r"], 
            "alpha_est" = x$estim["alpha"],
            "r_se" = x$se["r"],
            "alpha_se" = x$se["alpha"],
            "status" = "good")
        }
        
        )#,
        # data.frame("listName" = names(allDat$drop_fits))
  ),
  # #DA fits
  # cbind(allDat[,c("autocor", "index" , "prop_miss_act" , "autocor_act","LOO_rep")], 
  #       map_df(allDat$, function(x) {
  #         data.frame(
  #           "type" = "EM",
  #           "r_est" = x$estim["r"], 
  #           "alpha_est" = x$estim["alpha"],
  #           "r_se" = x$se["r"],
  #           "alpha_se" = x$se["alpha"],
  #           "status" = "good")
  #       }
  #       
  #       )#,
  #       # data.frame("listName" = names(allDat$drop_fits))
  # )
  #,
  #MI fits
  cbind(allDat[,c("autocor", "index" , "prop_miss_act" , "autocor_act","LOO_rep", "missingnessType")], 
        map_df(allDat$MI_fits, function(x) {
        #   data.frame(
        #     "type" = "MI",
        #     "r_est" = x$estim["r"], 
        #     "alpha_est" = x$estim["alpha"],
        #     "r_se" = x$se["r"],
        #     "alpha_se" = x$se["alpha"],
        #     "status" = "good")
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
                      "type" = "MultipleImputations",
                      "r_est" = x$estim["r"],
                      "alpha_est" = x$estim["alpha"],
                      "r_se" = x$se["r"],
                      "alpha_se" = x$se["alpha"],
                      "status" = "good")
                  }
          }
          
        
        )#,
        # data.frame("listName" = names(allDat$drop_fits))
  )
)

## get data for models fit to non-missing data, and add as a reference to data frame generated above

#no missingness fits
noMissDat <- cbind(allDat[,c("autocor", "index" , "prop_miss_act" , "autocor_act","LOO_rep", "missingnessType")], 
      map_df(allDat$noMiss_fits, function(x) {
        data.frame(
          "type" = "noMiss",
          "r_est" = x$estim["r"], 
          "alpha_est" = x$estim["alpha"],
          "r_se" = x$se["r"],
          "alpha_se" = x$se["alpha"],
          "status" = "good")
      }
      
      )#,
      # data.frame("listName" = names(allDat$drop_fits))
) %>% 
  rename(r_est_noMiss = r_est,
         alpha_est_noMiss = alpha_est, 
         r_se_noMiss = r_se, 
         alpha_se_noMiss = alpha_se)

ricDat_new <- ricDat %>% 
  left_join(noMissDat %>% select(-type, -status), by = c("autocor", "index", "prop_miss_act", "autocor_act", "LOO_rep", "missingnessType"))

# # remove NAs (from missingness Limit Reached issue)
# ricDat_new <- ricDat[ricDat$status != "missingnessLimitReached",]
# ricDat_new <- unique(ricDat_new)

## note, "autoCorr" and "propMiss" in ricDat are actual values, not input values
# get autocorrelation amount 
#ricDat_new$actAutoCorr <- as.numeric(str_split(ricDat_new$listName, pattern = "_", simplify = TRUE)[,9])


# prepare for figures -----------------------------------------------------

#make into long data.frame 
paramEstLong <- ricDat_new %>% 
  #select(SimNumber, r_sim, actAutoCorr, actPropMiss, type, r_est, r_se, status, r_paramDiff) %>% 
  pivot_longer(cols = c(r_est, alpha_est), 
               values_to = "paramEst", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_est", "")) %>% 
  select(-r_est_noMiss, -alpha_est_noMiss, -r_se_noMiss, -alpha_se_noMiss)

paramSimLong <- ricDat_new %>% 
  pivot_longer(cols = c(r_est_noMiss, alpha_est_noMiss), 
               values_to = "paramSim", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_sim", "")) %>% 
  select(-r_est, -alpha_est, -r_se, -alpha_se, -r_se_noMiss, -alpha_se_noMiss) %>% 
  unique()

paramSELong <- ricDat_new %>% 
  pivot_longer(cols = c(r_se, alpha_se), 
               values_to = "paramSE", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_se", "")) %>% 
  select(-r_est, -alpha_est,-r_est, -r_est_noMiss, -alpha_est_noMiss, -r_se_noMiss , -alpha_se_noMiss) %>% 
  # remove the values for Expectation Maximization, since we don't have SE for that method
  filter(type != "ExpectationMaximization")


ricDat_new_long <- left_join(paramEstLong, paramSimLong) %>% 
  left_join(paramSELong)

#calculate standardized parameter estimates
ricDat_new_long <- ricDat_new_long %>%
  mutate("paramDiff" = (paramEst - paramSim)/abs(paramSim),
         "paramDiff_abs" = abs(paramEst - paramSim)/abs(paramSim))

# filter for low and high autocor
ricDat_new_long[ricDat_new$missingnessType == "MCAR" & ricDat_new_long$autocor_act <=0.3 & !is.na(ricDat_new_long$autocor_act), "missingnessType"] <- "MCAR: Low AC"
ricDat_new_long[ricDat_new$missingnessType == "MCAR" & ricDat_new_long$autocor_act >0.3 & ricDat_new_long$autocor_act <0.6 & !is.na(ricDat_new_long$autocor_act), "missingnessType"] <- "MCAR: Med. AC"
ricDat_new_long[ricDat_new$missingnessType == "MCAR" & ricDat_new_long$autocor_act  >= 0.6 & !is.na(ricDat_new_long$autocor_act), "missingnessType"] <- "MCAR: High AC"

# remove NAs (from missingness Limit Reached issue)
ricDat_new_long <- ricDat_new_long[ricDat_new_long$status != "missingnessLimitReached",]



 ricDat_new_long$type <- factor(ricDat_new_long$type, levels = c(#"DataAugmentation", 
   "dropNA_cc", 
                                                                 "dropNA", "EM", "MultipleImputations"))
 
# for data w/ no missingness, replace NA in propMiss with 0
# ricDat_new_long[ricDat_new_long$listName == "y_noMiss", "propMiss"] <- 0
# ricDat_new_long[ricDat_new_long$listName == "y_noMiss", "missingnessType"] <- "MCAR: Med. AC"
# ricDat_new_long_NOMISS <- ricDat_new_long[ricDat_new_long$listName == "y_noMiss", ]
# ricDat_new_long_NOMISS$missingnessType <- "MNAR"
# ricDat_new_long <- ricDat_new_long %>% 
#   rbind(ricDat_new_long_NOMISS)
# ricDat_new_long_NOMISS_ForMI <-  ricDat_new_long[ricDat_new_long$listName == "y_noMiss", ]
# ricDat_new_long_NOMISS_ForMI <- ricDat_new_long_NOMISS_ForMI %>% 
#   filter(type == "dropNA") %>% 
#   mutate(type = "MultipleImputations")
# ricDat_new_long <- ricDat_new_long %>% 
#   rbind(ricDat_new_long_NOMISS_ForMI)
# reformat data 
figDat_lines <- ricDat_new_long %>% 
  mutate(autoCor = round(as.numeric(autocor_act), 1), 
         amtMiss = prop_miss_act,
         amtMiss = replace(amtMiss, prop_miss_act <=0.3 & prop_miss_act > 0, 0.2),
         amtMiss = replace(amtMiss, prop_miss_act > 0.3 & prop_miss_act <=0.5, 0.4),
         amtMiss = replace(amtMiss, prop_miss_act > 0.5, 0.6)#,
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
  mutate(type=fct_relevel(type,c("dropNA", "dropNA_cc" ,"MultipleImputations","EM"#,"DataAugmentation"
                                 )))


figDat_lines2<-figDat_lines%>% filter(missingnessType %in% c("MCAR: Med. AC", "MNAR"))


# parameter recovery bias -------------------------------------------------
param.labs <- c("r", "alpha")
names(param.labs) <- c("r", expression(alpha))
figDat_lines2 <- figDat_lines2 %>% 
  mutate(figDat_lines2 = replace(param, param == "r", "'r'")) %>% 
  mutate(figDat_lines2 = replace(param, param == "alpha", 'alpha')) %>% 
  mutate(amtMiss = as.numeric(amtMiss))

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
  # xlim(c(-0.03,0.65)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73","#8c8c8c"#, "#CC79A7"
                                 ),
                        labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputation","Expectation Maximization"#,# "Data Augmentation"
                                   ))
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
   #xlim(c(-0.03,0.65)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73","#8c8c8c"#, "#CC79A7"
   ),
   labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputation","Expectation Maximization"#,# "Data Augmentation"
   ))
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
  mutate(autocor_act = as.numeric(autocor_act),
         prop_miss_act = as.numeric(prop_miss_act)) %>% 
  filter(param != "sigma",
         param != "intercept",
         prop_miss_act <=.65) %>% 
  mutate(autocor_act = round(as.numeric(autocor_act), 1), 
         amtMiss = prop_miss_act,
         amtMiss = replace(amtMiss, prop_miss_act <=0.3 & prop_miss_act > 0, 0.2),
         amtMiss = replace(amtMiss, prop_miss_act > 0.3 & prop_miss_act <=0.5, 0.4),
         amtMiss = replace(amtMiss, prop_miss_act > 0.5, 0.6)
  ) %>% 
  group_by(missingnessType, type, param, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN) 



# figDat_cov <- figDat_cov %>% 
#   mutate(type=fct_relevel(type,c("dropNA", "CompleteCaseDropNA" ,"MultipleImputations","DataAugmentation"))) ## No expectation maximization for coverage###


figDat_cov2<-figDat_cov%>% filter(missingnessType %in% c("MCAR: Med. AC", "MNAR"))



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
   # xlim(c(-0.03,0.65)) + 
    #scale_colour_brewer(palette = "Dark2",labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputations","Data Augmentation" ))
    
    #ylim(c(0,.3)) +
    scale_color_discrete(type = c("#D55E00","#E69F00", "#009E73"#, "#CC79A7"
                                  ),
                         labels = c( "Data Deletion-Complete","Data Deletion-Simple","Multiple Imputation"#, "Data Augmentation"
                                     )
                                     )
  #                        
)

# put figures together long version----------------------------------------------------

poiss_paramRecovery_bias_MAR2<-poiss_paramRecovery_bias_MAR+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.text = element_text(size=7))+guides(color = guide_legend(nrow = 2))

poiss_paramRecovery_SE_MAR2<-poiss_paramRecovery_SE_MAR+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.text = element_text(size=7))+guides(color = guide_legend(nrow = 2))

poiss_paramRecovery_coverage_MAR2<-poiss_paramRecovery_coverage_MAR+theme(legend.text = element_text(size=7))+guides(color = guide_legend(nrow = 2))

(poiss_paramRecovMAR <- ggarrange(poiss_paramRecovery_bias_MAR2, poiss_paramRecovery_SE_MAR2, 
                                  poiss_paramRecovery_coverage_MAR2, common.legend = TRUE, ncol = 1))

## save results
png(file = "./figures/parameterRecoveryPoisson_Real_long.png", width = 5, height = 8, units = "in", res = 700)
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

