#/////////////////////
# Compile Beartooth outputs into one csv files
# (with real and simulated datasets)
# 22 August 2023
#/////////////////////
library(tidyverse)


# Make sure % missingness values are correct ----------------------------------------
# temp1 <- read_rds(file = "./data/missingDatasets/gauss_sim_randMiss_A.rds")
# temp2 <- read_rds(file = "./data/missingDatasets/gauss_sim_randMiss_B.rds")
# simDat_raw <- c(temp1, temp2)
# 
# # get true missingness and autocorrelation for each missing time series
# # reshape dataset list
# simDat_raw2 <- map(simDat_raw, function(x) {
#   temp <- data.frame("jdate" = paste0("day_",1:365), x$y$y_noMiss)
#   temp <- data.frame(temp[1:292,], x$y[2:16])
#   outDat <- list_rbind(apply(temp[2:16], MARGIN = 2, FUN = function(x) {
#     temp_2 <- cbind(temp[1],x) %>% 
#       select(jdate, 2) %>% 
#       pivot_wider( names_from = jdate, values_from = 2)
#   }
#   ))
#   rownames(outDat) <- names(temp)[2:16]
#   return(outDat)
# }
# )
# # assign the appropriate names 
# simDat_names <- map(simDat_raw, function(x) {
#   temp <- data.frame(as.vector(names(x$y))[2:16], "names" = 1:15) 
# }
# )  %>% 
#   list_rbind()
# 
# simDat_raw3 <- simDat_raw2 %>% 
#   list_rbind()
# # add back in names w/ amount missing and amount autocor
# simDat_raw3$sim_miss_Names <- simDat_names$as.vector.names.x.y..
# # add back in simulation name
# simDat_raw3$sim_Name_long <- unlist(apply(data.frame(names(simDat_raw)), 1, 
#                                           function(x) rep.int(x, times = 15) , simplify = FALSE) )
# # extract simulation name
# simDat_raw3$sim_Names <- str_sub(simDat_raw3$sim_Name_long, 
#                                  start = 7, 
#                                  end = sapply(
#                                    str_locate_all(simDat_raw3$sim_Name_long, pattern = "_"), 
#                                    function(x) x[2,2]
#                                  )-1
# )
# 
# # extract name autocorrelation
# simDat_raw3$names_autoCorr <-  str_split(simDat_raw3$sim_miss_Names, pattern =  "_", simplify = TRUE)[,4] %>% 
#   as.numeric()
# # extract name amount missingness
# simDat_raw3$names_amtMiss <- str_split(simDat_raw3$sim_miss_Names, pattern =  "_", simplify = TRUE)[,2] %>% 
#   as.numeric()
# # calculate true amount of missingness
# simDat_raw3$true_amtMiss <- simDat_raw3 %>% 
#   select(1:292) %>% 
#   apply(MARGIN = 1, function(x) {
#     round(sum(is.na(x))/292,2)
#   })
# 
# # calculate true amount of autocorrelation
# simDat_raw3$true_autoCorr <- simDat_raw3 %>% 
#   select(1:292) %>% 
#   apply(MARGIN = 1, function(x) {
#     # change ts to 1 (not missing) and 0 (missing)
#     x[which(!is.na(x))] <- 1
#     x[which(is.na(x))] <- 0
#     round(acf(x = unlist(x), plot = FALSE, lag.max = 5, na.action = na.pass)$acf[,,1][2],2)
#   })
# 
# ## values are correct! phew

# gauss_sim_MAR_arima models ----------------------------------------------
# read in first group of output files (stored outside of the Git)
# output from arima model runs
outData_A <- read.csv("./data/model_results/gauss_sim_randMiss_modelResults_A/AllParams_arima.csv")

# for (i in 1:1000) {
#   assign(x = "temp", 
#          value = read.csv(paste0("./data/model_results/gauss_sim_randMiss_modelResults_A/",fileNames_A[i])))
#   if (i == 1){
#     outData_A <- temp
#   } else {
#     outData_A <- rbind(outData_A, temp)
#   }
# }

## add back in parameter info
params <- readRDS("./data/missingDatasets/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_A_final <- left_join(outData_A, params, by = c("sim_no" = "SimNumber"))

## read in the group of output files (stored outside of the Git)
outData_B <- read.csv("./data/model_results/gauss_sim_randMiss_modelResults_B/AllParams_arima.csv")

# for (i in 1:1000) {
#   assign(x = "temp", 
#          value = read.csv(paste0("./data/model_results/gauss_sim_randMiss_modelResults_B/",fileNames_B[i])))
#   
#   if (i == 1){
#     outData_B <- temp
#   } else {
#     outData_B <- rbind(outData_B, temp)
#   }
# }

## add back in parameter info
params <- readRDS("./data/missingDatasets/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_B_final <- left_join(outData_B, params, by = c("sim_no" = "SimNumber"))

## combine A and B into one d.f and remove unnecessary columns
outData_MAR_arima <- rbind(outData_A_final, outData_B_final) %>% 
  #select(-curSim) %>% 
  mutate("2.5%" = NA, "50%" = NA, "97.5%" = NA) %>% 
  select("sim_no", "curSim", "missingprop_autocor", "missingness", "type", "parameters", "param_value", "param_se", 
         "2.5%", "50%", "97.5%", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim") %>% 
  mutate(parameters = replace(parameters, parameters == "ar1", "phi")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim", 
         "simName" = "sim_no")


# ggplot(testDat) + 
#   geom_boxplot(aes(x = parameters, y = param_value, col = type))

# gauss_sim_MNAR_arima models ---------------------------------------------------

# read in the group of output files (stored outside of the Git)
  outData_MNAR_arima <- read.csv("./data/model_results/gauss_sim_minMax_modelResults/AllParams_arima.csv")

# for (i in 1:length(fileNames_MNAR)) {
#   assign(x = "temp", 
#          value = read.csv(paste0("./data/model_results/gauss_sim_minMax_modelResults/",fileNames_MNAR[i])))
#   
#   if (i == 1){
#     outData_MNAR <- temp
#   } else {
#     outData_MNAR <- rbind(outData_MNAR, temp)
#   }
# }


## add back in parameter info
params <- readRDS("./data/missingDatasets/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_MNAR_final <- left_join(outData_MNAR_arima, params, by = c("sim_no" = "SimNumber"))

# change curSim to "simNumber" (the same thing, in this case)
outData_MNAR_arima <- outData_MNAR_final %>% 
  mutate(simName = curSim) %>%
  mutate("2.5%" = NA, "50%" = NA, "97.5%" = NA) %>% 
  select("simName", "curSim", "missingprop_autocor", "missingness", "type", "parameters", "param_value", "param_se", "2.5%", "50%", "97.5%", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim") %>%  
  mutate(parameters = replace(parameters, parameters == "ar1", "phi")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim")

#  gauss_sim_MAR_brms models ----------------------------------------------

brms_MAR_A<- read_csv("./data/model_results/gauss_sim_randMiss_modelResults_A/brmsResults/AllParams_brms.csv", 
                      show_col_types = FALSE) %>% 
  filter(parameter != "Intercept") #379999 X  10#
brms_MAR_B<- read_csv("./data/model_results/gauss_sim_randMiss_modelResults_B/brmsResults/AllParams_brms.csv", 
                      show_col_types = FALSE) %>% 
  filter(parameter != "Intercept")  # 379999 X 10#
# combine together
outData_MAR_brms <- rbind(brms_MAR_A, brms_MAR_B
                          ) %>% 
 # brms_MAR_A %>% 
  rename("param" = "parameter", "value" = "mean", "SE" = "sd") %>% 
  select("missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%", "run_no") %>% 
  filter( missingness != "missingness") %>%  # remove rows that have column names (??)
  mutate(run_no = as.numeric(run_no)) 

# change run number to simulation numbers
simDF <- data.frame("run_no" = 1:5000, 
                    "simName" = rep.int(1:1000, times = 5))
outData_MAR_brms <- outData_MAR_brms %>% 
  left_join(simDF) %>% 
  mutate(curSim = NA) %>% 
  select("simName", "curSim", "missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>% 
  left_join(params, by = c("simName" = "SimNumber")) %>% 
  mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
  mutate(param = replace(param, param == "b_light", "light")) %>% 
  mutate(param = replace(param, param == "b_discharge", "discharge")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim",
         "parameters" = "param", 
         "param_value" = "value", 
         "param_se" = "SE")

# # gauss_sim_MAR_brms models w/ Normal Prior on Phi  -----------------------
# brms_MAR_A_norm <- read.csv("./data/model_results/gauss_sim_MAR_A_brms_results_normPrior.csv")
# brms_MAR_B_norm <- read.csv('./data/model_results/gauss_sim_MAR_B_brms_results_normPrior.csv')
# # add in simulation data
# simDF <- data.frame("run_no" = 1:5000, 
#                     "simName" = rep.int(1:1000, times = 5))
# 
# # reformat data to be consistent w/ uniform prior on phi version
# outData_MAR_brms_norm <- brms_MAR_A_norm %>% 
#   rbind(brms_MAR_B_norm) %>% 
#   left_join(simDF) %>% 
#   select(-run_no) %>% 
#   rename("2.5%" = "X2.5.", "50%" = "X50.", "97.5%" = "X97.5.",  "value" = "mean", "param" = "parameter", "SE" = "sd") %>% 
#   select(simName, missingprop_autocor, missingness, type, param, value, SE, "2.5%", "50%", "97.5%") %>% 
#   left_join(params, by = c("simName" = "SimNumber")) %>% #add back in simulation parameters
#   mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
#   mutate(param = replace(param, param == "b_light", "light")) %>% 
#   mutate(param = replace(param, param == "b_discharge", "discharge")) %>% 
#   rename("intercept_sim" = "beta1_sim", 
#          "light_sim" = "beta2_sim", 
#          "discharge_sim" = "beta3_sim")
# 
# # gauss_sim_MAR_brms models w/ Normal Prior on Phi w/ no boundaries -----------------------
# brms_MAR_A_normNB <- read.csv("./data/model_results/gauss_sim_MAR_A_brms_results_normPriorNB.csv")
# brms_MAR_B_normNB <- read.csv('./data/model_results/gauss_sim_MAR_B_brms_results_normPriorNB.csv')
# # add in simulation data
# simDF <- data.frame("run_no" = 1:5000, 
#                     "simName" = rep.int(1:1000, times = 5))
# 
# # reformat data to be consistent w/ uniform prior on phi version
# outData_MAR_brms_normNB <- brms_MAR_A_normNB %>% 
#   rbind(brms_MAR_B_normNB) %>% 
#   left_join(simDF) %>% 
#   select(-run_no) %>% 
#   rename("2.5%" = "X2.5.", "50%" = "X50.", "97.5%" = "X97.5.",  "value" = "mean", "param" = "parameter", "SE" = "sd") %>% 
#   select(simName, missingprop_autocor, missingness, type, param, value, SE, "2.5%", "50%", "97.5%") %>% 
#   left_join(params, by = c("simName" = "SimNumber")) %>% #add back in simulation parameters
#   mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
#   mutate(param = replace(param, param == "b_light", "light")) %>% 
#   mutate(param = replace(param, param == "b_discharge", "discharge")) %>% 
#   rename("intercept_sim" = "beta1_sim", 
#          "light_sim" = "beta2_sim", 
#          "discharge_sim" = "beta3_sim")

# gauss_sim_MNAR_brms models ----------------------------------------------
brms_MNAR <- read_csv("./data/model_results/gauss_sim_minMax_modelResults/brmsResults/AllParams_brms.csv", show_col_types = FALSE) %>% 
  filter(parameter != "Intercept") # 80999 X 10 # 

# combine together
outData_MNAR_brms <- brms_MNAR %>% 
  rename("param" = "parameter", "value" = "mean", "SE" = "sd") %>% 
  select("missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%", "run_no") %>% 
  filter( missingness != "missingness") %>%  # remove rows that have column names (??)
  mutate(run_no = as.numeric(run_no))

# change run number to simulation numbers
simDF <- data.frame("run_no" = 1:5000, 
                    "simName" = rep.int(1:1000, times = 5))
outData_MNAR_brms <- outData_MNAR_brms %>% 
  left_join(simDF) %>% 
  mutate(curSim = NA) %>% 
  select("simName", "curSim", "missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>% 
  left_join(params, by = c("simName" = "SimNumber")) %>% 
  mutate(missingness = "MNAR") %>% # change "MAR" to "MNAR" (was a mistake)
  mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
  mutate(param = replace(param, param == "b_light", "light")) %>% 
  mutate(param = replace(param, param == "b_discharge", "discharge")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim",
         "parameters" = "param", 
         "param_value" = "value", 
         "param_se" = "SE")

# save model outputs for Gaussian simulation models ----------------------------

## combine all of the model results for gaussian simulated data
outData_gauss_sim <- rbind(outData_MAR_arima, outData_MNAR_arima, outData_MAR_brms, outData_MNAR_brms)

## clean up, and calculate simulation data
outData_gauss_sim <- unique(outData_gauss_sim)

# make columns for "autocor" and "missingness"
outData_gauss_sim$autoCor <- outData_gauss_sim$missingprop_autocor %>% 
  str_extract(pattern = "0.[0-9]+$") %>% 
  as.numeric()
outData_gauss_sim[outData_gauss_sim$missingness=="MNAR", "autoCor"] <- NA
outData_gauss_sim$amtMiss <- outData_gauss_sim$missingprop_autocor %>% 
  str_extract(pattern = "0.[0-9]+") %>% 
  as.numeric
outData_gauss_sim <- outData_gauss_sim %>% 
  mutate(param_value = as.numeric(param_value),
         param_se = as.numeric(param_se))
# 
outData_gauss_sim[outData_gauss_sim$missingness == "MAR" & 
                    is.na(outData_gauss_sim$autoCor), "autoCor"] <- 0

# remove values for models fitted to time series with no missingness (Doesn't work for all model approaches)
gauss_sim_figDat <- outData_gauss_sim#[outData_gauss_sim$missingprop_autocor != "y_noMiss",]

# make names of parameters consistent
gauss_sim_figDat <- gauss_sim_figDat %>% 
  mutate(parameters = replace(parameters, parameters  %in% c("Intercept", "intercept"), "intercept"), 
         parameters = replace(parameters, parameters  %in% c("xreg1"), "light"), 
         parameters = replace(parameters, parameters  %in% c("xreg2"), "discharge"))

simDat <- gauss_sim_figDat %>% 
  select(simName, phi_sim, intercept_sim, light_sim, discharge_sim) %>% 
  pivot_longer(cols = c(phi_sim, intercept_sim, light_sim, discharge_sim), 
               names_to = "param", 
               values_to = "param_simVal",
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  unique()


# # get sim parameters from the raw data to compare
# simDat_raw_params <- lapply(simDat_raw, FUN = function(x) {
#   return(data.frame("phi" = x[["sim_params"]]$phi,
#              "beta1" = x[["sim_params"]]$beta[1],
#              "beta2" = x[["sim_params"]]$beta[2],
#              "beta3" = x[["sim_params"]]$beta[3]
#   ))
# })
# simDat_raw_params <- simDat_raw_params %>% 
#   purrr::list_rbind() %>% 
#   unique()
# simDat_raw_params$simName <- c(1:nrow(simDat_raw_params))

# remove columns for simulation data
gauss_sim_figDat <- gauss_sim_figDat %>% 
  select(-phi_sim, -intercept_sim, -light_sim, -discharge_sim)
# add back in simulation parameter data
gauss_sim_figDat <- gauss_sim_figDat %>% 
  rename("param" = "parameters") %>% 
  left_join(simDat, by = c("simName", "param"))

# calculate the standardized difference between parameter estimates and simulated values
gauss_sim_figDat <- gauss_sim_figDat %>% 
  dplyr::mutate(paramDiff = ((param_value - param_simVal)/abs(param_simVal)),
                paramDiff_absDiff = abs((param_value - param_simVal)/abs(param_simVal))) 


saveRDS(gauss_sim_figDat, file = "./data/model_results/gauss_sim_ModelResults.rds")

# # save model outputs for Gaussian simulation models w/ brms that use Normal priors w/ no Bounds----------------------------

# tidy gaussian real data -------------------------------------------------

# read in data
# get MCAR data
gauss_real_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/au_sable/gauss_auSable_real_MAR_arima_FORECASTvals.csv") %>%
  select(missingprop_autocor:run_no)
gauss_real_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/auSable/gauss_auSable_real_MAR_brms_FORECASTvals.csv") %>%
  select(missingprop_autocor:sd, missingness:run_no) %>%
  rename(parameters = parameter,
         param_value = mean,
         param_se = sd) %>% 
  filter(parameters != "Intercept") %>% 
  unique()

gauss_real_figDat <- gauss_real_arima %>%
  rbind(gauss_real_brms) %>%
  mutate(amtMiss = as.numeric(str_extract(missingprop_autocor, pattern = "[0-9/.0-9]+")),
         autoCor = as.numeric(str_extract(missingprop_autocor, pattern = "[0-9/.0-9]+$")),
         parameters = replace(parameters, parameters == "xreg1", "light"),
         parameters = replace(parameters, parameters == "b_light","light"),
         parameters = replace(parameters, parameters == "b_Intercept","intercept"),
         parameters = replace(parameters, parameters == "b_Q","Q"),
         parameters = replace(parameters, parameters == "xreg2","Q"),
         parameters = replace(parameters, parameters == "ar1","phi")
  )

# get MNAR data
gauss_real_arima_MNAR <- read.csv("./data/model_results/gauss_real_MNAR_arima_modResults/au_sable/gauss_auSable_real_MNAR_arima_FORECASTvals.csv") %>%
  select(missingprop_autocor:run_no)
gauss_real_brms_MNAR <- read.csv("./data/model_results/gauss_real_MNAR_brms_modResults/auSable/brmsvals.csv") %>%
  select(missingprop_autocor:sd, missingness:run_no) %>%
  rename(parameters = parameter,
         param_value = mean,
         param_se = sd)%>% 
  filter(parameters != "Intercept") %>% 
  unique()

gauss_real_figDat_MNAR <- gauss_real_arima_MNAR %>%
  rbind(gauss_real_brms_MNAR) %>%
  mutate(autoCor = NA,
         amtMiss = as.numeric(str_extract(missingprop_autocor, pattern = "[0-9/.0-9]+$"))
  ) %>%
  mutate(amtMiss = replace(amtMiss, is.na(amtMiss), 0.44)) %>%
  mutate(
    parameters = replace(parameters, parameters == "xreg1", "light"),
    parameters = replace(parameters, parameters == "b_light","light"),
    parameters = replace(parameters, parameters == "b_Intercept","intercept"),
    parameters = replace(parameters, parameters == "Intercept","intercept"),
    parameters = replace(parameters, parameters == "b_discharge","Q"),
    parameters = replace(parameters, parameters == "xreg2","Q"),
    parameters = replace(parameters, parameters == "ar1","phi")
  )
## add the MNAR data to the MAR data
gauss_real_figDat <- gauss_real_figDat %>%
  rbind(gauss_real_figDat_MNAR)
## get the data from models fit to the 'complete' dataset (has a few NAs), which will be the parameters that those fit to the missing datasets will be compared to
gauss_real_RefParams_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/au_sable/CompleteDataset_arimavals.csv") %>%
  select(missingprop_autocor:run_no)
gauss_real_RefParams_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/auSable/CompleteDataset_brmsvals.csv") %>%
  select(missingprop_autocor:sd, missingness:run_no) %>%
  rename(parameters = parameter,
         param_value = mean,
         param_se = sd) %>%
  filter(parameters != "Intercept")

gauss_real_RefParams <- gauss_real_RefParams_arima %>%
  rbind(gauss_real_RefParams_brms) %>%
  mutate(
    parameters = replace(parameters, parameters == "xreg1", "light"),
    parameters = replace(parameters, parameters == "b_light","light"),
    parameters = replace(parameters, parameters == "b_Intercept","intercept"),
    parameters = replace(parameters, parameters == "b_Q","Q"),
    parameters = replace(parameters, parameters == "xreg2","Q"),
    parameters = replace(parameters, parameters == "ar1","phi")
  ) %>%
  rename(param_simVal = param_value)

## add the parameters from the models fit to the 'complete' dataset as reference
gauss_real_figDat <- gauss_real_figDat %>%
  left_join(gauss_real_RefParams %>% select(parameters, param_simVal, type)) %>%
  # calculate difference between model-derived parameters w/ missingness and no missingness
  dplyr::mutate(paramDiff = ((param_value - param_simVal)/abs(param_simVal)),
                paramDiff_absDiff = abs((param_value - param_simVal)/abs(param_simVal)))
# save for later
saveRDS(gauss_real_figDat, "./data/model_results/gauss_real_ModelResults.rds")


# 
# ## combine all of the model results for gaussian simulated data
# outData_gauss_sim_normPrior <- rbind(outData_MAR_arima, 
#                                      outData_MNAR_arima, 
#                                      outData_MAR_brms_normNB, outData_MNAR_brms)
# 
# ## clean up, and calculate simulation data
# outData_gauss_sim_normPrior <- unique(outData_gauss_sim_normPrior)
# 
# # make columns for "autocor" and "missingness"
# outData_gauss_sim_normPrior$autoCor <- outData_gauss_sim_normPrior$missingprop_autocor %>% 
#   str_extract(pattern = "0.[0-9]+$") %>% 
#   as.numeric()
# outData_gauss_sim_normPrior[outData_gauss_sim_normPrior$missingness=="MNAR", "autoCor"] <- NA
# outData_gauss_sim_normPrior$amtMiss <- outData_gauss_sim_normPrior$missingprop_autocor %>% 
#   str_extract(pattern = "0.[0-9]+") %>% 
#   as.numeric
# outData_gauss_sim_normPrior <- outData_gauss_sim_normPrior %>% 
#   mutate(value = as.numeric(value),
#          SE = as.numeric(SE))
# # 
# outData_gauss_sim_normPrior[outData_gauss_sim_normPrior$missingness == "MAR" & 
#                               is.na(outData_gauss_sim_normPrior$autoCor), "autoCor"] <- 0
# 
# # remove values for models fitted to time series with no missingness (Doesn't work for all model approaches)
# gauss_sim_figDat_normPrior <- outData_gauss_sim_normPrior[outData_gauss_sim_normPrior$missingprop_autocor != "y_noMiss",]
# 
# simDat <- gauss_sim_figDat_normPrior %>% 
#   select(simName, phi_sim, intercept_sim, light_sim, discharge_sim) %>% 
#   pivot_longer(cols = c(phi_sim, intercept_sim, light_sim, discharge_sim), 
#                names_to = "param", 
#                values_to = "param_simVal",
#                names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
#   unique()
# 
# # remove columns for simulation data
# gauss_sim_figDat_normPrior <- gauss_sim_figDat_normPrior %>% 
#   select(-phi_sim, -intercept_sim, -light_sim, -discharge_sim)
# 
# gauss_sim_figDat_normPrior <- gauss_sim_figDat_normPrior %>% 
#   left_join(simDat, by = c("simName", "param"))
# 
# # calculate the standardized difference between parameter estimates and simulated values
# gauss_sim_figDat_normPrior <- gauss_sim_figDat_normPrior %>% 
# 
# 
#   dplyr::mutate(paramDiff = ((value - param_simVal)/abs(param_simVal)),
#                 paramDiff_absDiff = abs((value - param_simVal)/abs(param_simVal))) 
# 
# 
# saveRDS(gauss_sim_figDat_normPrior, file = "./data/model_results/gauss_sim_ModelResults_normPrior.rds")

# # gauss_real_MAR_arima models ----------------------------------------------
# # read in output file
# gauss_real_MAR_arima <- read.csv("./data/model_results/gauss_real_MAR_arima.csv")
# 
# ## combine A and B into one d.f and remove unnecessary columns
# outData_real_MAR_arima <- gauss_real_MAR_arima %>% 
#   select(-CurSim, -missingnessVersion) %>% 
#   mutate("2.5%" = NA, "50%" = NA, "97.5%" = NA) %>% 
#   select("missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>% 
#   mutate(param = replace(param, param == "ar1", "phi"))
# 
# # gauss_real_MNAR_arima models ---------------------------------------------------
# gauss_real_MNAR_arima <- read.csv("./data/model_results/gauss_real_MNAR_arima.csv")
# 
# # change curSim to "simNumber" (the same thing)
# outData_real_MNAR_arima <- gauss_real_MNAR_arima %>% 
#   select( -missingnessVersion) %>% 
#   mutate("2.5%" = NA, "50%" = NA, "97.5%" = NA) %>% 
#   select("missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>%  
#   mutate(param = replace(param, param == "ar1", "phi")) 
# 
# #  gauss_sim_MAR_brms models ----------------------------------------------
# gauss_real_MAR_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_results.csv")
# 
# # combine together
# outData_real_MAR_brms <- gauss_real_MAR_brms %>% 
#   rename("param" = "parameter", "value" = "mean", "SE" = "sd", "2.5%" = "X2.5.", "50%" = "X50.", "97.5%" = "X97.5.") %>% 
#   select("missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>% 
#   filter( missingness != "missingness") %>%  # remove rows that have column names (??)
#   mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
#   mutate(param = replace(param, param == "b_light", "light")) %>% 
#   mutate(param = replace(param, param == "b_Q", "discharge")) 
# 
# # gauss_real_MNAR_brms models ----------------------------------------------
# gauss_real_MNAR_brms <- read.csv("./data/model_results/gauss_real_MNAR_brms_modresults.csv")
# 
# # combine together
# outData_real_MNAR_brms <- gauss_real_MNAR_brms %>% 
#   rename("param" = "parameter", "value" = "mean", "SE" = "sd", "2.5%" = "X2.5.", "50%" = "X50.", "97.5%" = "X97.5.") %>% 
#   select("missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>% 
#   filter( missingness != "missingness") %>%  # remove rows that have column names (??)
#   mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
#   mutate(param = replace(param, param == "b_light", "light")) %>% 
#   mutate(param = replace(param, param == "b_Q", "discharge")) 
# 
# # save model outputs ------------------------------------------------------
# 
# ## combine all of the model results for gaussian real data
# outData_gauss_real <- rbind(outData_real_MAR_arima, outData_real_MNAR_arima, outData_real_MAR_brms, outData_real_MNAR_brms)
# 
# saveRDS(outData_gauss_real, file = "./data/model_results/gauss_real_ModelResults.rds")

