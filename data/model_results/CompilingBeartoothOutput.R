#/////////////////////
# Compile Beartooth outputs into one csv files
# (with real and simulated datasets)
# 22 August 2023
#/////////////////////
library(tidyverse)

# gauss_sim_MAR_arima models ----------------------------------------------
# read in first group of output files (stored outside of the Git)
fileNames_A <- list.files("../BeartoothOutputs/gauss_sim_MAR_arima_modResults_A/")

for (i in 1:length(fileNames_A)) {
  assign(x = "temp", 
         value = read.csv(paste0("../BeartoothOutputs/gauss_sim_MAR_arima_modResults_A/",fileNames_A[i])))
  if (i == 1){
    outData_A <- temp
  } else {
    outData_A <- rbind(outData_A, temp)
  }
}

## add back in parameter info
params <- readRDS("./data/missingDatasets/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_A_final <- left_join(outData_A, params, by = c("simName" = "SimNumber"))

## read in the group of output files (stored outside of the Git)
fileNames_B <- list.files("../BeartoothOutputs/gauss_sim_MAR_arima_modResults_B/")

for (i in 1:length(fileNames_B)) {
  assign(x = "temp", 
         value = read.csv(paste0("../BeartoothOutputs/gauss_sim_MAR_arima_modResults_B/",fileNames_B[i])))
  
  if (i == 1){
    outData_B <- temp
  } else {
    outData_B <- rbind(outData_B, temp)
  }
}

## add back in parameter info
params <- readRDS("./data/missingDatasets/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_B_final <- left_join(outData_B, params, by = c("simName" = "SimNumber"))

## combine A and B into one d.f and remove unnecessary columns
outData_MAR_arima <- rbind(outData_A_final, outData_B_final) %>% 
  select(-CurSim, -missingnessVersion) %>% 
  mutate("2.5%" = NA, "50%" = NA, "97.5%" = NA) %>% 
  select("simName", "missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim") %>% 
  mutate(param = replace(param, param == "ar1", "phi")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim")


# gauss_sim_MNAR_arima models ---------------------------------------------------

# read in the group of output files (stored outside of the Git)
fileNames_MNAR <- list.files("../BeartoothOutputs/gauss_sim_MNAR_arima_modResults/")

for (i in 1:length(fileNames_MNAR)) {
  assign(x = "temp", 
         value = read.csv(paste0("../BeartoothOutputs/gauss_sim_MNAR_arima_modResults/",fileNames_MNAR[i])))
  
  if (i == 1){
    outData_MNAR <- temp
  } else {
    outData_MNAR <- rbind(outData_MNAR, temp)
  }
}


## add back in parameter info
params <- readRDS("./data/missingDatasets/forBeartooth/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_MNAR_final <- left_join(outData_MNAR, params, by = c("CurSim" = "SimNumber"))

# change curSim to "simNumber" (the same thing)
outData_MNAR_arima <- outData_MNAR_final %>% 
  rename("simName" = "CurSim") %>%
  mutate("2.5%" = NA, "50%" = NA, "97.5%" = NA) %>% 
  select("simName", "missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim") %>%  
  mutate(param = replace(param, param == "ar1", "phi")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim")


#  gauss_sim_MAR_brms models ----------------------------------------------

brms_MAR_A<- read_csv("data/model_results/00_combined_gauss_sim_randMiss_A.csv", 
                      show_col_types = FALSE) #379999 X  10#
brms_MAR_B<- read_csv("data/model_results/00_combined_gauss_sim_randMiss_B.csv", 
                      show_col_types = FALSE) # 379999 X 10#
# combine together
outData_MAR_brms <- rbind(brms_MAR_A, brms_MAR_B) %>% 
  rename("param" = "parameter", "value" = "mean", "SE" = "sd") %>% 
  select("missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%", "run_no") %>% 
  filter( missingness != "missingness") %>%  # remove rows that have column names (??)
  mutate(run_no = as.numeric(run_no)) 

# change run number to simulation numbers
simDF <- data.frame("run_no" = 1:5000, 
                    "simName" = rep.int(1:1000, times = 5))
outData_MAR_brms <- outData_MAR_brms %>% 
  left_join(simDF) %>% 
  select("simName", "missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>% 
  left_join(params, by = c("simName" = "SimNumber")) %>% 
  mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
  mutate(param = replace(param, param == "b_light", "light")) %>% 
  mutate(param = replace(param, param == "b_discharge", "discharge")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim")

# gauss_sim_MNAR_brms models ----------------------------------------------
brms_MNAR <- read_csv("data/model_results/00_combined_gauss_sim_minMaxMiss.csv", show_col_types = FALSE) # 80999 X 10 # 

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
  select("simName", "missingprop_autocor", "missingness", "type", "param", "value", "SE", "2.5%", "50%", "97.5%") %>% 
  left_join(params, by = c("simName" = "SimNumber")) %>% 
  mutate(missingness = "MNAR") %>% # change "MAR" to "MNAR" (was a mistake)
  mutate(param = replace(param, param == "b_Intercept", "intercept")) %>% 
  mutate(param = replace(param, param == "b_light", "light")) %>% 
  mutate(param = replace(param, param == "b_discharge", "discharge")) %>% 
  rename("intercept_sim" = "beta1_sim", 
         "light_sim" = "beta2_sim", 
         "discharge_sim" = "beta3_sim")


# save model outputs ------------------------------------------------------

## combine all of the model results for gaussian simulated data
outData_gauss_sim <- rbind(outData_MAR_arima, outData_MNAR_arima, outData_MAR_brms, outData_MNAR_brms)

saveRDS(outData_gauss_sim, file = "./data/model_results/gauss_sim_ModelResults.rds")


