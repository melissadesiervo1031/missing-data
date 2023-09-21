#/////////////////////
# Compile Beartooth outputs into one csv files
# (with real and simulated datasets)
# 22 August 2023
#/////////////////////

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
params <- readRDS("./data/missingDatasets/forBeartooth/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_A_final <- left_join(outData_A, params, by = c("simName" = "SimNumber"))
saveRDS(outData_A_final, file = "./data/BeartoothOutputData/gaussSim_MAR_A_arimaOut.rds")

# read in the group of output files (stored outside of the Git)
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
params <- readRDS("./data/missingDatasets/forBeartooth/gauss_sim_params.rds")
names(params) <- c("SimNumber", "phi_sim", "beta1_sim", "beta2_sim", "beta3_sim")

outData_B_final <- left_join(outData_B, params, by = c("simName" = "SimNumber"))
saveRDS(outData_B_final, file = "./data/BeartoothOutputData/gaussSim_MAR_B_arimaOut.rds")



# for MNAR arima output ---------------------------------------------------


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
saveRDS(outData_MNAR_final, file = "./data/BeartoothOutputData/gaussSim_MNAR_arimaOut.rds")

