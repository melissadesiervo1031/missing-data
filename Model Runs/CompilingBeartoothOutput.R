#/////////////////////
# Compile Beartooth outputs into one csv files
# (with real and simulated datasets)
# 22 August 2023
#/////////////////////

# read in first group of output files (stored outside of the Git)
fileNames_A <- list.files("../BeartoothOutputs/gauss_sim_randMiss_modResults_A_NEWOLD//")

for (i in 1:length(fileNames_A)) {
  assign(x = "temp", 
         value = read.csv(paste0("../BeartoothOutputs/gauss_sim_randMiss_modResults_A_NEWOLD/",fileNames_A[i])))
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

# # iterations that ran
# fileNumbers <- data.frame("iteration" = sort(as.numeric(str_extract(string = fileNames_A, pattern = "[0-9]*"))),
#                           "source" = "BeartoothIndex")
# 
# allNumbers <- data.frame("iteration" = as.numeric(1:5000),
#                          "source" = "actualIndex")
# 
# test <- full_join(fileNumbers, allNumbers, by = "iteration")
# # get missing iterations
# missingIterations <- test[is.na(test$source.x),"iteration"]

# read in second group of output files (stored outside of the Git)
fileNames_AA <- list.files("../BeartoothOutputs/gauss_sim_randMiss_modResults_AA/")

for (i in 1:length(fileNames_AA)) {
  assign(x = "temp", 
         value = read.csv(paste0("../BeartoothOutputs/gauss_sim_randMiss_modResults_AA/",fileNames_AA[i])))
  
  if (i == 1){
    outData_AA <- temp
  } else {
    outData_AA <- rbind(outData_AA, temp)
  }
}


# read in the group of output files (stored outside of the Git)
fileNames_B <- list.files("../BeartoothOutputs/gauss_sim_randMiss_modResults_B/")

for (i in 1:length(fileNames_B)) {
  assign(x = "temp", 
         value = read.csv(paste0("../BeartoothOutputs/gauss_sim_randMiss_modResults_B/",fileNames_B[i])))
  
  if (i == 1){
    outData_B <- temp
  } else {
    outData_B <- rbind(outData_B, temp)
  }
}

