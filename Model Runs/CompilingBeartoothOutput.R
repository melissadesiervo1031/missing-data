#/////////////////////
# Compile Beartooth outputs into one csv files
# (with real and simulated datasets)
# 22 August 2023
#/////////////////////

# read in first group of output files (stored outside of the Git)
fileNames_AAA <- list.files("../BeartoothOutputs/gauss_sim_randMiss_modResults_AAA/")

for (i in 1:length(fileNames_AAA)) {
  assign(x = "temp", 
         value = read.csv(paste0("../BeartoothOutputs/gauss_sim_randMiss_modResults_AAA/",fileNames_AAA[i])))
  
  if (i == 1){
    outData_AAA <- temp
  } else {
    outData_AAA <- rbind(outData_AAA, temp)
  }
}

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



grep(pattern = "^[:num:]*", x = fileNames_A)
