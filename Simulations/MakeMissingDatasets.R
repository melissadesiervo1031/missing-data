#/////////////////////
# Generates datasets with various levels of missingness
# (with real and simulated datasets)
# 24 April 2023
#/////////////////////


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(lubridate)

# Load Functions ----------------------------------------------------------
# import the "makeMissing" function
source("./Functions/missing_data_functions.R")

### Note: throughout the simulations of added missingness, only add missingness to response variable (i.e. GPP)

# Gaussian Simulated Data -------------------------------------------------
# read in data
gauss_sim <- readRDS("./data/gauss_ar1_0miss_datasets.rds")

## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness
  tempOutList <- lapply(X = gauss_sim, FUN = function(x) 
    list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", 
                           autoCorr = inputAutocor[i]),
         "sim_params" <- x$sim_params)
  )
  # name the elements of the list with the amount of missingness 
  for (j in 1:length(tempOutList)) {
    names(tempOutList[[j]]) <- c("y", "sim_params")
    tempOutList[[j]]$y <- c(list("y_noMiss" = gauss_sim[[j]]$y), tempOutList[[j]]$y)
  }
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("gauss_sim_randMiss_autoCorr_0") , 
           value = tempOutList)
  } else {
    assign(x = paste0("gauss_sim_randMiss_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutList)
  }
}

# missing in min and max of data
gauss_sim_minMaxMiss <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "minMax"), 

       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_minMaxMiss)) {
  names(gauss_sim_minMaxMiss[[i]]) <- c("y", "sim_params")
  gauss_sim_minMaxMiss[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_minMaxMiss[[i]]$y)

}

# Gaussian Real Data (Pine River GPP data) --------------------------------
# read in data
dat <- read.csv('./data/NWIS_MissingTS_subset.csv')
mdat <- read.csv('data/NWIS_MissingTSinfo_subset.csv')

id <- mdat$site_name[4]
pr <- dat %>% filter(site_name == id) %>% select(date, GPP, light, Q, GPP.upper,GPP.lower) %>% mutate(Jdate= yday(date), light.rel = light/max(light))
gauss_real <- as.data.frame(pr)

## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness
  tempOutDf <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, 
                                           typeMissing = "random", 
                                           autoCorr = inputAutocor[i]))
  
  # name the elements of the Df with the amount of missingness 
  names(tempOutDf) <- paste0("GPP_",names(tempOutDf))
  tempOutDf <- cbind(gauss_real, tempOutDf)
  
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("gauss_real_randMiss_autoCorr_0") , 
           value = tempOutDf)
  } else {
    assign(x = paste0("gauss_real_randMiss_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutDf)
  
  }
}

## missing in min and max of data
gauss_real_minMaxMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, 
                                                   typeMissing = "minMax"))


names(gauss_real_minMaxMiss_TEMP) <- paste0("GPP_",names(gauss_real_minMaxMiss_TEMP))

gauss_real_minMaxMiss <- cbind(gauss_real, gauss_real_minMaxMiss_TEMP)


# Poisson Simulated Data --------------------------------------------------
# read in data
pois_sim <- readRDS("./data/ricker_0miss_datasets.rds")

# make missing data types

## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness
  tempOutList <- lapply(X = pois_sim, FUN = function(x) 
    list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = inputAutocor[i]), 
         "sim_params" <- x$sim_params)
  )
  # name the elements of the list with the amount of missingness 
  for (j in 1:length(tempOutList)) {
    names(tempOutList[[j]]) <- c("y", "sim_params")
    tempOutList[[j]]$y <- c(list("y_noMiss" = pois_sim[[j]]$y), tempOutList[[j]]$y)
  }
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("pois_sim_randMiss_autoCorr_0") , 
           value = tempOutList)
  } else {
    assign(x = paste0("pois_sim_randMiss_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutList)
  }
}

# missing in min and max of data
pois_sim_minMaxMiss <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "minMax"), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_minMaxMiss)) {
  names(pois_sim_minMaxMiss[[i]]) <- c("y", "sim_params")
  pois_sim_minMaxMiss[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_minMaxMiss[[i]]$y)
}

# Poisson Real Data (Wytham Great Tit data) -------------------------------
# read in data
pois_real <- read.csv('data/Wytham_tits.csv')

## make missing data types
## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness
  tempOutDf <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, 
                                         typeMissing = "random", 
                                         autoCorr = inputAutocor[i]))
  
  # name the elements of the Df with the amount of missingness 
  names(tempOutDf) <- paste0("Broods_",names(tempOutDf))
  tempOutDf <- cbind(pois_real, tempOutDf)
  
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("pois_real_randMiss_autoCorr_0") , 
           value = tempOutDf)
  } else {
    assign(x = paste0("pois_real_randMiss_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutDf)
  }
}

# missing in min and max of data
pois_real_minMaxMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, 
                                                        typeMissing = "minMax"))


names(pois_real_minMaxMiss_TEMP) <- paste0("Broods_",names(pois_real_minMaxMiss_TEMP))

pois_real_minMaxMiss <- cbind(pois_real, pois_real_minMaxMiss_TEMP)


# Store missing data  -----------------------------------------------------
# all datasets will be stored in "data/missingDatasets/"

# if it doesn't exist, make a folder to hold the datasets
if (dir.exists("./data/missingDatasets") == FALSE) {
  dir.create("./data/missingDatasets")
}

## store simulated Gaussian data (are stored in a list, each elemnt of the list 
# is a simulation run. Within each simulation run, the $y element contains 16 
# elements that have the response variable ranging from no missing data to the 
# highest proportion of missing data. The $sim_params element contains the 
# parameters used to generate that simulated dataset)
saveRDS(gauss_sim_randMiss_autoCorr_0, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_0.rds")
saveRDS(gauss_sim_randMiss_autoCorr_10, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_10.rds")
saveRDS(gauss_sim_randMiss_autoCorr_20, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_20.rds")
saveRDS(gauss_sim_randMiss_autoCorr_30, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_30.rds")
saveRDS(gauss_sim_randMiss_autoCorr_40, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_40.rds")
saveRDS(gauss_sim_randMiss_autoCorr_50, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_50.rds")
saveRDS(gauss_sim_randMiss_autoCorr_60, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_60.rds")
saveRDS(gauss_sim_randMiss_autoCorr_70, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_70.rds")
saveRDS(gauss_sim_randMiss_autoCorr_80, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_80.rds")
saveRDS(gauss_sim_randMiss_autoCorr_90, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_90.rds")
saveRDS(gauss_sim_minMaxMiss, file = "./data/missingDatasets/gauss_sim_minMaxMiss.rds")


## store simulated Poisson data (are stored in a list, each elemnt of the list 
# is a simulation run. Within each simulation run, the $y element contains 16 
# elements that have the response variable ranging from no missing data to the 
# highest proportion of missing data. The $sim_params element contains the 
# parameters used to generate that simulated dataset)
saveRDS(pois_sim_randMiss_autoCorr_0, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_0.rds")
saveRDS(pois_sim_randMiss_autoCorr_10, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_10.rds")
saveRDS(pois_sim_randMiss_autoCorr_20, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_20.rds")
saveRDS(pois_sim_randMiss_autoCorr_30, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_30.rds")
saveRDS(pois_sim_randMiss_autoCorr_40, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_40.rds")
saveRDS(pois_sim_randMiss_autoCorr_50, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_50.rds")
saveRDS(pois_sim_randMiss_autoCorr_60, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_60.rds")
saveRDS(pois_sim_randMiss_autoCorr_70, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_70.rds")
saveRDS(pois_sim_randMiss_autoCorr_80, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_80.rds")
saveRDS(pois_sim_randMiss_autoCorr_90, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_90.rds")
saveRDS(pois_sim_minMaxMiss, file = "./data/missingDatasets/pois_sim_minMaxMiss.rds")


## store real Gaussian data (a data frame with columns added for increasing 
# amounts of missingness in the response variable)
saveRDS(gauss_real_randMiss_autoCorr_0, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_0.rds")
saveRDS(gauss_real_randMiss_autoCorr_10, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_10.rds")
saveRDS(gauss_real_randMiss_autoCorr_20, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_20.rds")
saveRDS(gauss_real_randMiss_autoCorr_30, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_30.rds")
saveRDS(gauss_real_randMiss_autoCorr_40, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_40.rds")
saveRDS(gauss_real_randMiss_autoCorr_50, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_50.rds")
saveRDS(gauss_real_randMiss_autoCorr_60, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_60.rds")
saveRDS(gauss_real_randMiss_autoCorr_70, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_70.rds")
saveRDS(gauss_real_randMiss_autoCorr_80, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_80.rds")
saveRDS(gauss_real_randMiss_autoCorr_90, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_90.rds")
saveRDS(gauss_real_minMaxMiss, file = "./data/missingDatasets/gauss_real_minMaxMiss.rds")

## store real Poisson data (a data frame with columns added for increasing 
# amounts of missingness in the response variable)
saveRDS(pois_real_randMiss_autoCorr_0, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_0.rds")
saveRDS(pois_real_randMiss_autoCorr_10, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_10.rds")
saveRDS(pois_real_randMiss_autoCorr_20, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_20.rds")
saveRDS(pois_real_randMiss_autoCorr_30, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_30.rds")
saveRDS(pois_real_randMiss_autoCorr_40, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_40.rds")
saveRDS(pois_real_randMiss_autoCorr_50, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_50.rds")
saveRDS(pois_real_randMiss_autoCorr_60, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_60.rds")
saveRDS(pois_real_randMiss_autoCorr_70, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_70.rds")
saveRDS(pois_real_randMiss_autoCorr_80, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_80.rds")
saveRDS(pois_real_randMiss_autoCorr_90, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_90.rds")
saveRDS(pois_real_minMaxMiss, file = "./data/missingDatasets/pois_real_minMaxMiss.rds")


# prepare datastes for Beartooth runs -------------------------------------

## bind missing datasets (Gauss sim) into sets of 5000 nested lists, rather than 1000
# for each "gauss_sim" list, name the sublist for each simulation
names(gauss_sim_randMiss_autoCorr_0) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_0")
names(gauss_sim_randMiss_autoCorr_10) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_10")
names(gauss_sim_randMiss_autoCorr_20) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_20")
names(gauss_sim_randMiss_autoCorr_30) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_30")
names(gauss_sim_randMiss_autoCorr_40) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_40")
names(gauss_sim_randMiss_autoCorr_50) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_50")
names(gauss_sim_randMiss_autoCorr_60) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_60")
names(gauss_sim_randMiss_autoCorr_70) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_70")
names(gauss_sim_randMiss_autoCorr_80) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_80")
names(gauss_sim_randMiss_autoCorr_90) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_90")

gauss_sim_randMiss_A <- c(gauss_sim_randMiss_autoCorr_0, 
                             gauss_sim_randMiss_autoCorr_10,
                             gauss_sim_randMiss_autoCorr_20,
                             gauss_sim_randMiss_autoCorr_30,
                             gauss_sim_randMiss_autoCorr_40
                             )
gauss_sim_randMiss_B <- c(gauss_sim_randMiss_autoCorr_50,
                          gauss_sim_randMiss_autoCorr_60,
                          gauss_sim_randMiss_autoCorr_70,
                          gauss_sim_randMiss_autoCorr_80,
                          gauss_sim_randMiss_autoCorr_90)

# Remove the y_no miss from all nested lists (Amelia doesn't like NO missing values) 
# for first chunk of data
test <- gauss_sim_randMiss_A

for (i in 1:length(gauss_sim_randMiss_A)) {
  # remove the "y_noMiss" sub_list
  test[[i]] <- gauss_sim_randMiss_A[[i]]$y[2:16]
  
}

gauss_sim_randMiss_A <- test

# for second chunk of data
test <- gauss_sim_randMiss_B

for (i in 1:length(gauss_sim_randMiss_B)) {
  # remove the "y_noMiss" sub_list
  test[[i]] <- gauss_sim_randMiss_B[[i]]$y[2:16]
  
}

gauss_sim_randMiss_B <- test

# pull out param values from the simulation (sim_pars) and stick them in the identifier 
# along w/ missing prop and autocor? (If that's a pain, just a separate file that has 
# all 5000 simulations and their associated parameters would be helpful. We should only 
# have 1000 unique sets of simulation parameters) 

gauss_sim_params <- data.frame(
  "SimNumber" = vector(mode = "integer", length = 1000),
  "phi" = vector(mode = "double", length = 1000),
  "beta1" = vector(mode = "double", length = 1000), 
  "beta2" = vector(mode = "double", length = 1000),
  "beta3" = vector(mode = "double", length = 1000)
)

for (i in 1:length(gauss_sim)) {
  gauss_sim_params[i, "SimNumber"] <- i
  gauss_sim_params[i,"phi"] <- gauss_sim[[i]]$sim_params$phi
  gauss_sim_params[i,"beta1"] <- gauss_sim[[i]]$sim_params$beta[1]
  gauss_sim_params[i,"beta2"] <- gauss_sim[[i]]$sim_params$beta[2]
  gauss_sim_params[i,"beta3"] <- gauss_sim[[i]]$sim_params$beta[3]
}

# save missing datasets for Beartooth -------------------------------------
# if it doesn't exist, make a folder to hold the datasets
if (dir.exists("./data/missingDatasets/forBeartooth") == FALSE) {
  dir.create("./data/missingDatasets/forBeartooth")
}

## store simulated Gaussian data (are stored in a list, each elemnt of the list 
# is a simulation run. Within each simulation run, the $y element contains 15 
# elements that have the response variable ranging from the lowest amount of 
# missing data to the highest proportion of missing data. 
saveRDS(gauss_sim_randMiss_A, file = "./data/missingDatasets/forBeartooth/gauss_sim_randMiss_A.rds")
saveRDS(gauss_sim_randMiss_B, file = "./data/missingDatasets/forBeartooth/gauss_sim_randMiss_B.rds")

## save a data.frame that has the parameters used to run each simulation (1000 
# rows, each corresponding to a simulation run)
saveRDS(gauss_sim_params, file = "./data/missingDatasets/forBeartooth/gauss_sim_params.rds")
