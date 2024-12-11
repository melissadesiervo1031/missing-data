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

# Gaussian Real Data (Au Sable River) --------------------------------
# read in data
gauss_auSable <- read.csv("./data/au_sable_river_prepped.csv")
## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
#inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)
inputAutocor <- c(0, .25, .25, .25, .5, .5, .5, .75, .75, .75)
inputPropMiss <- c(.2, .2, .2, .2, .2, .4, .4, .4, .4, 4., .6, .6, .6, .6, .6)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness
  tempOutDf <- as.data.frame(makeMissing(timeSeries = gauss_auSable$GPP, 
                                         typeMissing = "random", 
                                         autoCorr = inputAutocor[i], 
                                         propMiss = inputPropMiss))
  
  # name the elements of the Df with the amount of missingness 
  names(tempOutDf) <- paste0("GPP_",names(tempOutDf))
  tempOutDf <- cbind(gauss_auSable, tempOutDf)
  
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("gauss_auSable_randMiss_autoCorr_0") , 
           value = tempOutDf)
  } else {
    assign(x = paste0("gauss_auSable_randMiss_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutDf)
    
  }
}

# put all of the data into one list
autcorVector <- c(0, 25, 25, 25, 50, 50, 50, 75, 75, 75)
gauss_auSable_randMiss_list <- vector(mode = "list", length = 10) 
names(gauss_auSable_randMiss_list) <- paste0("gauss_auSable_randMiss_autoCor_", autcorVector)
for (i in 1:length(inputAutocor)) {
  # get the correct autocorrelation/missing data data.frame
  tempDF <- get(x = paste0("gauss_auSable_randMiss_autoCorr_", autcorVector[i]))
  # change into a list element
  gauss_auSable_randMiss_list[[i]]$y <- tempDF %>% 
    rename_with(~ str_replace(string = names(tempDF), pattern = "GPP_", replacement = "")) %>% 
    rename("y" = "GPP" ) %>% 
    select(-date, -ER, -light, -Q) %>% 
    as.list()
  gauss_auSable_randMiss_list[[i]]$sim_params$date <- tempDF %>% 
    select(date)
  gauss_auSable_randMiss_list[[i]]$sim_params$light <- tempDF %>% 
    select(light)
  gauss_auSable_randMiss_list[[i]]$sim_params$Q <- tempDF %>% 
    select(Q)
}

gauss_auSable_randMiss <- gauss_auSable_randMiss_list

## missing in min and max of data
## for the MNAR data, we need to do slightly different propMissing values, since the same proportion will cutoff the same exact values 

inputPropMiss <- c(.18, .19, .2, .21, .22, .38, .39, .4, .41, 42., .58, .59, .6, .61, .62)
gauss_auSable_minMaxMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_auSable$GPP, 
                                                           typeMissing = "minMax", 
                                                           propMiss = inputPropMiss))


#names(gauss_auSable_minMaxMiss_TEMP) <- paste0("GPP_",names(gauss_auSable_minMaxMiss_TEMP))

gauss_auSable_minMaxMiss <- cbind(gauss_auSable, gauss_auSable_minMaxMiss_TEMP)

# change "broods" to "y" in the data.frame
names(gauss_auSable_minMaxMiss)[2] <- "y" 

# transform into a list (for consistency w/ simulated data)
gauss_auSable_minMaxMiss_list <- vector(mode = "list", length = 1) 
gauss_auSable_minMaxMiss_list[[1]]$y <- as.list(gauss_auSable_minMaxMiss[c(2,9:ncol(gauss_auSable_minMaxMiss))])
gauss_auSable_minMaxMiss_list[[1]]$sim_params$date <- gauss_auSable_minMaxMiss$date 
gauss_auSable_minMaxMiss_list[[1]]$sim_params$light <- gauss_auSable_minMaxMiss$light
gauss_auSable_minMaxMiss_list[[1]]$sim_params$light.rel <- gauss_auSable_minMaxMiss$light.rel
gauss_auSable_minMaxMiss_list[[1]]$sim_params$Q <- gauss_auSable_minMaxMiss$Q

gauss_auSable_minMaxMiss <- gauss_auSable_minMaxMiss_list

# Gaussian Real Data (Badger Mill Creek) --------------------------------
# read in data
gauss_badger <- read.csv("./data/badger_mill_creek_prepped.csv")
## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
#inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)
inputAutocor <- c(0, .25, .25, .25, .5, .5, .5, .75, .75, .75)
inputPropMiss <- c(.2, .2, .2, .2, .2, .4, .4, .4, .4, 4., .6, .6, .6, .6, .6)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness
  tempOutDf <- as.data.frame(makeMissing(timeSeries = gauss_badger$GPP, 
                                         typeMissing = "random", 
                                         autoCorr = inputAutocor[i], 
                                         propMiss = inputPropMiss))
  
  # name the elements of the Df with the amount of missingness 
  names(tempOutDf) <- paste0("GPP_",names(tempOutDf))
  tempOutDf <- cbind(gauss_badger, tempOutDf)
  
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("gauss_badger_randMiss_autoCorr_0") , 
           value = tempOutDf)
  } else {
    assign(x = paste0("gauss_badger_randMiss_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutDf)
    
  }
}

# put all of the data into one list
autcorVector <- c(0, 25, 25, 25, 50, 50, 50, 75, 75, 75)
gauss_badger_randMiss_list <- vector(mode = "list", length = 10) 
names(gauss_badger_randMiss_list) <- paste0("gauss_badger_randMiss_autoCor_", autcorVector)
for (i in 1:length(inputAutocor)) {
  # get the correct autocorrelation/missing data data.frame
  tempDF <- get(x = paste0("gauss_badger_randMiss_autoCorr_", autcorVector[i]))
  # change into a list element
  gauss_badger_randMiss_list[[i]]$y <- tempDF %>% 
    rename_with(~ str_replace(string = names(tempDF), pattern = "GPP_", replacement = "")) %>% 
    rename("y" = "GPP" ) %>% 
    select(-date, -ER, -light, -Q) %>% 
    as.list()
  gauss_badger_randMiss_list[[i]]$sim_params$date <- tempDF %>% 
    select(date)
  gauss_badger_randMiss_list[[i]]$sim_params$light <- tempDF %>% 
    select(light)
  gauss_badger_randMiss_list[[i]]$sim_params$Q <- tempDF %>% 
    select(Q)
}

gauss_badger_randMiss <- gauss_badger_randMiss_list

## missing in min and max of data
## for the MNAR data, we need to do slightly different propMissing values, since the same proportion will cutoff the same exact values 

inputPropMiss <- c(.18, .19, .2, .21, .22, .38, .39, .4, .41, 42., .58, .59, .6, .61, .62)
gauss_badger_minMaxMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_badger$GPP, 
                                                           typeMissing = "minMax", 
                                                           propMiss = inputPropMiss))
#names(gauss_badger_minMaxMiss_TEMP) <- paste0("GPP_",names(gauss_badger_minMaxMiss_TEMP))

gauss_badger_minMaxMiss <- cbind(gauss_badger, gauss_badger_minMaxMiss_TEMP)

# change "broods" to "y" in the data.frame
names(gauss_badger_minMaxMiss)[2] <- "y" 

# transform into a list (for consistency w/ simulated data)
gauss_badger_minMaxMiss_list <- vector(mode = "list", length = 1) 
gauss_badger_minMaxMiss_list[[1]]$y <- as.list(gauss_badger_minMaxMiss[c(2,9:ncol(gauss_badger_minMaxMiss))])
gauss_badger_minMaxMiss_list[[1]]$sim_params$date <- gauss_badger_minMaxMiss$date 
gauss_badger_minMaxMiss_list[[1]]$sim_params$light <- gauss_badger_minMaxMiss$light
gauss_badger_minMaxMiss_list[[1]]$sim_params$light.rel <- gauss_badger_minMaxMiss$light.rel
gauss_badger_minMaxMiss_list[[1]]$sim_params$Q <- gauss_badger_minMaxMiss$Q

gauss_badger_minMaxMiss <- gauss_badger_minMaxMiss_list

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

# put all of the data into one list
autcorVector <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)
pois_real_randMiss_list <- vector(mode = "list", length = 10) 
names(pois_real_randMiss_list) <- paste0("pois_real_randMiss_autoCor_", autcorVector)
for (i in 1:length(inputAutocor)) {
  # get the correct autocorrelation/missing data data.frame
  tempDF <- get(x = paste0("pois_real_randMiss_autoCorr_", autcorVector[i]))
  # change into a list element
  pois_real_randMiss_list[[i]]$y <- tempDF %>% 
    rename_with(~ str_replace(string = names(tempDF), pattern = "Broods_", replacement = "")) %>% 
    rename("y" = "Broods" ) %>% 
    select(-Year) %>% 
    as.list()
  pois_real_randMiss_list[[i]]$sim_params <- NA
}
pois_real_randMiss <- pois_real_randMiss_list

## missing in min and max of data
pois_real_minMaxMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, 
                                                       typeMissing = "minMax"))


names(pois_real_minMaxMiss_TEMP) <- paste0("Broods_",names(pois_real_minMaxMiss_TEMP))

pois_real_minMaxMiss <- cbind(pois_real, pois_real_minMaxMiss_TEMP)

# change "broods" to "y" in the data.frame
names(pois_real_minMaxMiss) <- str_replace(string = names(pois_real_minMaxMiss), pattern = "Broods_", replacement = "")
names(pois_real_minMaxMiss)[2] <- "y" 

# transform into a list (for consistency w/ simulated data)
pois_real_minMaxMiss_list <- vector(mode = "list", length = 1) 
pois_real_minMaxMiss_list[[1]]$y <- as.list(pois_real_minMaxMiss[2:ncol(pois_real_minMaxMiss)])
pois_real_minMaxMiss_list[[1]]$sim_params <- NA

pois_real_minMaxMiss <- pois_real_minMaxMiss_list

# Store missing data  -----------------------------------------------------
# all datasets will be stored in "data/missingDatasets/"

# if it doesn't exist, make a folder to hold the datasets
if (dir.exists("./data/missingDatasets") == FALSE) {
  dir.create("./data/missingDatasets")
}

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
#test <- gauss_sim_randMiss_A
test <- vector(mode = "list", length = 5000)
for (i in 1:length(gauss_sim_randMiss_A)) {
  # remove the "y_noMiss" sub_list
  test[[i]]$y <- gauss_sim_randMiss_A[[i]]$y[2:16]
  test[[i]]$sim_params$X <- gauss_sim_randMiss_A[[i]]$sim$X
}
# rename w/ accurate names
names(test) <- names(gauss_sim_randMiss_A)
gauss_sim_randMiss_A <- test

# for second chunk of data
test <- gauss_sim_randMiss_B

for (i in 1:length(gauss_sim_randMiss_B)) {
  # remove the "y_noMiss" sub_list
  test[[i]]$y <- gauss_sim_randMiss_B[[i]]$y[2:16]
  test[[i]]$sim_params$X <- gauss_sim_randMiss_B[[i]]$sim$X
}
# rename w/ accurate names
names(test) <- names(gauss_sim_randMiss_B)
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


## do the same for the Ricker data

## bind missing datasets (pois sim) into sets of 5000 nested lists, rather than 1000
# for each "pois_sim" list, name the sublist for each simulation
names(pois_sim_randMiss_autoCorr_0) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_0")
names(pois_sim_randMiss_autoCorr_10) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_10")
names(pois_sim_randMiss_autoCorr_20) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_20")
names(pois_sim_randMiss_autoCorr_30) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_30")
names(pois_sim_randMiss_autoCorr_40) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_40")
names(pois_sim_randMiss_autoCorr_50) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_50")
names(pois_sim_randMiss_autoCorr_60) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_60")
names(pois_sim_randMiss_autoCorr_70) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_70")
names(pois_sim_randMiss_autoCorr_80) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_80")
names(pois_sim_randMiss_autoCorr_90) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_90")

pois_sim_randMiss_A <- c(pois_sim_randMiss_autoCorr_0, 
                         pois_sim_randMiss_autoCorr_10,
                         pois_sim_randMiss_autoCorr_20,
                         pois_sim_randMiss_autoCorr_30,
                         pois_sim_randMiss_autoCorr_40
)
pois_sim_randMiss_B <- c(pois_sim_randMiss_autoCorr_50,
                         pois_sim_randMiss_autoCorr_60,
                         pois_sim_randMiss_autoCorr_70,
                         pois_sim_randMiss_autoCorr_80,
                         pois_sim_randMiss_autoCorr_90)

# Remove the y_no miss from all nested lists (Amelia doesn't like NO missing values) 
# for first chunk of data
#test <- pois_sim_randMiss_A
test <- vector(mode = "list", length = 5000)
for (i in 1:length(pois_sim_randMiss_A)) {
  # remove the "y_noMiss" sub_list
  test[[i]]$y <- pois_sim_randMiss_A[[i]]$y[2:16]
  test[[i]]$sim_params$X <- pois_sim_randMiss_A[[i]]$sim$X
}
# rename w/ accurate names
names(test) <- names(pois_sim_randMiss_A)
pois_sim_randMiss_A <- test

# for second chunk of data
test <- pois_sim_randMiss_B

for (i in 1:length(pois_sim_randMiss_B)) {
  # remove the "y_noMiss" sub_list
  test[[i]]$y <- pois_sim_randMiss_B[[i]]$y[2:16]
  test[[i]]$sim_params$X <- pois_sim_randMiss_B[[i]]$sim$X
}
# rename w/ accurate names
names(test) <- names(pois_sim_randMiss_B)
pois_sim_randMiss_B <- test

# pull out param values from the simulation (sim_pars) and stick them in the identifier 
# along w/ missing prop and autocor? (If that's a pain, just a separate file that has 
# all 5000 simulations and their associated parameters would be helpful. We should only 
# have 1000 unique sets of simulation parameters) 

pois_sim_params <- data.frame(
  "SimNumber" = vector(mode = "integer", length = 1000),
  "r" = vector(mode = "double", length = 1000),
  "alpha" = vector(mode = "double", length = 1000), 
  "N0" = vector(mode = "double", length = 1000)
)

for (i in 1:length(pois_sim)) {
  pois_sim_params[i, "SimNumber"] <- i
  pois_sim_params[i,"r"] <- pois_sim[[i]]$sim_params$r
  pois_sim_params[i,"alpha"] <- pois_sim[[i]]$sim_params$alpha
  pois_sim_params[i,"N0"] <- pois_sim[[i]]$sim_params$N0
}

# save missing datasets for Beartooth -------------------------------------
## store simulated Poisson data (are stored in a list, each element of the list 
# is a simulation run. Within each simulation run, the $y element contains 15 
# elements that have the response variable ranging from the lowest amount of 
# missing data to the highest proportion of missing data. 
saveRDS(pois_sim_randMiss_A, file = "./data/missingDatasets/pois_sim_randMiss_A.rds")
saveRDS(pois_sim_randMiss_B, file = "./data/missingDatasets/pois_sim_randMiss_B.rds")

## save a data.frame that has the parameters used to run each simulation (1000 
# rows, each corresponding to a simulation run)
saveRDS(pois_sim_params, file = "./data/missingDatasets/pois_sim_params.rds")

## save poisson simulated MNAR data
saveRDS(pois_sim_minMaxMiss, file = "./data/missingDatasets/pois_sim_MinMaxMiss.rds")

## save poisson real MNAR data
saveRDS(pois_real_minMaxMiss, file = "./data/missingDatasets/pois_real_MinMaxMiss.rds")
## save poisson real MAR data
saveRDS(pois_real_randMiss, file = "./data/missingDatasets/pois_real_randMiss.rds")



## store simulated Gaussian data (are stored in a list, each element of the list 
# is a simulation run. Within each simulation run, the $y element contains 15 
# elements that have the response variable ranging from the lowest amount of 
# missing data to the highest proportion of missing data. 
saveRDS(gauss_sim_randMiss_A, file = "./data/missingDatasets/gauss_sim_randMiss_A.rds")
saveRDS(gauss_sim_randMiss_B, file = "./data/missingDatasets/gauss_sim_randMiss_B.rds")

## save a data.frame that has the parameters used to run each simulation (1000 
# rows, each corresponding to a simulation run)
saveRDS(gauss_sim_params, file = "./data/missingDatasets/gauss_sim_params.rds")

## save gaussian simulated MNAR data
saveRDS(gauss_sim_minMaxMiss, file = "./data/missingDatasets/gauss_sim_MinMaxMiss.rds")

## save gaussian real MNAR data
saveRDS(gauss_auSable_minMaxMiss, file = "./data/missingDatasets/gauss_real_auSable_MinMaxMiss.rds")
## save gaussian real MAR data
saveRDS(gauss_auSable_randMiss, file = "./data/missingDatasets/gauss_real_auSable_randMiss.rds")

## save gaussian real MNAR data
saveRDS(gauss_badger_minMaxMiss, file = "./data/missingDatasets/gauss_real_badger_MinMaxMiss.rds")
## save gaussian real MAR data
saveRDS(gauss_badger_randMiss, file = "./data/missingDatasets/gauss_real_badger_randMiss.rds")
