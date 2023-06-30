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

# make missing data types
# missing at random (autocorrelation = 0.01)
gauss_sim_randMiss_autoCorr_01 <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = 0.01), 

       "sim_params" <- x$sim_params)
)


for (i in 1:length(gauss_sim_randMiss_autoCorr_01)) {
  names(gauss_sim_randMiss_autoCorr_01[[i]]) <- c("y", "sim_params")
  gauss_sim_randMiss_autoCorr_01[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_randMiss_autoCorr_01[[i]]$y)
}

# missing at random (medium autocorrelation = .25)
gauss_sim_randMiss_autoCorr_25 <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .25), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_randMiss_autoCorr_25)) {
  names(gauss_sim_randMiss_autoCorr_25[[i]]) <- c("y", "sim_params")
  gauss_sim_randMiss_autoCorr_25[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_randMiss_autoCorr_25[[i]]$y)
}

# missing at random (medium autocorrelation = .50)
gauss_sim_randMiss_autoCorr_50 <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .50), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_randMiss_autoCorr_50)) {
  names(gauss_sim_randMiss_autoCorr_50[[i]]) <- c("y", "sim_params")
  gauss_sim_randMiss_autoCorr_50[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_randMiss_autoCorr_50[[i]]$y)
}

# missing at random (medium autocorrelation = .75)
gauss_sim_randMiss_autoCorr_75 <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .75), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_randMiss_autoCorr_75)) {
  names(gauss_sim_randMiss_autoCorr_75[[i]]) <- c("y", "sim_params")
  gauss_sim_randMiss_autoCorr_75[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_randMiss_autoCorr_75[[i]]$y)
}

# missing at random (medium autocorrelation = .90)
gauss_sim_randMiss_autoCorr_90 <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .90), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_randMiss_autoCorr_90)) {
  names(gauss_sim_randMiss_autoCorr_90[[i]]) <- c("y", "sim_params")
  gauss_sim_randMiss_autoCorr_90[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_randMiss_autoCorr_90[[i]]$y)
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
gauss_real <-as.data.frame(pr)

# missing at random (autocorrelation = .01)
gauss_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, typeMissing = "random", autoCorr = .01))
names(gauss_real_randMiss_TEMP) <- paste0("GPP_",names(gauss_real_randMiss_TEMP))

gauss_real_randMiss_autoCorr_01 <- cbind(gauss_real, gauss_real_randMiss_TEMP)

# missing at random (autocorrelation = .25)
gauss_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, typeMissing = "random", autoCorr = .25))
names(gauss_real_randMiss_TEMP) <- paste0("GPP_",names(gauss_real_randMiss_TEMP))

gauss_real_randMiss_autoCorr_25 <- cbind(gauss_real, gauss_real_randMiss_TEMP)

# missing at random (autocorrelation = .50)
gauss_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, typeMissing = "random", autoCorr = .50))
names(gauss_real_randMiss_TEMP) <- paste0("GPP_",names(gauss_real_randMiss_TEMP))

gauss_real_randMiss_autoCorr_50 <- cbind(gauss_real, gauss_real_randMiss_TEMP)


# missing at random (autocorrelation = .75)
gauss_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, typeMissing = "random", autoCorr = .75))
names(gauss_real_randMiss_TEMP) <- paste0("GPP_",names(gauss_real_randMiss_TEMP))

gauss_real_randMiss_autoCorr_75 <- cbind(gauss_real, gauss_real_randMiss_TEMP)

# missing at random (autocorrelation = .90)
gauss_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, typeMissing = "random", autoCorr = .90))
names(gauss_real_randMiss_TEMP) <- paste0("GPP_",names(gauss_real_randMiss_TEMP))

gauss_real_randMiss_autoCorr_90 <- cbind(gauss_real, gauss_real_randMiss_TEMP)

# missing in min and max of data
gauss_real_minMaxMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, 
                                                   typeMissing = "minMax"))


names(gauss_real_minMaxMiss_TEMP) <- paste0("GPP_",names(gauss_real_minMaxMiss_TEMP))

gauss_real_minMaxMiss <- cbind(gauss_real, gauss_real_minMaxMiss_TEMP)


# Poisson Simulated Data --------------------------------------------------
# read in data
pois_sim <- readRDS("./data/ricker_0miss_datasets.rds")

# make missing data types
# missing at random (autocorrelation = .01)
pois_sim_randMiss_autoCorr_01 <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .01), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_randMiss_autoCorr_01)) {
  names(pois_sim_randMiss_autoCorr_01[[i]]) <- c("y", "sim_params")
  pois_sim_randMiss_autoCorr_01[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_randMiss_autoCorr_01[[i]]$y)
}

# missing at random (autocorrelation = .25)
pois_sim_randMiss_autoCorr_25 <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .25), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_randMiss_autoCorr_25)) {
  names(pois_sim_randMiss_autoCorr_25[[i]]) <- c("y", "sim_params")
  pois_sim_randMiss_autoCorr_25[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_randMiss_autoCorr_25[[i]]$y)
}

# missing at random (autocorrelation = .50)
pois_sim_randMiss_autoCorr_50 <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .50), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_randMiss_autoCorr_50)) {
  names(pois_sim_randMiss_autoCorr_50[[i]]) <- c("y", "sim_params")
  pois_sim_randMiss_autoCorr_50[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_randMiss_autoCorr_50[[i]]$y)
}

# missing at random (autocorrelation = .75)
pois_sim_randMiss_autoCorr_75 <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .75), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_randMiss_autoCorr_75)) {
  names(pois_sim_randMiss_autoCorr_75[[i]]) <- c("y", "sim_params")
  pois_sim_randMiss_autoCorr_75[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_randMiss_autoCorr_75[[i]]$y)
}

# missing at random (autocorrelation = .90)
pois_sim_randMiss_autoCorr_90 <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random", autoCorr = .90), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_randMiss_autoCorr_90)) {
  names(pois_sim_randMiss_autoCorr_90[[i]]) <- c("y", "sim_params")
  pois_sim_randMiss_autoCorr_90[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_randMiss_autoCorr_90[[i]]$y)
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
# missing at random (autocorrelation = .01)
pois_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, typeMissing = "random", autoCorr = .01))
names(pois_real_randMiss_TEMP) <- paste0("Broods_",names(pois_real_randMiss_TEMP))

pois_real_randMiss_autoCorr_01 <- cbind(pois_real, pois_real_randMiss_TEMP)

# missing at random (autocorrelation = .25)
pois_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, typeMissing = "random", autoCorr = .25))
names(pois_real_randMiss_TEMP) <- paste0("Broods_",names(pois_real_randMiss_TEMP))

pois_real_randMiss_autoCorr_25 <- cbind(pois_real, pois_real_randMiss_TEMP)

# missing at random (autocorrelation = .50)
pois_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, typeMissing = "random", autoCorr = .50))
names(pois_real_randMiss_TEMP) <- paste0("Broods_",names(pois_real_randMiss_TEMP))

pois_real_randMiss_autoCorr_50 <- cbind(pois_real, pois_real_randMiss_TEMP)

# missing at random (autocorrelation = .75)
pois_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, typeMissing = "random", autoCorr = .75))
names(pois_real_randMiss_TEMP) <- paste0("Broods_",names(pois_real_randMiss_TEMP))

pois_real_randMiss_autoCorr_75 <- cbind(pois_real, pois_real_randMiss_TEMP)

# missing at random (autocorrelation = .90)
pois_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, typeMissing = "random", autoCorr = .90))
names(pois_real_randMiss_TEMP) <- paste0("Broods_",names(pois_real_randMiss_TEMP))

pois_real_randMiss_autoCorr_90 <- cbind(pois_real, pois_real_randMiss_TEMP)

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
saveRDS(gauss_sim_randMiss_autoCorr_01, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_01.rds")
saveRDS(gauss_sim_randMiss_autoCorr_25, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_25.rds")
saveRDS(gauss_sim_randMiss_autoCorr_50, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_50.rds")
saveRDS(gauss_sim_randMiss_autoCorr_75, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_75.rds")
saveRDS(gauss_sim_randMiss_autoCorr_90, file = "./data/missingDatasets/gauss_sim_randMiss_autoCorr_90.rds")
saveRDS(gauss_sim_minMaxMiss, file = "./data/missingDatasets/gauss_sim_minMaxMiss.rds")


## store simulated Poisson data (are stored in a list, each elemnt of the list 
# is a simulation run. Within each simulation run, the $y element contains 16 
# elements that have the response variable ranging from no missing data to the 
# highest proportion of missing data. The $sim_params element contains the 
# parameters used to generate that simulated dataset)
saveRDS(pois_sim_randMiss_autoCorr_01, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_01.rds")
saveRDS(pois_sim_randMiss_autoCorr_25, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_25.rds")
saveRDS(pois_sim_randMiss_autoCorr_50, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_50.rds")
saveRDS(pois_sim_randMiss_autoCorr_75, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_75.rds")
saveRDS(pois_sim_randMiss_autoCorr_90, file = "./data/missingDatasets/pois_sim_randMiss_autoCorr_90.rds")
saveRDS(pois_sim_minMaxMiss, file = "./data/missingDatasets/pois_sim_minMaxMiss.rds")


## store real Gaussian data (a data frame with columns added for increasing 
# amounts of missingness in the response variable)
saveRDS(gauss_real_randMiss_autoCorr_01, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_01.rds")
saveRDS(gauss_real_randMiss_autoCorr_25, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_25.rds")
saveRDS(gauss_real_randMiss_autoCorr_50, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_50.rds")
saveRDS(gauss_real_randMiss_autoCorr_75, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_75.rds")
saveRDS(gauss_real_randMiss_autoCorr_90, file = "./data/missingDatasets/gauss_real_randMiss_autoCorr_90.rds")
saveRDS(gauss_real_minMaxMiss, file = "./data/missingDatasets/gauss_real_minMaxMiss.rds")

## store real Poisson data (a data frame with columns added for increasing 
# amounts of missingness in the response variable)
saveRDS(pois_real_randMiss_autoCorr_01, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_01.rds")
saveRDS(pois_real_randMiss_autoCorr_25, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_25.rds")
saveRDS(pois_real_randMiss_autoCorr_50, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_50.rds")
saveRDS(pois_real_randMiss_autoCorr_75, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_75.rds")
saveRDS(pois_real_randMiss_autoCorr_90, file = "./data/missingDatasets/pois_real_randMiss_autoCorr_90.rds")
saveRDS(pois_real_minMaxMiss, file = "./data/missingDatasets/pois_real_minMaxMiss.rds")

