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
source("./missing_data_functions.R")

### Note: throughout the simulations of added missingness, only add missingness to response variable (i.e. GPP)

# Gaussian Simulated Data -------------------------------------------------
# read in data
gauss_sim <- readRDS("./data/gauss_ar1_0miss_datasets.rds")

# make missing data types
# missing at random
gauss_sim_randMiss <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random"), 

       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_randMiss)) {
  names(gauss_sim_randMiss[[i]]) <- c("y", "sim_params")
  gauss_sim_randMiss[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_randMiss[[i]]$y)
}


# missing in evenly spaced chunks
gauss_sim_evenChunkMiss <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "evenChunks", chunkSize = 3), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_evenChunkMiss)) {
  names(gauss_sim_evenChunkMiss[[i]]) <- c("y", "sim_params")
  gauss_sim_evenChunkMiss[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_evenChunkMiss[[i]]$y)
}

# missing in randomly spaced chunks
gauss_sim_randChunkMiss <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "randChunks", chunkSize = 4), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_randChunkMiss)) {
  names(gauss_sim_randMiss[[i]]) <- c("y", "sim_params")
  gauss_sim_randChunkMiss[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_randChunkMiss[[i]]$y)
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

# missing at random
gauss_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, typeMissing = "random"))
names(gauss_real_randMiss_TEMP) <- paste0("GPP_",names(gauss_real_randMiss_TEMP))

gauss_real_randMiss <- cbind(gauss_real, gauss_real_randMiss_TEMP)

# missing in evenly spaced chunks
gauss_real_evenChunkMiss_TEMP <-  as.data.frame(makeMissing(timeSeries = gauss_real$GPP, 
                                                      typeMissing = "evenChunks", chunkSize = 4))


names(gauss_real_evenChunkMiss_TEMP) <- paste0("GPP_",names(gauss_real_evenChunkMiss_TEMP))

gauss_real_evenChunkMiss <- cbind(gauss_real, gauss_real_evenChunkMiss_TEMP)


# missing in randomly spaced chunks
gauss_real_randChunkMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, 
                                                        typeMissing = "randChunks", chunkSize = 4))


names(gauss_real_randChunkMiss_TEMP) <- paste0("GPP_",names(gauss_real_randChunkMiss_TEMP))

gauss_real_randChunkMiss <- cbind(gauss_real, gauss_real_randChunkMiss_TEMP)


# missing in min and max of data
gauss_real_minMaxMiss_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_real$GPP, 
                                                   typeMissing = "minMax"))


names(gauss_real_minMaxMiss_TEMP) <- paste0("GPP_",names(gauss_real_minMaxMiss_TEMP))

gauss_real_minMaxMiss <- cbind(gauss_real, gauss_real_minMaxMiss_TEMP)


# Poisson Simulated Data --------------------------------------------------
# read in data
pois_sim <- readRDS("./data/ricker_0miss_datasets.rds")

# make missing data types
# missing at random
pois_sim_randMiss <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "random"), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_randMiss)) {
  names(pois_sim_randMiss[[i]]) <- c("y", "sim_params")
  pois_sim_randMiss[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_randMiss[[i]]$y)
}

# missing in evenly spaced chunks
pois_sim_evenChunkMiss <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "evenChunks", chunkSize = 3), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_evenChunkMiss)) {
  names(pois_sim_evenChunkMiss[[i]]) <- c("y", "sim_params")
  pois_sim_evenChunkMiss[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_evenChunkMiss[[i]]$y)
}

# missing in randomly spaced chunks
pois_sim_randChunkMiss <- lapply(X = pois_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y, typeMissing = "randChunks", chunkSize = 3), 
       "sim_params" <- x$sim_params)
)

for (i in 1:length(pois_sim_randChunkMiss)) {
  names(pois_sim_randMiss[[i]]) <- c("y", "sim_params")
  pois_sim_randChunkMiss[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_randChunkMiss[[i]]$y)
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
# missing at random
pois_real_randMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, typeMissing = "random"))
names(pois_real_randMiss_TEMP) <- paste0("Broods_",names(pois_real_randMiss_TEMP))

pois_real_randMiss <- cbind(pois_real, pois_real_randMiss_TEMP)

# missing in evenly spaced chunks
pois_real_evenChunkMiss_TEMP <-  as.data.frame(makeMissing(timeSeries = pois_real$Broods, 
                                                            typeMissing = "evenChunks", chunkSize = 2))


names(pois_real_evenChunkMiss_TEMP) <- paste0("Broods_",names(pois_real_evenChunkMiss_TEMP))

pois_real_evenChunkMiss <- cbind(pois_real, pois_real_evenChunkMiss_TEMP)


# missing in randomly spaced chunks
pois_real_randChunkMiss_TEMP <- as.data.frame(makeMissing(timeSeries = pois_real$Broods, 
                                                           typeMissing = "randChunks", chunkSize = 3))


names(pois_real_randChunkMiss_TEMP) <- paste0("Broods_",names(pois_real_randChunkMiss_TEMP))

pois_real_randChunkMiss <- cbind(pois_real, pois_real_randChunkMiss_TEMP)


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
saveRDS(gauss_sim_randMiss, file = "./data/missingDatasets/gauss_sim_randMiss.rds")
saveRDS(gauss_sim_evenChunkMiss, file = "./data/missingDatasets/gauss_sim_evenChunkMiss.rds")
saveRDS(gauss_sim_randChunkMiss, file = "./data/missingDatasets/gauss_sim_randChunkMiss.rds")
saveRDS(gauss_sim_minMaxMiss, file = "./data/missingDatasets/gauss_sim_minMaxMiss.rds")


## store simulated Poisson data (are stored in a list, each elemnt of the list 
# is a simulation run. Within each simulation run, the $y element contains 16 
# elements that have the response variable ranging from no missing data to the 
# highest proportion of missing data. The $sim_params element contains the 
# parameters used to generate that simulated dataset)
saveRDS(pois_sim_randMiss, file = "./data/missingDatasets/pois_sim_randMiss.rds")
saveRDS(pois_sim_evenChunkMiss, file = "./data/missingDatasets/pois_sim_evenChunkMiss.rds")
saveRDS(pois_sim_randChunkMiss, file = "./data/missingDatasets/pois_sim_randChunkMiss.rds")
saveRDS(pois_sim_minMaxMiss, file = "./data/missingDatasets/pois_sim_minMaxMiss.rds")


## store real Gaussian data (a data frame with columns added for increasing 
# amounts of missingness in the response variable)
saveRDS(gauss_real_randMiss, file = "./data/missingDatasets/gauss_real_randMiss.rds")
saveRDS(gauss_real_evenChunkMiss, file = "./data/missingDatasets/gauss_real_evenChunkMiss.rds")
saveRDS(gauss_real_randChunkMiss, file = "./data/missingDatasets/gauss_real_randChunkMiss.rds")
saveRDS(gauss_real_minMaxMiss, file = "./data/missingDatasets/gauss_real_minMaxMiss.rds")

## store real Poisson data (a data frame with columns added for increasing 
# amounts of missingness in the response variable)
saveRDS(pois_real_randMiss, file = "./data/missingDatasets/pois_real_randMiss.rds")
saveRDS(pois_real_evenChunkMiss, file = "./data/missingDatasets/pois_real_evenChunkMiss.rds")
saveRDS(pois_real_randChunkMiss, file = "./data/missingDatasets/pois_real_randChunkMiss.rds")
saveRDS(pois_real_minMaxMiss, file = "./data/missingDatasets/pois_real_minMaxMiss.rds")

