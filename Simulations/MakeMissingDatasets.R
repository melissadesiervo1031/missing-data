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

# Gaussian Simulated Data  - curtailed for forecasting -------------------------------------------------
# read in data
gauss_sim <- readRDS("./data/gauss_ar1_0miss_datasets.rds")

## make missing data types for increasing levels of autocorrelation

# read in data


## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)
inputPropMiss <- c(.18, .19, .2, .21, .22, .38, .39, .4, .41, .42, .58, .59, .6, .61, .62)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness
  # remove 20% of observations from the end as a part of generating missing values
  tempOutList <- lapply(X = gauss_sim, FUN = function(x) 
    list("y" = makeMissing(timeSeries = x$y[1:292], typeMissing = "random", 
                           autoCorr = inputAutocor[i],
                           propMiss = inputPropMiss),
         "sim_params" <- x$sim_params)
  )
  # name the elements of the list with the amount of missingness 
  for (j in 1:length(tempOutList)) {
    names(tempOutList[[j]]) <- c("y", "sim_params")
    tempOutList[[j]]$y <- c(list("y_noMiss" = gauss_sim[[j]]$y), tempOutList[[j]]$y)
  }
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("gauss_sim_randMiss_trim_autoCorr_0") , 
           value = tempOutList)
  } else {
    assign(x = paste0("gauss_sim_randMiss_trim_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutList)
  }
}

# for (i in 1:length(gauss_sim_randMiss_trim_autoCorr_10)) {
#   MissFromName <- as.numeric(str_extract(names(gauss_sim_randMiss_trim_autoCorr_10[[i]]$y), pattern = "0.[0-9]++"))
#   MissFromData <- sapply(gauss_sim_randMiss_trim_autoCorr_10[[i]]$y,
#                          FUN = function(x) {
#                            sum(is.na(x))/length(x)
#                          }
#   )
#   temp <- data.frame("missFromName" = MissFromName,
#                      "missFromData" = MissFromData,
#                      "diff" = MissFromName - MissFromData,
#                      "i" = i)
#   if (i == 1) {
#     outDat <- temp
#   } else {
#     outDat <- rbind(outDat, temp)
#   }
# }
# plot(density(outDat[!is.na(outDat$diff), "diff"]))

# missing in min and max of data
gauss_sim_minMaxMiss_trim <- lapply(X = gauss_sim, FUN = function(x) 
  list("y" = makeMissing(timeSeries = x$y[1:292], typeMissing = "minMax", type =  "Gaussian", 
                         propMiss = inputPropMiss ), 
       
       "sim_params" <- x$sim_params)
)

for (i in 1:length(gauss_sim_minMaxMiss_trim)) {
  names(gauss_sim_minMaxMiss_trim[[i]]) <- c("y", "sim_params")
  gauss_sim_minMaxMiss_trim[[i]]$y <- c(list("y_noMiss" = gauss_sim[[i]]$y), gauss_sim_minMaxMiss_trim[[i]]$y)
  
}

# # check missingness values included in names
# for (i in 1:length(gauss_sim_minMaxMiss_trim)) {
#   MissFromName <- as.numeric(str_extract(names(gauss_sim_minMaxMiss_trim[[i]]$y), pattern = "[0-9]\\.[0-9]+$"))
#   MissFromData <- sapply(gauss_sim_minMaxMiss_trim[[i]]$y,
#                           FUN = function(x) {
#                             sum(is.na(x))/length(x)
#                           }
#                          )
#   temp <- data.frame("missFromName" = MissFromName,
#              "missFromData" = MissFromData,
#              "diff" = MissFromName - MissFromData,
#              "i" = i)
#   if (i == 1) {
#     outDat <- temp
#   } else {
#     outDat <- rbind(outDat, temp)
#   }
# }
# plot(density(outDat$diff))

# Gaussian Real Data (Au Sable River)  - curtailed for forecasting--------------------------------
# read in data
gauss_auSable <- read.csv("./data/au_sable_river_prepped.csv")

# cut off the last 10 values before going to missingness procedure so we don't need to recalculate
gauss_auSable <- gauss_auSable[1:731,]

## make missing data types
## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.25, .50, .75)
reps <- 50
tolerance <- 0.05
listNames <- vector()
for (i in 1:(reps*length(inputAutocor))) {
  # calculate missing vectors with increasing levels of missingness
  
  # loop for making sure propMiss and inputAutocor are within a short distance of desired values (for both + or - 0.05 from desired value)
  check=F
  while(!check){
    propMiss=c(0.2,0.4,0.6)
    tempOutDf <- as.data.frame(makeMissing(timeSeries = gauss_auSable$GPP, propMiss = propMiss,
                                           typeMissing = "random", 
                                           autoCorr = inputAutocor[((i-1)%%length(inputAutocor))+1]))
    
    
    miss_actual <- as.numeric(sub(".*propMissAct_([0-9\\.]+).*", "\\1", names(tempOutDf)))
    print(miss_actual)
    autocor_actual <- as.numeric(sub(".*autoCorr_([0-9\\.]+).*", "\\1", names(tempOutDf)))
    print(autocor_actual)
    
    # check
    if(abs(miss_actual[1]-propMiss[1])<=tolerance&abs(miss_actual[2]-propMiss[2])<=tolerance&abs(miss_actual[3]-propMiss[3])<=tolerance&
       abs(autocor_actual[1]-inputAutocor[((i-1)%%length(inputAutocor))+1])<=tolerance&abs(autocor_actual[2]-inputAutocor[((i-1)%%length(inputAutocor))+1])<=tolerance&abs(autocor_actual[2]-inputAutocor[((i-1)%%length(inputAutocor))+1])<=tolerance){
      check=T
    } else {
      print(paste0("retrying on i=",i))
    }
  }
  
  # name the elements of the Df with the amount of missingness 
  names(tempOutDf) <- paste0("GPP_",names(tempOutDf))
  tempOutDf <- cbind(gauss_auSable, tempOutDf)
  
  # rename the output list to reflect the input autocorrelation
  # if (((i-1)%%length(inputAutocor))+1 == 1) {
  #   assign(x = paste0("gauss_auSable_randMiss_trim_autoCorr_0_i",rep(1:50, each = 4)[i]) , 
  #          value = tempOutDf)
  # } else 
  {
    assign(x = paste0("gauss_auSable_randMiss_trim_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[((i-1)%%length(inputAutocor))+1], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0"),"_i",rep(1:50, each = 4)[i]) , 
           value = tempOutDf)
    listNames[i] <- paste0("gauss_auSable_randMiss_trim_autoCorr_", 
                           str_pad(str_extract_all(string = inputAutocor[((i-1)%%length(inputAutocor))+1], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0"),"_i",rep(1:50, each = 4)[i])
  }
}


# put all of the data into one list
gauss_auSable_randMiss_trim_list <- vector(mode = "list", length = 150) 
names(gauss_auSable_randMiss_trim_list) <- listNames
for (i in 1:length(listNames)) {
  # get the correct autocorrelation/missing data data.frame
  tempDF <- get(x = listNames[i])
  # change into a list element
  gauss_auSable_randMiss_trim_list[[i]]$y <- tempDF %>% 
    rename_with(~ str_replace(string = names(tempDF), pattern = "GPP_", replacement = "")) %>% 
    rename("y" = "GPP" ) %>% 
    dplyr::select(-date, -ER, -light, -Q) %>% 
    as.list()
  gauss_auSable_randMiss_trim_list[[i]]$sim_params$date <- tempDF %>% 
    dplyr::select(date)
  gauss_auSable_randMiss_trim_list[[i]]$sim_params$light <- tempDF %>% 
    dplyr::select(light)
  gauss_auSable_randMiss_trim_list[[i]]$sim_params$Q <- tempDF %>% 
    dplyr::select(Q)
}

gauss_auSable_randMiss_trim <- gauss_auSable_randMiss_trim_list


## missing in min and max of data
## for the MNAR data, we need to do slightly different propMissing values, since the same proportion will cutoff the same exact values 

inputPropMiss <- c(.18, .19, .2, .21, .22, .38, .39, .4, .41, .42, .58, .59, .6, .61, .62)
gauss_auSable_minMaxMiss_trim_TEMP <- as.data.frame(makeMissing(timeSeries = gauss_auSable$GPP, 
                                                                typeMissing = "minMax", 
                                                                propMiss = inputPropMiss, type =  "Gaussian"))


#names(gauss_auSable_minMaxMiss_trim_TEMP) <- paste0("GPP_",names(gauss_auSable_minMaxMiss_trim_TEMP))

gauss_auSable_minMaxMiss_trim <- cbind(gauss_auSable, gauss_auSable_minMaxMiss_trim_TEMP)

# change "broods" to "y" in the data.frame
names(gauss_auSable_minMaxMiss_trim)[2] <- "y" 

# transform into a list (for consistency w/ simulated data)
gauss_auSable_minMaxMiss_trim_list <- vector(mode = "list", length = 1) 
gauss_auSable_minMaxMiss_trim_list[[1]]$y <- as.list(gauss_auSable_minMaxMiss_trim[c(2,9:ncol(gauss_auSable_minMaxMiss_trim))])
gauss_auSable_minMaxMiss_trim_list[[1]]$sim_params$date <- gauss_auSable_minMaxMiss_trim$date 
gauss_auSable_minMaxMiss_trim_list[[1]]$sim_params$light <- gauss_auSable_minMaxMiss_trim$light
gauss_auSable_minMaxMiss_trim_list[[1]]$sim_params$light.rel <- gauss_auSable_minMaxMiss_trim$light.rel
gauss_auSable_minMaxMiss_trim_list[[1]]$sim_params$Q <- gauss_auSable_minMaxMiss_trim$Q

gauss_auSable_minMaxMiss_trim <- gauss_auSable_minMaxMiss_trim_list

# Poisson Simulated Data - curtailed for forecasting --------------------------------------------------
pois_sim <- readRDS("./data/ricker_0miss_datasets.rds")

# make missing data types

## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)

# remove 20% of observations
# read in data (pois_sim)
# trim to remove the last %20 of observations from each time series

## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.0, .10, .20, .30, .40, .50, .60, .70, .80, .90)

inputPropMiss <- c(.18, .19, .2, .21, .22, .38, .39, .4, .41, .42, .58, .59, .6, .61, .62)

for (i in 1:length(inputAutocor)) {
  # calculate missing vectors with increasing levels of missingness

  tempOutList <- lapply(X = pois_sim, FUN = function(x)  {
    last_int <- round(length(x$y)*.8) # identifying the last point of the time series w/ 20% removed for forecasting
    list("y" = makeMissing(timeSeries = x$y[1:last_int], typeMissing = "random", autoCorr = inputAutocor[i], 
                           propMiss = inputPropMiss),
         "sim_params" <- x$sim_params)
  }
  )
  # name the elements of the list with the amount of missingness 
  for (j in 1:length(tempOutList)) {
    names(tempOutList[[j]]) <- c("y", "sim_params")
    tempOutList[[j]]$y <- c(list("y_noMiss" = pois_sim[[j]]$y), tempOutList[[j]]$y)
  }
  # rename the output list to reflect the input autocorrelation
  if (i == 1) {
    assign(x = paste0("pois_sim_randMiss_autoCorr_0_20PercTrimmed_") , 
           value = tempOutList)
  } else {
    assign(x = paste0("pois_sim_randMiss_autoCorr_20PercTrimmed_", 
                      str_pad(str_extract_all(string = inputAutocor[i], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0")) , 
           value = tempOutList)
  }
}

# missing in min and max of data
pois_sim_minMaxMiss_20PercTrimmed <- lapply(X = pois_sim, FUN = function(x) {
  last_int <- round(length(x$y)*.8)
  list("y" = makeMissing(timeSeries = x$y[1:last_int], typeMissing = "minMax", type = "Poisson", 
                         propMiss = inputPropMiss), 
       "sim_params" <- x$sim_params)
}
)

for (i in 1:length(pois_sim_minMaxMiss_20PercTrimmed)) {
  names(pois_sim_minMaxMiss_20PercTrimmed[[i]]) <- c("y", "sim_params")
  pois_sim_minMaxMiss_20PercTrimmed[[i]]$y <- c(list("y_noMiss" = pois_sim[[i]]$y), pois_sim_minMaxMiss_20PercTrimmed[[i]]$y)
}

# confirm accuracy of missingness values 
# for (i in 1:length(pois_sim_minMaxMiss_20PercTrimmed)) {
#   MissFromName <- as.numeric(str_extract(names(pois_sim_minMaxMiss_20PercTrimmed[[i]]$y), pattern = "0.[0-9]++"))
#   MissFromData <- sapply(pois_sim_minMaxMiss_20PercTrimmed[[i]]$y,
#                          FUN = function(x) {
#                            sum(is.na(x))/length(x)
#                          }
#   )
#   temp <- data.frame("missFromName" = MissFromName,
#                      "missFromData" = MissFromData,
#                      "diff" = MissFromName - MissFromData,
#                      "i" = i)
#   if (i == 1) {
#     outDat <- temp
#   } else {
#     outDat <- rbind(outDat, temp)
#   }
# }
# plot(density(outDat[!is.na(outDat$diff), "diff"]))

# Poisson Real Data (Protist data) - full dataset -------------------------------
# for these data, because there are 10 different replicates and we want the missingness to be across all replicates, if we want 10% missingness across the replicates, we introduce 1% missingness to each replicate
# read in data
bursaria <- read.csv('data/Bursaria.csv')

#  Data formatting 

bursaria <- bursaria %>% 
  janitor::clean_names() %>%
  mutate(
    date = lubridate::mdy(date),
    patch = factor(patch)
  )

# grouping by patch, then adding a number of days column
bursaria_patchlist <- purrr::map(
  unique(bursaria$patch),
  ~ {
    filter(bursaria, patch == .x) %>%
      arrange(date) %>%
      mutate(
        days = as.numeric(date - date[1], units = "days")
      )
  }
)

# now offset each so that we have a response that is time t and explanatory variable that is
# time t - 1
bursaria_patchlist_diff <- purrr::map(
  bursaria_patchlist,
  ~ {
    df2 <- tibble(
      count_tm1 = .x$number[-nrow(.x)]
    )
    cbind(.x[-1, ], df2)
  }
)

# recombine into one dataframe
bursaria_diff <- Reduce(rbind, bursaria_patchlist_diff)
#bursaria2 <- Reduce(rbind, bursaria_patchlist)

# cut off the last 10 values before going to missingness procedure so we don't need to recalculate
#pois_real=pois_real[1:49,]

## make missing data types
## make missing data types for increasing levels of autocorrelation
# possible autocorrelation vector
inputAutocor <- c(.25, .50, .75)
reps <- 50
tolerance <- .1#0.05
patchNames <- unique(bursaria_diff$patch)
for (i in 1:(reps*length(inputAutocor))) {
  # calculate missing vectors with increasing levels of missingness
  
  # loop for making sure propMiss and inputAutocor are within a short distance of desired values (for both + or - 0.05 from desired value)
  check=F
  while(!check){
    propMiss=c(0.2,0.4,0.6)
    # for each replicate, remove propMiss percent of the observations
    
    tempOutList <- lapply(X = 1:10, FUN = function(x) {
      temp <- cbind(bursaria_diff[bursaria_diff$patch == patchNames[x],c("patch", "date", "number", "days")],
                    
        as.data.frame(makeMissing(timeSeries = bursaria_diff[bursaria_diff$patch == patchNames[x],"number"], 
                                propMiss = propMiss,
                                typeMissing = "random", 
                                autoCorr = inputAutocor[((i-1)%%length(inputAutocor))+1])))
      # save autocorrelation values
      temp$miss1_autocorr <- as.numeric(sub(".*autoCorr_([0-9\\.]+).*", "\\1" ,names(temp)[5]))
      temp$miss2_autocorr <- as.numeric(sub(".*autoCorr_([0-9\\.]+).*", "\\1" ,names(temp)[6]))
      temp$miss3_autocorr <- as.numeric(sub(".*autoCorr_([0-9\\.]+).*", "\\1" ,names(temp)[7]))
      names(temp)[5:7] <- c("propMiss_0.20", "propMiss_0.40", "propMiss_0.60")
      return(temp)
    }
    ) %>% 
      purrr::list_rbind()
    # calculate the amount of missingness and autocorrelation across all patches
    names(tempOutList)[5] <- paste0("propMissAct_",
      round(sum(is.na(tempOutList$propMiss_0.20))/length(tempOutList$propMiss_0.20),2), # calculate total missingness across patches
      "_autoCorr_",
      # calculate average autocorrelation across patches 
      round(mean(unique(tempOutList$miss1_autocorr), na.rm = TRUE),2)
    )
    names(tempOutList)[6] <- paste0("propMissAct_",
                                    round(sum(is.na(tempOutList$propMiss_0.40))/length(tempOutList$propMiss_0.40),2), # calculate total missingness across patches
                                    "_autoCorr_",
                                    # calculate average autocorrelation across patches 
                                    round(mean(unique(tempOutList$miss2_autocorr), na.rm = TRUE),2)
    )
    names(tempOutList)[7] <- paste0("propMissAct_",
                                    round(sum(is.na(tempOutList$propMiss_0.60))/length(tempOutList$propMiss_0.60),2), # calculate total missingness across patches
                                    "_autoCorr_",
                                    # calculate average autocorrelation across patches 
                                   round(mean(unique(tempOutList$miss3_autocorr), na.rm = TRUE),2)
    )
    
    
    miss_actual <- as.numeric(sub(".*propMissAct_([0-9\\.]+).*", "\\1", names(tempOutList)[5:7]))
    print(miss_actual)
    autocor_actual <- as.numeric(sub(".*autoCorr_([0-9\\.]+).*", "\\1", names(tempOutList)[5:7]))
    print(autocor_actual)
    
    # check
    if(abs(miss_actual[1]-propMiss[1])<=tolerance&abs(miss_actual[2]-propMiss[2])<=tolerance&abs(miss_actual[3]-propMiss[3])<=tolerance&
       abs(autocor_actual[1]-inputAutocor[((i-1)%%length(inputAutocor))+1])<=tolerance&abs(autocor_actual[2]-inputAutocor[((i-1)%%length(inputAutocor))+1])<=tolerance&abs(autocor_actual[3]-inputAutocor[((i-1)%%length(inputAutocor))+1])<=tolerance){
      check=T
    } else {
      print(paste0("retrying on i=",i))
    }
  }
  
  # name the elements of the Df with the amount of missingness 
  # names(tempOutList) <- paste0("Broods_",names(tempOutList))
  # tempOutDf <- cbind(pois_real, tempOutDf)
  # 
  # rename the output list to reflect the input autocorrelation
  # if (((i-1)%%length(inputAutocor))+1 == 1) {
  #   assign(x = paste0("pois_real_randMiss_autoCorr_0_i",rep(1:50, each = 4)[i]) , 
  #          value = tempOutDf)
  # } else 
    
    assign(x = paste0("pois_real_randMiss_autoCorr_", 
                      str_pad(str_extract_all(string = inputAutocor[((i-1)%%length(inputAutocor))+1], pattern = "\\d+" , simplify = TRUE)[,2], width = 2, side = "right", pad = "0"),"_i",rep(1:50, each = 4)[i]) , 
           value = tempOutList)
}

# put all of the data into one list
autcorVector <- c(25, 50, 75)
pois_real_randMiss_list <- vector(mode = "list", length = length(inputAutocor)*reps)
names(pois_real_randMiss_list) <- paste0("pois_real_randMiss_autoCor_", rep(autcorVector,reps),"_i",rep(1:reps,each=3))
for (i in 1:(reps*length(inputAutocor))) {
  # get the correct autocorrelation/missing data data.frame
  tempDF <- get(x = paste0("pois_real_randMiss_autoCorr_", autcorVector[((i-1)%%length(inputAutocor))+1],"_i",rep(1:50, each = 4)[i]))
  # change into a list element
  pois_real_randMiss_list[[i]]$y <- tempDF %>%
    #rename_with(~ str_replace(string = names(tempDF), pattern = "Broods_", replacement = "")) %>%
    #rename("y" = "count" ) %>%
    #dplyr::select(-Year) %>%
    dplyr::select(patch, number, days, 5, 6, 7) %>% 
    as.list()
  pois_real_randMiss_list[[i]]$sim_params <- NA
}
pois_real_randMiss_full <- pois_real_randMiss_list
# 
# # confirm accuracy of missingness values
# for (i in 1:length(pois_real_randMiss_full)) {
#   MissFromName <- as.numeric(str_extract(names(pois_real_randMiss_full[[i]]$y), pattern = "0.[0-9]++"))
#   MissFromData <- sapply(pois_real_randMiss_full[[i]]$y,
#                          FUN = function(x) {
#                            sum(is.na(x))/length(x)
#                          }
#   )
#   temp <- data.frame("missFromName" = MissFromName,
#                      "missFromData" = MissFromData,
#                      "diff" = MissFromName - MissFromData,
#                      "i" = i)
#   if (i == 1) {
#     outDat <- temp
#   } else {
#     outDat <- rbind(outDat, temp)
#   }
# }
# plot(density(outDat[!is.na(outDat$diff), "diff"]))

# MNAR missingness 
## missing in min and max of data
  # calculate missing vectors with increasing levels of missingness
  
    propMiss=c(.18, .19, .2, .21, .22, .38, .39, .4, .41, .42, .58, .59, .6, .61, .62)
    # for each replicate, remove propMiss percent of the observations
        tempOutList <- cbind(bursaria_diff[,c("patch", "date", "number", "days")],
                    as.data.frame(makeMissing(timeSeries = bursaria_diff[,"number"], 
                                              propMiss = propMiss,
                                              typeMissing = "minMax", type = "Poisson")))
      names(tempOutList)[5:19] <- c("propMiss_0.18", "propMiss_0.19", "propMiss_0.2", "propMiss_0.21", "propMiss_0.22", "propMiss_0.38", 
                            "propMiss_0.39", "propMiss_0.4", "propMiss_0.41", "propMiss_0.42", "propMiss_0.58", 
                            "propMiss_0.59", "propMiss_0.6", "propMiss_0.61", "propMiss_0.62")

    # calculate the amount of missingness and autocorrelation across all patches
      
    names(tempOutList)[5:19] <- paste0("propMissAct_",apply(tempOutList[,5:19], MARGIN = 2, FUN = function(x) {
      round(sum(is.na(x))/length(x),2)
      }), 
      "_i", seq(1:15))
   
  # name the elements of the Df with the amount of missingness 
  # names(tempOutList) <- paste0("Broods_",names(tempOutList))
  # tempOutDf <- cbind(pois_real, tempOutDf)
  # 
  # rename the output list to reflect the input autocorrelation
  # if (((i-1)%%length(inputAutocor))+1 == 1) {
  #   assign(x = paste0("pois_real_randMiss_autoCorr_0_i",rep(1:50, each = 4)[i]) , 
  #          value = tempOutDf)
  # } else 

# put all of the data into an organized list
pois_real_minMaxMiss_full <- list(
  "y" = tempOutList %>% 
    dplyr::select(patch, number, days, 5:19) %>% 
    as.list(), 
  "sim_params" = NA
)


# Store missing data  -----------------------------------------------------
# all datasets will be stored in "data/missingDatasets/"

# if it doesn't exist, make a folder to hold the datasets
if (dir.exists("./data/missingDatasets") == FALSE) {
  dir.create("./data/missingDatasets")
}

# prepare datastes for Beartooth runs -------------------------------------
## do the same for the Gaussian sim data (truncated dataset)
## bind missing datasets (Gauss sim) into sets of 5000 nested lists, rather than 1000
# for each "gauss_sim" list, name the sublist for each simulation
names(gauss_sim_randMiss_trim_autoCorr_0) <-  paste0("gauss_sim",1:1000, "_randMiss_autoCorr_0")
names(gauss_sim_randMiss_trim_autoCorr_10) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_10")
names(gauss_sim_randMiss_trim_autoCorr_20) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_20")
names(gauss_sim_randMiss_trim_autoCorr_30) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_30")
names(gauss_sim_randMiss_trim_autoCorr_40) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_40")
names(gauss_sim_randMiss_trim_autoCorr_50) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_50")
names(gauss_sim_randMiss_trim_autoCorr_60) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_60")
names(gauss_sim_randMiss_trim_autoCorr_70) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_70")
names(gauss_sim_randMiss_trim_autoCorr_80) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_80")
names(gauss_sim_randMiss_trim_autoCorr_90) <- paste0("gauss_sim",1:1000, "_randMiss_autoCorr_90")

gauss_sim_randMiss_trim_A <- c(gauss_sim_randMiss_trim_autoCorr_0, 
                               gauss_sim_randMiss_trim_autoCorr_10,
                               gauss_sim_randMiss_trim_autoCorr_20,
                               gauss_sim_randMiss_trim_autoCorr_30,
                               gauss_sim_randMiss_trim_autoCorr_40
)
gauss_sim_randMiss_trim_B <- c(gauss_sim_randMiss_trim_autoCorr_50,
                               gauss_sim_randMiss_trim_autoCorr_60,
                               gauss_sim_randMiss_trim_autoCorr_70,
                               gauss_sim_randMiss_trim_autoCorr_80,
                               gauss_sim_randMiss_trim_autoCorr_90)


gauss_sim_params <- data.frame(
  "SimNumber" = vector(mode = "integer", length = 1000),
  "phi" = vector(mode = "integer", length = 1000),
  "beta1" = vector(mode = "double", length = 1000),
  "beta2" = vector(mode = "double", length = 1000),
  "beta3" = vector(mode = "double", length = 1000)#,
  #"X" = vector(mode = "list", length = 1000)
)

for (i in 1:length(gauss_sim)) {
  gauss_sim_params[i, "SimNumber"] <- i
  gauss_sim_params[i,"phi"] <- gauss_sim[[i]]$sim_params$phi
  gauss_sim_params[i,"beta1"] <- gauss_sim[[i]]$sim_params$beta[1]
  gauss_sim_params[i,"beta2"] <- gauss_sim[[i]]$sim_params$beta[2]
  gauss_sim_params[i,"beta3"] <- gauss_sim[[i]]$sim_params$beta[3]
  #gauss_sim_params[i,"X"] <- list(gauss_sim[[i]]$sim_params$X)
}
# # Remove the y_no miss from all nested lists (Amelia doesn't like NO missing values) 
# # for first chunk of data
# #test <- gauss_sim_randMiss_trim_A
# test <- vector(mode = "list", length = 5000)
# for (i in 1:length(gauss_sim_randMiss_trim_A)) {
#   # remove the "y_noMiss" sub_list
#   test[[i]]$y <- gauss_sim_randMiss_trim_A[[i]]$y[2:16]
#   test[[i]]$sim_params$X <- gauss_sim_randMiss_trim_A[[i]]$sim$X
# }
# # rename w/ accurate names
# names(test) <- names(gauss_sim_randMiss_trim_A)
# gauss_sim_randMiss_trim_A <- test
# 
# # for second chunk of data
# test <- gauss_sim_randMiss_trim_B
# 
# for (i in 1:length(gauss_sim_randMiss_trim_B)) {
#   # remove the "y_noMiss" sub_list
#   test[[i]]$y <- gauss_sim_randMiss_trim_B[[i]]$y[2:16]
#   test[[i]]$sim_params$X <- gauss_sim_randMiss_trim_B[[i]]$sim$X
# }
# # rename w/ accurate names
# names(test) <- names(gauss_sim_randMiss_trim_B)
# gauss_sim_randMiss_trim_B <- test
# 

## do the same for the Ricker data (full dataset)
## bind missing datasets (pois sim) into sets of 5000 nested lists, rather than 1000
# for each "pois_sim" list, name the sublist for each simulation
names(pois_sim_randMiss_autoCorr_0_20PercTrimmed_) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_0")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_10) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_10")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_20) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_20")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_30) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_30")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_40) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_40")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_50) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_50")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_60) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_60")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_70) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_70")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_80) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_80")
names(pois_sim_randMiss_autoCorr_20PercTrimmed_90) <- paste0("pois_sim",1:1000, "_randMiss_autoCorr_90")

pois_sim_randMiss_trim_A <- c(pois_sim_randMiss_autoCorr_0_20PercTrimmed_, 
                              pois_sim_randMiss_autoCorr_20PercTrimmed_10,
                              pois_sim_randMiss_autoCorr_20PercTrimmed_20,
                              pois_sim_randMiss_autoCorr_20PercTrimmed_30,
                              pois_sim_randMiss_autoCorr_20PercTrimmed_40
)
pois_sim_randMiss_trim_B <- c(pois_sim_randMiss_autoCorr_20PercTrimmed_50,
                              pois_sim_randMiss_autoCorr_20PercTrimmed_60,
                              pois_sim_randMiss_autoCorr_20PercTrimmed_70,
                              pois_sim_randMiss_autoCorr_20PercTrimmed_80,
                              pois_sim_randMiss_autoCorr_20PercTrimmed_90)

# # Remove the y_no miss from all nested lists (Amelia doesn't like NO missing values) 
# # for first chunk of data
# #test <- pois_sim_randMiss_A
# test <- vector(mode = "list", length = 5000)
# for (i in 1:length(pois_sim_randMiss_trim_A)) {
#   # remove the "y_noMiss" sub_list
#   test[[i]]$y <- pois_sim_randMiss_trim_A[[i]]$y[2:16]
#   test[[i]]$sim_params$X <- pois_sim_randMiss_trim_A[[i]]$sim$X
# }
# # rename w/ accurate names
# names(test) <- names(pois_sim_randMiss_trim_A)
# pois_sim_randMiss_trim_A <- test
# 
# # for second chunk of data
# test <- pois_sim_randMiss_trim_B
# 
# for (i in 1:length(pois_sim_randMiss_trim_B)) {
#   # remove the "y_noMiss" sub_list
#   test[[i]]$y <- pois_sim_randMiss_trim_B[[i]]$y[2:16]
#   test[[i]]$sim_params$X <- pois_sim_randMiss_trim_B[[i]]$sim$X
# }
# # rename w/ accurate names
# names(test) <- names(pois_sim_randMiss_trim_B)
# pois_sim_randMiss_trim_B <- test


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
saveRDS(pois_sim_randMiss_trim_A, file = "./data/missingDatasets/pois_sim_randMiss_A.rds")
saveRDS(pois_sim_randMiss_trim_B, file = "./data/missingDatasets/pois_sim_randMiss_B.rds")

## save a data.frame that has the parameters used to run each simulation (1000 
# rows, each corresponding to a simulation run)
saveRDS(pois_sim_params, file = "./data/missingDatasets/pois_sim_params.rds")

## save poisson simulated MNAR data
saveRDS(pois_sim_minMaxMiss_20PercTrimmed, file = "./data/missingDatasets/pois_sim_MinMaxMiss.rds")

## save poisson real MNAR data (not trimmed, since we're going to do a LOO appraoch)
saveRDS(pois_real_minMaxMiss_full, file = "./data/missingDatasets/pois_real_MinMaxMiss.rds")
## save poisson real MAR data
saveRDS(pois_real_randMiss_full, file = "./data/missingDatasets/pois_real_randMiss.rds")


## store simulated Gaussian data (are stored in a list, each element of the list 
# is a simulation run. Within each simulation run, the $y element contains 15 
# elements that have the response variable ranging from the lowest amount of 
# missing data to the highest proportion of missing data. 
saveRDS(gauss_sim_randMiss_trim_A, file = "./data/missingDatasets/gauss_sim_randMiss_A.rds")
saveRDS(gauss_sim_randMiss_trim_B, file = "./data/missingDatasets/gauss_sim_randMiss_B.rds")

## save gaussian simulated MNAR data
saveRDS(gauss_sim_minMaxMiss_trim, file = "./data/missingDatasets/gauss_sim_minMaxMiss.rds")

## save a data.frame that has the parameters used to run each simulation (1000 
# rows, each corresponding to a simulation run)
saveRDS(gauss_sim_params, file = "./data/missingDatasets/gauss_sim_params.rds")

## save gaussian real MNAR data
saveRDS(gauss_auSable_minMaxMiss_trim, file = "./data/missingDatasets/gauss_real_auSable_MinMaxMiss.rds")
## save gaussian real MAR data
saveRDS(gauss_auSable_randMiss_trim, file = "./data/missingDatasets/gauss_real_auSable_randMiss.rds")

