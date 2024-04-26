## this script trims those simulated population datasets that went extinct

# read in missing data
poiss_A <- readRDS("./data/missingDatasets/pois_sim_randMiss_A.rds")
poiss_B <- readRDS("./data/missingDatasets/pois_sim_randMiss_B.rds")

 
# get those datasets that include a 0, and trim them to end at the year when the 0 occurs
extinctDat_A <- lapply(poiss_A, function(x) {
  temp <- sapply(x$y, simplify = FALSE, function(y) {
    # is there a 0 in the string?
    if (sum(y==0, na.rm = TRUE) > 0) { 
      # is the 0 NOT at the end of the string?
      if (which(y == 0)[1] < length(y)) { 
        # if there is a zero that is not at the end of the string, then end the ts at that value
        y_new <- y[1:which(y==0)[1]]
        return(y_new)
      }
    } 
  })
  # assign names with accurate values for propMiss and autoCorr
  names(temp) <- sapply(temp, function(z) {
    ifelse(is.null(z),
           yes = paste0("no_zeros"),
           no = 
    paste0("propMissAct_",round(sum(is.na(z))/ length(z),2), 
           "_autoCorr_",round(acf(x = z, plot = FALSE, lag.max = 5, na.action = na.pass)$acf[2],2))
    )
  })
  
  if (sum(names(temp) %in% c("no_zeros")) > 0) {
    temp2 <- NULL
  } else {
    # remove any duplicates (arise if the trimmed extinct ts are identical)
    temp <- temp[!duplicated(temp)]
    # put into an output list
    temp2 <- list("y" = temp)
  }
  return(temp2)
  })

extinctDat_B <- lapply(poiss_B, function(x) {
  temp <- sapply(x$y,  simplify = FALSE, function(y) {
    # is there a 0 in the string?
    if (sum(y==0, na.rm = TRUE) > 0) { 
      # is the 0 NOT at the end of the string?
      if (which(y == 0)[1] < length(y)) { 
        # if there is a zero that is not at the end of the string, then end the ts at that value
        y_new <- y[1:which(y==0)[1]]
        return(y_new)
      }
    } 
  })
  # assign names with accurate values for propMiss and autoCorr
  names(temp) <- sapply(temp, function(z) {
    ifelse(is.null(z),
           yes = paste0("no_zeros"),
           no = 
             paste0("propMissAct_",round(sum(is.na(z))/ length(z),2), 
                    "_autoCorr_",round(acf(x = z, plot = FALSE, lag.max = 5, na.action = na.pass)$acf[2],2))
    )
  })
  
  if (sum(names(temp) %in% c("no_zeros")) > 0) {
    temp2 <- NULL
  } else {
    # remove any duplicates (arise if the trimmed extinct ts are identical)
    temp <- temp[!duplicated(temp)]
    # put into an output list
    temp2 <- list("y" = temp)
  }
  return(temp2)
})
  

# keep only the outputs with zeros
extinct_B_temp <- extinctDat_B %>% 
  keep(function(x)  !is.null(x)) 

extinct_A_temp <- extinctDat_A %>% 
  keep(function(x)  !is.null(x))

extinctDat <- append(extinct_A_temp, extinct_B_temp)   

# write to file
saveRDS(extinctDat, file = "./data/missingDatasets/pois_sim_randMiss_extinctions.rds")
