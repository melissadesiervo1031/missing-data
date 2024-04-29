# Load packages ## 
#make sure these are already in the folder on supercomputer where I need them ##

.libPaths("/pfs/tc1/home/astears/R/x86_64-pc-linux-gnu-library/4.2")

library(tidyverse)
#library(forecast) ## it hates this package...run with lowercase arima# 
library(Amelia)

# This script will run 3 ARIMA functions (drop missing, Kalman, Multiple imputations 
#over a nested list with increasing prop missing, over 1000+ simulations ###)

#CurSim = like a loop ##

CurSim <- commandArgs(trailingOnly = TRUE) #Look at command line arguments only after the R script
CurSim <- as.numeric(CurSim)
CurSim <- CurSim + 1 # since the Slurm array is 0 indexed

## read in the autocor_01 list ##

gauss_sim_MNAR <- readRDS("/project/modelscape/analyses/MissingTS/gauss_sim_minMaxMiss.rds")
#gauss_sim_MNAR <- readRDS("./data/missingDatasets/gauss_sim_minMaxMiss.rds")
# make file for output beforehand in supercomputer folder 
# will put them all together after all run, using the command line
OutFile <- paste("gauss_sim_minMax_modResults/", CurSim, "arimavals.csv", sep = "")

#########################################################################################
### MY ARIMA FUNCTIONS #####
##########################################################################################

### Function that will drop missing values and then fit model using ARIMA ###

fit_arima_dropmissing <- function(sim_list, sim_pars){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$X[,2], 
                                                           discharge = sim_pars$X[,3]))
  ## drop the missing values ###
  
  sim_missing_list_drop <- lapply(seq_along(simmissingdf), function(j) {
    drop_na(simmissingdf[[j]])
  })
  
  
  # fit arima models to list of datasets
  Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
    xreg1<-sim_missing_list_drop [[j]][["light"]]
    xreg2<-sim_missing_list_drop [[j]][["discharge"]]
    modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
    arimacoefsdrop<-modeldrop$coef
    names(arimacoefsdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
    arimasesdrop<-sqrt(diag(vcov(modeldrop)))
    names(arimasesdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
    list(arimacoefsdrop=arimacoefsdrop, arimasesdrop=arimasesdrop)
    
    return(list(arima_pars = arimacoefsdrop,
                arima_errors = arimasesdrop,
                sim_params = sim_pars))
  })
}

# fit complete case drop missing
fit_arima_dropmissing_CC <- function(sim_list, sim_pars){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$X[,2], 
                                                           discharge = sim_pars$X[,3]))
  # remove data in a "complete case" way
  # compile into sliced dataframe
  sim_missing_list_drop <- map(simmissingdf, function(x) {
    temp <- data.frame(
      yt = x[2:nrow(x),],
      ytm1 = x[1:(nrow(x)-1),]
    )
    # drop incomplete cases
    x[complete.cases(temp),]
  }
  )
  
  # fit arima models to list of datasets
  Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
    xreg1<-sim_missing_list_drop [[j]][["light"]]
    xreg2<-sim_missing_list_drop [[j]][["discharge"]]
    modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
    arimacoefsdrop<-modeldrop$coef
    names(arimacoefsdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
    arimasesdrop<-sqrt(diag(vcov(modeldrop)))
    names(arimasesdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
    list(arimacoefsdrop=arimacoefsdrop, arimasesdrop=arimasesdrop)
    
    return(list(arima_pars = arimacoefsdrop,
                arima_errors = arimasesdrop,
                sim_params = sim_pars))
  })
}
  


### Function that will have missing values as NA and then fit model using ARIMA w/ Kalman filter ###

fit_arima_Kalman <- function(sim_list, sim_pars){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$X[,2], 
                                                           discharge = sim_pars$X[,3]))
  
  ## fit ARIMA with the missing values as NAS . Applies KALMAN FILTER###
  
  
  ArimaoutputNAs <- lapply(seq_along(simmissingdf), function(j) {
    xreg1<-simmissingdf [[j]][["light"]]
    xreg2<-simmissingdf [[j]][["discharge"]]
    modelNAs <- arima(simmissingdf[[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
    arimacoefsNAs<-modelNAs$coef
    names(arimacoefsNAs) <- c("ar1", "intercept", "xreg1", "xreg2")
    arimasesNAs<-sqrt(diag(vcov(modelNAs)))
    names(arimasesNAs) <- c("ar1", "intercept", "xreg1", "xreg2")
    list(arimacoefsNAs=arimacoefsNAs, arimasesNAs=arimasesNAs)
    
    return(list(arima_pars = arimacoefsNAs,
                arima_errors = arimasesNAs,
                sim_params = sim_pars))
  })
}


###### 

### Function that will impute missing values w/ AMELIA and then fit model using ARIMA ###

fit_arima_MI <- function(sim_list, sim_pars, imputationsnum){
  
  days<-seq(1, 365)
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(days= days,
                                                           GPP = X, 
                                                           light = sim_pars$X[,2], 
                                                           discharge = sim_pars$X[,3]))
  
  amelia1sim <-lapply(X = simmissingdf  , FUN = function(X)   amelia(X, ts="days", m=imputationsnum, lags="GPP")) ## lags by 1 day ##
  
  
  ##nested list of dataframes that just has the imputations###
  amelias11sim<-map(amelia1sim , ~.[["imputations"]])
  
  ##forloop## gives us the model parameters and errors for all the ARIMA models on imputed datasets
  
  modelparamlistsim=list()
  modelerrorlistsim=list()
  
  for (i in seq_along(amelias11sim)) {
    a=list()
    aa=list()
    for (j in seq_along(amelias11sim[[i]])) {
      xreg1<-amelias11sim [[i]][[j]][["light"]]
      xreg2<-amelias11sim [[i]][[j]][["discharge"]]
      tempobj=arima(amelias11sim[[i]][[j]]$GPP, order = c(1,0,0), xreg = matrix(c(xreg1, xreg2), ncol = 2))
      arimacoefs<-tempobj$coef
      names(arimacoefs) <- c("ar1", "intercept", "xreg1", "xreg2")
      arimases<-sqrt(diag(vcov(tempobj)))
      names(arimases) <- c("ar1", "intercept", "xreg1", "xreg2")
      name <- paste('imp',seq_along((amelias11sim)[[i]])[[j]],sep='')
      a[[name]] <- arimacoefs
      aa[[name]]<-arimases
    }
    name1 <- names(amelias11sim)[[i]]
    modelparamlistsim[[name1]] <- a
    modelerrorlistsim[[name1]] <- aa
  }
  
  
  ### Averages the models together back to 1 model per missing data prop ##
  
  listcoefsessim<-mapply(function(X,Y) {
    list(mi.meld(data.frame(X), data.frame(Y), byrow=FALSE))
  }, X=modelparamlistsim, Y=modelerrorlistsim)
  
  return(list(paramlistsim<-map(listcoefsessim , ~.["q.mi"]),
              selistsim<-map(listcoefsessim , ~.["se.mi"])))
  
}





#####################################################
#### MODEL RUN ARIMA DROP ##############
#########################################################

arima_drop_MNAR<- fit_arima_dropmissing(gauss_sim_MNAR[[CurSim]]$y,gauss_sim_MNAR[[CurSim]]$sim_params)


########### formatting for figure #############

names(arima_drop_MNAR) <- names(gauss_sim_MNAR[[CurSim]][["y"]])

modeldropparamlist<-purrr::map(arima_drop_MNAR , ~.["arima_pars"])
modeldropSElist<-purrr::map(arima_drop_MNAR , ~.["arima_errors"])

modeldropparamlist2 <- lapply(modeldropparamlist, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElist2 <- lapply(modeldropSElist, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdf <- map_df(modeldropparamlist2, ~as.data.frame(.x), .id="missingprop_autocor")
modeldropSEdf <- map_df(modeldropSElist2, ~as.data.frame(.x), .id="missingprop_autocor")

modeldropdf<-modeldropparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  select(missingprop_autocor, ar1, intercept, light, discharge)  %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion Simple")

modeldropSEdf<-modeldropSEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2)%>%   select(missingprop_autocor, ar1, intercept, light, discharge) %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion Simple")

## long form ##

paramdroplong <- gather(modeldropdf, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElong <- gather(modeldropSEdf, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplong2<-merge(paramdroplong, paramdropSElong)

#####################################################
#### MODEL RUN ARIMA DROP --complete case ##############
#########################################################

arima_drop_CC_MNAR<- fit_arima_dropmissing_CC(gauss_sim_MNAR[[CurSim]]$y,gauss_sim_MNAR[[CurSim]]$sim_params)


########### formatting for figure #############

names(arima_drop_CC_MNAR) <- names(gauss_sim_MNAR[[CurSim]][["y"]])

modeldropCCparamlist<-purrr::map(arima_drop_CC_MNAR , ~.["arima_pars"])
modeldropCCSElist<-purrr::map(arima_drop_CC_MNAR , ~.["arima_errors"])

modeldropCCparamlist2 <- lapply(modeldropCCparamlist, function(x) as.data.frame(do.call(rbind, x)))
modeldropCCSElist2 <- lapply(modeldropCCSElist, function(x) as.data.frame(do.call(rbind, x)))


modeldropCCparamdf <- map_df(modeldropCCparamlist2, ~as.data.frame(.x), .id="missingprop_autocor")
modeldropCCSEdf <- map_df(modeldropCCSElist2, ~as.data.frame(.x), .id="missingprop_autocor")

modeldropCCdf<-modeldropCCparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  select(missingprop_autocor, ar1, intercept, light, discharge)  %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion CC")

modeldropCCSEdf<-modeldropCCSEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2)%>%   select(missingprop_autocor, ar1, intercept, light, discharge) %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion CC")

## long form ##

paramdropCClong <- gather(modeldropCCdf, param, value, ar1:discharge, factor_key=TRUE)

paramdropCCSElong <- gather(modeldropCCSEdf, param, SE, ar1:discharge, factor_key=TRUE)

paramdropCClong2<-merge(paramdropCClong, paramdropCCSElong)


#####################################################
#### MODEL RUN KALMAN FILTER ##############
#########################################################


arima_kalman_MNAR<- fit_arima_Kalman(gauss_sim_MNAR[[CurSim]]$y,gauss_sim_MNAR[[CurSim]]$sim_params)


## pull out and label what we need ###

names(arima_kalman_MNAR) <- names(gauss_sim_MNAR[[CurSim]][["y"]])

modelNAparamlist<-purrr::map(arima_kalman_MNAR , ~.["arima_pars"])
modelNASElist<-purrr::map(arima_kalman_MNAR , ~.["arima_errors"])

modelNAparamlist2 <- lapply(modelNAparamlist, function(x) as.data.frame(do.call(rbind, x)))
modelNASElist2 <- lapply(modelNASElist, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdf <- map_df(modelNAparamlist2, ~as.data.frame(.x), .id="missingprop_autocor")
modelNASEdf <- map_df(modelNASElist2, ~as.data.frame(.x), .id="missingprop_autocor")


modelNAdf<-modelNAparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(missingprop_autocor, ar1, intercept, light, discharge) %>% mutate(missingness="MNAR") %>% mutate(type="Kalman filter")

modelNASEdf<-modelNASEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(missingprop_autocor, ar1, intercept, light, discharge) %>% mutate(missingness="MNAR") %>% mutate(type="Kalman filter")


## long form ##

paramNAlong <- gather(modelNAdf, param, value, ar1:discharge, factor_key=TRUE)

paramNASElong <- gather(modelNASEdf, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlong2<-merge(paramNAlong, paramNASElong)


#######################################################################################################


#####################################################
#### MODEL RUN MULTIPLE IMPUTATIONS  ##############
#########################################################

arima_mi_MNAR<- fit_arima_MI(gauss_sim_MNAR[[CurSim]]$y,gauss_sim_MNAR[[CurSim]]$sim_params, imputationsnum=5)

##pulls out parameters and ses ##

paramlistsim<-arima_mi_MNAR[[1]]

selistsim<-arima_mi_MNAR[[2]]

avgparamdf <- map_df(paramlistsim, ~as.data.frame(.x), .id="missingprop_autocor")
avglSEdf <- map_df(selistsim, ~as.data.frame(.x), .id="missingprop_autocor")


avgparamdf2<-avgparamdf %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>%  select(missingprop_autocor,  ar1, intercept, light, discharge)  %>% mutate(missingness="MNAR") %>% mutate(type="Multiple imputations")

avglSEdf2<-avglSEdf  %>% dplyr::rename(ar1=se.mi.ar1, intercept=se.mi.intercept, light=se.mi.xreg1, discharge=se.mi.xreg2) %>%  select(missingprop_autocor,  ar1, intercept, light, discharge)   %>% mutate(missingness="MNAR") %>% mutate(type="Multiple imputations")


paramMIlong <- gather(avgparamdf2, param, value, ar1:discharge, factor_key=TRUE)

paramMISElong <- gather(avglSEdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramMIlong2<-merge(paramMIlong,paramMISElong)

#############################################################################################################

###################################################
#### COMBINE ALL MODEL RUNS AND SAVE #########
#################################################

###
paramarimaall<-rbind(paramdroplong2, paramNAlong2, paramMIlong2)


# 
Output <- matrix(data=NA, nrow=nrow(paramarimaall), ncol=ncol(paramarimaall))



# Save the results of the current script's simulation to the appropriate column of output
Output<- paramarimaall

# add in the "simulation number" for this iteration (which is stored in the name of the data list element)

Output2<-cbind(CurSim, Output)


# Write the output to the folder which will contain all output files as separate csv
#    files with a single line of data.
write.csv(Output2, file = OutFile, row.names = FALSE)


# Once the job finishes, you can use the following command from within the folder
#    containing all single line csv files to compile them into a single csv file:
#     awk '(NR == 1) || (FNR > 1)' *vals.csv > AllResults.csv
# The * is a wildcard character so the input to this will match any file within
#    your current folder that ends with vals.csv regardless of the rest of the filename.
#    These will then all be combined into a single file (AllResults.csv). The order
#    will be based on how linux orders the file names within the directory, so it 
#    might not match the original order of your parameter input file, but all the
#    entries will be there and it can be sorted later. Alternatively, you can name
#    your output files in a way in which the proper order will be enforced (e.g.,
#    if you will have a total of 100 jobs, you can name them all with 3 digits like
#    001_vals.csv, 002_vals.csv, etc.)
# Once you have combined all the single line csv files into your master results file,
#    you can remove them using the wildcard character again (e.g., rm *vals.csv)



