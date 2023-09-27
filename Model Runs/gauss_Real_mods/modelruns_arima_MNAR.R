# Load packages ## 
library(tidyverse)
#library(forecast) ## it hates this package...run with lowercase arima# 
library(Amelia)

## read in the data ##

gauss_real_MinMaxMiss <- readRDS("./data/missingDatasets/gauss_real_MinMaxMiss.rds")

# make file for output beforehand in supercomputer folder 
# will put them all together after all run, using the command line
OutFile <- paste("./data/model_results/gauss_real_MNAR_arima.csv")

CurSim <- 1

# Drop missing + arima function ---------------------------------------------------

fit_arima_dropmissing <- function(sim_list, sim_pars){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$light, 
                                                           discharge = sim_pars$Q)) # Q is discharge
  ## drop the missing values ###
  
  sim_missing_list_drop <- lapply(seq_along(simmissingdf), function(j) {
    drop_na(simmissingdf[[j]])
  })
  
  # fit arima models to list of datasets
  Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
    xreg1<-sim_missing_list_drop [[j]][["light"]]
    xreg2<-sim_missing_list_drop [[j]][["discharge"]]
    modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
    arimacoefsdrop <-c(modeldrop$coef, modeldrop$sigma2)
    names(arimacoefsdrop) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
    arimasesdrop<-sqrt(diag(vcov(modeldrop)))
    names(arimasesdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
    list(arimacoefsdrop=arimacoefsdrop, arimasesdrop=arimasesdrop)
    
    return(list(arima_pars = arimacoefsdrop,
                arima_errors = arimasesdrop,
                sim_params = sim_pars))
  })
}

# arima + Kalman filter function ------------------------------------------

fit_arima_Kalman <- function(sim_list, sim_pars){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$light, 
                                                           discharge = sim_pars$Q))
  
  ## fit ARIMA with the missing values as NAS . Applies KALMAN FILTER###
  
  
  ArimaoutputNAs <- lapply(seq_along(simmissingdf), function(j) {
    xreg1<-simmissingdf [[j]][["light"]]
    xreg2<-simmissingdf [[j]][["discharge"]]
    modelNAs <- arima(simmissingdf[[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
    arimacoefsNAs <- c(modelNAs$coef, modelNAs$sigma2)
    names(arimacoefsNAs) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
    arimasesNAs<-sqrt(diag(vcov(modelNAs)))
    names(arimasesNAs) <- c("ar1", "intercept", "xreg1", "xreg2")
    list(arimacoefsNAs=arimacoefsNAs, arimasesNAs=arimasesNAs)
    
    return(list(arima_pars = arimacoefsNAs,
                arima_errors = arimasesNAs,
                sim_params = sim_pars))
  })
}

# multiple imputations + Arima function -----------------------------------

fit_arima_MI <- function(sim_list, sim_pars, imputationsnum){
  
  days<-seq(1, 366)
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(days= days,
                                                           GPP = X, 
                                                           light = sim_pars$light, 
                                                           discharge = sim_pars$Q))
  
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
      arimacoefs<-c(tempobj$coef, tempobj$sigma2)
      names(arimacoefs) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
      arimases<-sqrt(diag(vcov(tempobj)))
      names(arimases) <- c("ar1", "intercept", "xreg1", "xreg2")
      name <- paste('imp',seq_along((amelias11sim)[[i]])[[j]],sep='')
      a[[name]] <- arimacoefs
      aa[[name]]<-arimases
    }
    #name1 <- names(amelias11sim)[[i]]
    modelparamlistsim[[i]] <- a
    modelerrorlistsim[[i]] <- aa
  }
  
  
  ### Averages the models together back to 1 model per missing data prop ##
  
  listcoefsessim<-mapply(function(X,Y) {
    list(mi.meld(data.frame(X), data.frame(Y), byrow=FALSE))
  }, X=lapply(modelparamlistsim, function(x) lapply(x, function (x) x[1:4])), 
  Y=modelerrorlistsim)
  # rename coef list
  names(listcoefsessim) <- names(amelias11sim)
  
  # put the mean sigma for each imputation back into the listcoefsessim list
  sigmas_temp <- sapply(modelparamlistsim, function(x) lapply(x, function(x) x[5]))
  sigmas <- apply(sigmas_temp, MARGIN = 2, function(x) mean(as.numeric(x), na.rm = TRUE))
  
  # make return values
  #paramlistsim <- map(listcoefsessim , ~.["q.mi"])
  
  paramlistsim <- lapply(seq(1:16), function(x) 
    matrix(c(listcoefsessim[[x]]$q.mi, sigmas[x]),
           nrow = 1, byrow = TRUE, 
           dimnames = 
             list(c(NULL),c("ar1", "intercept", "xreg1", "xreg2", "sigma")))
  )
  names(paramlistsim) <- names(listcoefsessim)
  
  selistsim <- lapply(listcoefsessim, function(x) x$se.mi)
  
  return(list(paramlistsim, 
              selistsim
  ))
  
}

# Run models with drop missing + arima ------------------------------------

arima_drop_MNAR <- fit_arima_dropmissing(gauss_real_MinMaxMiss[[CurSim]]$y,gauss_real_MinMaxMiss[[CurSim]]$sim_params)

## formatting for figure

names(arima_drop_MNAR) <- names(gauss_real_MinMaxMiss[[CurSim]][["y"]])

modeldropparamlist<-purrr::map(arima_drop_MNAR , ~.["arima_pars"])
modeldropSElist<-purrr::map(arima_drop_MNAR , ~.["arima_errors"])

modeldropparamlist2 <- lapply(modeldropparamlist, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElist2 <- lapply(modeldropSElist, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdf <- map_df(modeldropparamlist2, ~as.data.frame(.x), .id="missingprop_autocor")
modeldropSEdf <- map_df(modeldropSElist2, ~as.data.frame(.x), .id="missingprop_autocor")

modeldropdf<-modeldropparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  select(missingprop_autocor, ar1, intercept, light, discharge, sigma)  %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion")

modeldropSEdf<-modeldropSEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2)%>%   select(missingprop_autocor, ar1, intercept, light, discharge) %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion")

## long form ##

paramdroplong <- gather(modeldropdf, param, value, ar1:sigma, factor_key=TRUE)
paramdroplong$missingnessVersion <- rep.int(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"), times = 5)

paramdropSElong <- gather(modeldropSEdf, param, SE, ar1:discharge, factor_key=TRUE)
paramdropSElong$missingnessVersion <- rep.int(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"), times = 4)

paramdroplong2 <- full_join(paramdroplong, paramdropSElong)


# Run models w/ Kalman filter + arima fxn ---------------------------------

arima_kalman_MNAR<- fit_arima_Kalman(gauss_real_MinMaxMiss[[CurSim]]$y,gauss_real_MinMaxMiss[[CurSim]]$sim_params)

## pull out and label what we need ###

names(arima_kalman_MNAR) <- names(gauss_real_MinMaxMiss[[CurSim]][["y"]])

modelNAparamlist<-purrr::map(arima_kalman_MNAR , ~.["arima_pars"])
modelNASElist<-purrr::map(arima_kalman_MNAR , ~.["arima_errors"])

modelNAparamlist2 <- lapply(modelNAparamlist, function(x) as.data.frame(do.call(rbind, x)))
modelNASElist2 <- lapply(modelNASElist, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdf <- map_df(modelNAparamlist2, ~as.data.frame(.x), .id="missingprop_autocor")
modelNASEdf <- map_df(modelNASElist2, ~as.data.frame(.x), .id="missingprop_autocor")


modelNAdf<-modelNAparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(missingprop_autocor, ar1, intercept, light, discharge, sigma) %>% mutate(missingness="MNAR") %>% mutate(type="Kalman filter")

modelNASEdf<-modelNASEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(missingprop_autocor, ar1, intercept, light, discharge) %>% mutate(missingness="MNAR") %>% mutate(type="Kalman filter")

## long form ##

paramNAlong <- gather(modelNAdf, param, value, ar1:sigma, factor_key=TRUE)
paramNAlong$missingnessVersion <- rep.int(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"), times = 5)

paramNASElong <- gather(modelNASEdf, param, SE, ar1:discharge, factor_key=TRUE)
paramNASElong$missingnessVersion <- rep.int(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"), times = 4)

paramNAlong2<-full_join(paramNAlong, paramNASElong)


# Run models w/ Multiple Imputations --------------------------------------

arima_mi_MNAR <-  fit_arima_MI(gauss_real_MinMaxMiss[[CurSim]]$y,gauss_real_MinMaxMiss[[CurSim]]$sim_params, imputationsnum=5)

##pulls out parameters and ses ##

paramlistsim<-arima_mi_MNAR[[1]]

selistsim<-arima_mi_MNAR[[2]]

avgparamdf <- map_df(paramlistsim, ~as.data.frame(.x), .id="missingprop_autocor")
avglSEdf <- map_df(selistsim, ~as.data.frame(.x), .id="missingprop_autocor")


avgparamdf2 <- avgparamdf %>% 
  dplyr::rename(light=xreg1, discharge=xreg2) %>%  
  select(missingprop_autocor,  ar1, intercept, light, discharge, sigma)  %>% 
  mutate(missingness="MNAR") %>% mutate(type="Multiple imputations")

avglSEdf2 <-avglSEdf  %>% dplyr::rename(light=xreg1, discharge=xreg2) %>%  
  select(missingprop_autocor,  ar1, intercept, light, discharge)   %>% 
  mutate(missingness="MNAR") %>% mutate(type="Multiple imputations")


paramMIlong <- gather(avgparamdf2, param, value, ar1:sigma, factor_key=TRUE)
paramMIlong$missingnessVersion <- rep.int(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"), times = 5)

paramMISElong <- gather(avglSEdf2, param, SE, ar1:discharge, factor_key=TRUE)
paramMISElong$missingnessVersion <- rep.int(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"), times = 4)

paramMIlong2 <-full_join(paramMIlong,paramMISElong)



# Combine output data and save --------------------------------------------

paramarimaall<-rbind(paramdroplong2, paramNAlong2, paramMIlong2)

Output <- matrix(data=NA, nrow=nrow(paramarimaall), ncol=ncol(paramarimaall))

# Save the results of the current script's simulation to the appropriate column of output
Output<- paramarimaall


# Write the output to the folder which will contain all output files as separate csv
#    files with a single line of data.
write.csv(Output, file = OutFile, row.names = FALSE)

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

