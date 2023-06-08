# Load packages
library(here)
library(stats)
library(forecast)
library(tidyverse)
library(lubridate)
library(Amelia)



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
          tempobj=Arima(amelias11sim[[i]][[j]]$GPP, order = c(1,0,0), xreg = matrix(c(simmissingdf [[j]][["light"]],simmissingdf [[j]][["discharge"]]), ncol = 2))
          arimacoefs<-tempobj$coef
          arimases<-sqrt(diag(vcov(tempobj)))
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



#example code using this function:
#gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
#GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

#arima_MI <- fit_arima_MI(GPP_sim_MAR$y,GPP_sim_MAR$sim_params, imputationsnum=5)
