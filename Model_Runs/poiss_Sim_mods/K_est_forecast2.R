# K estimation and new forecasts
library(here)
library(tidyverse)
source(here("Simulations/ricker_data_sims.R"))

# copied from modelruns_ricker.R
forecast_rmse=function(ralpha,trueTS){
  r=ralpha[1]
  alpha=ralpha[2]
  l_ts=length(trueTS)
  # we assume that the 1st value in trueTS is the given starting value for all
  pred_TS=numeric(l_ts-1) 
  pred_TS[1]=trueTS[1]
  for(i in 1:(l_ts-1)){
    pred_TS[i+1]=pred_TS[i]*exp(r-alpha*pred_TS[i])
    pred_TS[i+1]=max(pred_TS[i+1],0)
  }
  
  RMSE=sqrt(mean((pred_TS[2:l_ts] - trueTS[2:l_ts])^2))
  return(RMSE)
}

#############################################################################
#############################################################################
# MCAR A and B

ricDat_tempA <- readRDS("./data/model_results/RickerA_resultTableRev1.rds")
ricDat_tempB <- readRDS("./data/model_results/RickerB_resultTableRev1.rds")

ricDat_temp <- rbind(ricDat_tempA, ricDat_tempB)
dat0=readRDS(here("data/ricker_0miss_datasets.rds"))

ricDat_temp$K_true=ricDat_temp$r/ricDat_temp$alpha
ricDat_temp$K_drop=NA
ricDat_temp$K_cc=NA
ricDat_temp$K_EM=NA
ricDat_temp$K_DA=NA
ricDat_temp$K_drop_err=NA
ricDat_temp$K_cc_err=NA
ricDat_temp$K_EM_err=NA
ricDat_temp$K_DA_err=NA

for(i in 1:nrow(ricDat_temp)){ # calculate and record estimated r/alpha for each method
  print(paste("we are at i=",i))
  K_drop_est=tryCatch({
    ricDat_temp$drop_fits[[i]]$estim[1]/ricDat_temp$drop_fits[[i]]$estim[2]
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  if(length(K_drop_est)==0){K_drop_est=NA}
  ricDat_temp$K_drop[i]=K_drop_est
  ricDat_temp$K_drop_err[i]=K_drop_est-ricDat_temp$K_true[i]
  K_cc_est=tryCatch({
      ricDat_temp$cc_fits[[i]]$estim[1]/ricDat_temp$cc_fits[[i]]$estim[2]
    },
    error=function(e){
      message("an error occurred",e$message)
      return(NA)},
    warning=function(w){
      message("a warning occurred",w$message)
    }
  )

  if(length(K_cc_est)==0){K_cc_est=NA}
  ricDat_temp$K_cc[i]=K_cc_est
  ricDat_temp$K_cc_err[i]=K_cc_est-ricDat_temp$K_true[i]
  K_EM_est=tryCatch({
    ricDat_temp$EM_fits[[i]]$estim[1]/ricDat_temp$EM_fits[[i]]$estim[2]
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  
  if(length(K_EM_est)==0){K_EM_est=NA}
  ricDat_temp$K_EM[i]=K_EM_est
  ricDat_temp$K_EM_err[i]=K_EM_est-ricDat_temp$K_true[i]
  K_DA_est=tryCatch({
    ricDat_temp$DA_fits[[i]]$estim[1]/ricDat_temp$DA_fits[[i]]$estim[2]
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  if(length(K_DA_est)==0){K_DA_est=NA}
  ricDat_temp$K_DA[i]=K_DA_est
  ricDat_temp$K_DA_err[i]=K_DA_est-ricDat_temp$K_true[i]
  
}

# tabulate error in K across methods and missingness
dat1=tapply(ricDat_temp$K_drop_err,ricDat_temp$propMiss,FUN=median,na.rm=T)
plot(as.numeric(names(dat1)),dat1,type="b",xlab="Proportion missing",
     ylab="median error in r/alpha",ylim=c(-7,5),main="MCAR simulated data")
dat2=tapply(ricDat_temp$K_cc_err,ricDat_temp$propMiss,FUN=median,na.rm=T)
lines(as.numeric(names(dat2)),dat2,type="b",col="red")
dat3=tapply(ricDat_temp$K_EM_err,ricDat_temp$propMiss,FUN=median,na.rm=T)
lines(as.numeric(names(dat3)),dat3,type="b",col="blue")
dat4=tapply(ricDat_temp$K_DA_err,ricDat_temp$propMiss,FUN=median,na.rm=T)
lines(as.numeric(names(dat4)),dat4,type="b",col="green")

dat1=tapply(ricDat_temp$K_drop_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.75)
lines(as.numeric(names(dat1)),dat1,type="l",col="black",lty=2)
dat2=tapply(ricDat_temp$K_cc_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.75)
lines(as.numeric(names(dat2)),dat2,type="l",col="red",lty=2)
dat3=tapply(ricDat_temp$K_EM_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.75)
lines(as.numeric(names(dat3)),dat3,type="l",col="blue",lty=2)
dat4=tapply(ricDat_temp$K_DA_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.75)
lines(as.numeric(names(dat4)),dat4,type="l",col="green",lty=2)

dat1=tapply(ricDat_temp$K_drop_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.25)
lines(as.numeric(names(dat1)),dat1,type="l",col="black",lty=2)
dat2=tapply(ricDat_temp$K_cc_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.25)
lines(as.numeric(names(dat2)),dat2,type="l",col="red",lty=2)
dat3=tapply(ricDat_temp$K_EM_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.25)
lines(as.numeric(names(dat3)),dat3,type="l",col="blue",lty=2)
dat4=tapply(ricDat_temp$K_DA_err,ricDat_temp$propMiss,FUN=quantile,na.rm=T,probs=0.25)
lines(as.numeric(names(dat4)),dat4,type="l",col="green",lty=2)


text(0.2,5,"drop",col="black")
text(0.2,4,"cc",col="red")
text(0.2,-5,"EM",col="blue")
text(0.2,-6,"DA",col="green")

# new forecasts from small



# in the other data sets we have time series of length 60, and we were fitting to the first 48, leaving 12 points held out
ricker_forecast_sims=list()
for(i in 1:max(ricDat_temp$SimNumber,na.rm=T)){
  ricker_forecast_sims[[i]]=ricker_sim(12, # number of points
                                       r=ricDat_temp$r[which(ricDat_temp$SimNumber==i)[1]], # r from simulation
                                       alpha=ricDat_temp$alpha[which(ricDat_temp$SimNumber==i)[1]], # alpha from simulation
                                       N0=5)
}

ricDat_temp$forecast_2_drop=NA
ricDat_temp$forecast_2_cc=NA
ricDat_temp$forecast_2_EM=NA
ricDat_temp$forecast_2_DA=NA

for(i in 1:nrow(ricDat_temp)){
  print(paste("working on", i))
  simnum=ricDat_temp$SimNumber[i]
  trueTS=ricker_forecast_sims[[simnum]]
  
  # drop
  ricDat_temp$forecast_2_drop[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$drop_fits[[i]]$estim[["r"]],ricDat_temp$drop_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  
  # cc
  ricDat_temp$forecast_2_cc[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$cc_fits[[i]]$estim[["r"]],ricDat_temp$cc_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )

  # EM
  ricDat_temp$forecast_2_EM[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$EM_fits[[i]]$estim[["r"]],ricDat_temp$EM_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  

  # DA
  ricDat_temp$forecast_2_DA[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$DA_fits[[i]]$estim[["r"]],ricDat_temp$DA_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
}

saveRDS(ricDat_temp[1:75000,],"./data/model_results/RickerA_resultTableRev1.rds")
saveRDS(ricDat_temp[75001:15000,],"./data/model_results/RickerB_resultTableRev1.rds")
saveRDS(ricker_forecast_sims,"./data/ricker_forecast_sim.rds")

#############################################################################
#############################################################################
# MCAR MI A and B

ricker_forecast_sims=readRDS("./data/ricker_forecast_sim.rds")
ricDat_MIA<-readRDS("./data/model_results/RickerA_resultTableMIRev1.rds")
ricDat_MIB<-readRDS("./data/model_results/RickerB_resultTableMIRev1.rds")
ricDat_temp=rbind(ricDat_MIA,ricDat_MIB)

ricDat_temp$K_true=ricDat_temp$r/ricDat_temp$alpha
ricDat_temp$K_MI=NA
ricDat_temp$K_MI_err=NA

for(i in 1:nrow(ricDat_temp)){ # calculate and record estimated r/alpha for each method
  print(paste("we are at i=",i))
  K_MI_est=tryCatch({
    as.numeric(ricDat_temp$estim_r[i])/as.numeric(ricDat_temp$estim_alpha[i])
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  if(length(K_MI_est)==0){K_MI_est=NA}
  ricDat_temp$K_MI[i]=K_MI_est
  ricDat_temp$K_MI_err[i]=K_MI_est-ricDat_temp$K_true[i]
 
  
}


ricDat_temp$forecast_2_MI=NA

for(i in 1:nrow(ricDat_temp)){
  print(paste("working on", i))
  simnum=ricDat_temp$SimNumber[i]
  trueTS=ricker_forecast_sims[[simnum]]
  
  ricDat_temp$forecast_2_MI[i]=tryCatch({
    forecast_rmse(c(as.numeric(ricDat_temp$estim_r[i]),as.numeric(ricDat_temp$estim_alpha[i])),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
    return(NA)
  }
  )
  
  
}


saveRDS(ricDat_temp[1:75000,],"./data/model_results/RickerA_resultTableMIRev1.rds")
saveRDS(ricDat_temp[75001:15000,],"./data/model_results/RickerB_resultTableMIRev1.rds")


#############################################################################
#############################################################################
# MNAR 
ricDat_temp=readRDS("./data/model_results/RickerMinMaxMissRev1.rds")
ricker_forecast_sims=readRDS("./data/ricker_forecast_sim.rds")

ricDat_temp$K_true=ricDat_temp$r/ricDat_temp$alpha
ricDat_temp$K_drop=NA
ricDat_temp$K_cc=NA
ricDat_temp$K_EM=NA
ricDat_temp$K_DA=NA
ricDat_temp$K_drop_err=NA
ricDat_temp$K_cc_err=NA
ricDat_temp$K_EM_err=NA
ricDat_temp$K_DA_err=NA

for(i in 1:nrow(ricDat_temp)){ # calculate and record estimated r/alpha for each method
  print(paste("we are at i=",i))
  K_drop_est=tryCatch({
    ricDat_temp$drop_fits[[i]]$estim[1]/ricDat_temp$drop_fits[[i]]$estim[2]
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  if(length(K_drop_est)==0){K_drop_est=NA}
  ricDat_temp$K_drop[i]=K_drop_est
  ricDat_temp$K_drop_err[i]=K_drop_est-ricDat_temp$K_true[i]
  K_cc_est=tryCatch({
    ricDat_temp$cc_fits[[i]]$estim[1]/ricDat_temp$cc_fits[[i]]$estim[2]
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  
  if(length(K_cc_est)==0){K_cc_est=NA}
  ricDat_temp$K_cc[i]=K_cc_est
  ricDat_temp$K_cc_err[i]=K_cc_est-ricDat_temp$K_true[i]
  K_EM_est=tryCatch({
    ricDat_temp$EM_fits[[i]]$estim[1]/ricDat_temp$EM_fits[[i]]$estim[2]
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  
  if(length(K_EM_est)==0){K_EM_est=NA}
  ricDat_temp$K_EM[i]=K_EM_est
  ricDat_temp$K_EM_err[i]=K_EM_est-ricDat_temp$K_true[i]
  K_DA_est=tryCatch({
    ricDat_temp$DA_fits[[i]]$estim[1]/ricDat_temp$DA_fits[[i]]$estim[2]
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  if(length(K_DA_est)==0){K_DA_est=NA}
  ricDat_temp$K_DA[i]=K_DA_est
  ricDat_temp$K_DA_err[i]=K_DA_est-ricDat_temp$K_true[i]
  
}


ricDat_temp$forecast_2_drop=NA
ricDat_temp$forecast_2_cc=NA
ricDat_temp$forecast_2_EM=NA
ricDat_temp$forecast_2_DA=NA

for(i in 1:nrow(ricDat_temp)){
  print(paste("working on", i))
  simnum=ricDat_temp$SimNumber[i]
  trueTS=ricker_forecast_sims[[simnum]]
  
  # drop
  ricDat_temp$forecast_2_drop[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$drop_fits[[i]]$estim[["r"]],ricDat_temp$drop_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  
  # cc
  ricDat_temp$forecast_2_cc[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$cc_fits[[i]]$estim[["r"]],ricDat_temp$cc_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  
  # EM
  ricDat_temp$forecast_2_EM[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$EM_fits[[i]]$estim[["r"]],ricDat_temp$EM_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  
  
  # DA
  ricDat_temp$forecast_2_DA[i]=tryCatch({
    forecast_rmse(c(ricDat_temp$DA_fits[[i]]$estim[["r"]],ricDat_temp$DA_fits[[i]]$estim[["alpha"]]),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
}

saveRDS(ricDat_temp,"./data/model_results/RickerMinMaxMissRev1.rds")

#############################################################################
#############################################################################
# MNAR MI
ricDat_temp=readRDS("./data/model_results/RickerMinMaxMiss_MIRev1.rds")
ricker_forecast_sims=readRDS("./data/ricker_forecast_sim.rds")


ricDat_temp$K_true=ricDat_temp$r/ricDat_temp$alpha
ricDat_temp$K_MI=NA
ricDat_temp$K_MI_err=NA

for(i in 1:nrow(ricDat_temp)){ # calculate and record estimated r/alpha for each method
  print(paste("we are at i=",i))
  K_MI_est=tryCatch({
    as.numeric(ricDat_temp$estim_r[i])/as.numeric(ricDat_temp$estim_alpha[i])
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
  }
  )
  if(length(K_MI_est)==0){K_MI_est=NA}
  ricDat_temp$K_MI[i]=K_MI_est
  ricDat_temp$K_MI_err[i]=K_MI_est-ricDat_temp$K_true[i]
  
  
}


ricDat_temp$forecast_2_MI=NA

for(i in 1:nrow(ricDat_temp)){
  print(paste("working on", i))
  simnum=ricDat_temp$SimNumber[i]
  trueTS=ricker_forecast_sims[[simnum]]
  
  ricDat_temp$forecast_2_MI[i]=tryCatch({
    forecast_rmse(c(as.numeric(ricDat_temp$estim_r[i]),as.numeric(ricDat_temp$estim_alpha[i])),trueTS)
  },
  error=function(e){
    message("an error occurred",e$message)
    return(NA)},
  warning=function(w){
    message("a warning occurred",w$message)
    return(NA)
  }
  )
  
  
}


saveRDS(ricDat_temp,"./data/model_results/RickerMinMaxMiss_MIRev1.rds")

