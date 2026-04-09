
library(here)
library(tidyverse)
source(here("Functions/ricker_MI_function.R"))
source(here("Functions/ricker_drop_function.R"))
source(here("Functions/ricker_count_MCMC.R"))
source(here("Functions/ricker_count_EM.R"))
source(here("Functions/ricker_count_likelihood_functions.R"))

forecast_rmse=function(ralpha,trueTS){
  r=ralpha[1]
  alpha=ralpha[2]
  l_ts=length(trueTS)
  # we assume that the 1st value in trueTS is the given starting value for all
  pred_TS=numeric(l_ts-1) 
  pred_TS[1]=trueTS[1]
  for(i in 1:(l_ts-1)){
    pred_TS[i+1]=pred_TS[i]*exp(r-alpha*pred_TS[i])
  }
  
  RMSE=sqrt(mean((pred_TS[2:l_ts] - trueTS[2:l_ts])^2))
  return(RMSE)
}

data_MNAR=readRDS(here("data/missingDatasets/pois_real_minMaxMiss.rds"))

data_MAR=readRDS(here("data/missingDatasets/pois_real_randMiss.rds"))

# Data comes in 10 reps, each one from a different "patch"
# 1) Need to carefully do the offsetting for modeling
# 2) data MAR is longer (has autocorrelation) whereas data MNAR is shorter
# 3) Our goal is to try leave one out prediction for the reps PLUS compare r/alpha estimates with those for non-missing timeseries
# 4) How did Dusty deal with 2 vs 3 days? for now he ignored them and treated them as regular intervals

# STEPS
# 1) careful offsetting
# 2) model training for FULL and for missing sets
# 3) model estimation and prediction
# 4) iterate for LOO
# 5) iterate for all data sets

# Issue with preslicing data with the several replicates model:
# Need to input presliced data so that we never predict the first time point of a timeseries
# However, this may mess up drop only- we don't want to drop NAs that are between time series
# Need to add extra functionality to take presliced data and deal with its complications
# Methods this will mess with and what to fix:
# DROP- don't link together points between different time series
# CC- no issue- all these cases are removed
# MI- feed data sets separately and then run model outside
# EM- Unsure
# DA- Unsure

#offset:

data_MAR$y$patch # use this to offset
patches=levels(data_MAR$y$patch)

autocor_result=numeric(0)
autocor_i=numeric(0)
autocor_result2=numeric(0)
missing_result=numeric(0)
drop_fits=list()
cc_fits=list()
DA_fits=list()
EM_fits=list()
#MI_fits=list()
noMiss_fits=list()
rmse_drop=numeric(0)
rmse_cc=numeric(0)
#rmse_MI=numeric(0)
rmse_noMiss=numeric(0)
rmse_DA=numeric(0)
rmse_EM=numeric(0)
LOO_rep=numeric(0)

start_time=Sys.time()
for(j in 1:150){ # iterate over autocor
  print(paste("we are starting j=",j))
  print(paste("it has been",Sys.time()-start_time))
  offset_ref=data_MAR[[j]]$y$patch # this is the patches reference
  patches=levels(offset_ref) # this is a vector of the different patches
  
  for(k in 4:6){ # iterate over missingness
    print(paste("we are starting k=",k))
    print(paste("it has been",Sys.time()-start_time))
    all_patch_dat=matrix(nrow=0,ncol=3)
   
    for(i in 1:length(patches)){ # within each random set separate the 10 patches for alignment
      
      patch_i=which(data_MAR[[j]]$y$patch==patches[i])
      
      #print(data_MAR[[j]]$y[[k]][patch_i])
      #print(length(data_MAR[[j]]$y[[k]][patch_i]))
      # offset data
      yt=c(data_MAR[[j]]$y[[k]][patch_i],NA)
      ytm1=c(NA,data_MAR[[j]]$y[[k]][patch_i])
      # remove leading and trailing NAs
      if(is.na(ytm1[1])){
        cut1=min(which(!is.na(ytm1)))
      } else {
        cut1=1
      }
      if(is.na(yt[length(yt)])){
        cut2=max(which(!is.na(yt)))
      } else {
        cut2=(length(ytm1)-1)
      }
      yt=yt[cut1:cut2]
      ytm1=ytm1[cut1:cut2]
      
      patchv=rep(patches[i],length(ytm1))
      # combine offset data
      patch_dat=cbind(ytm1,yt,patchv)
      
      all_patch_dat=rbind(all_patch_dat,patch_dat)
    }
    # offset data for all patches
    all_patch_dat=as.data.frame(all_patch_dat)
    
    for(l in 1:length(patches)){ # iterate over patches for LOO length(patches)
      print(paste("we are starting l=",l))
      print(paste("it has been",Sys.time()-start_time))
      # create LOO data sets by excluding patches[l]
      LOO_patch_dat=all_patch_dat[-which(all_patch_dat$patchv==patches[l]),]
      LOO_patch_dat$ytm1=as.numeric(LOO_patch_dat$ytm1)
      LOO_patch_dat$yt=as.numeric(LOO_patch_dat$yt)
      
      # fit actual model
      fits_cc=fit_ricker_cc(LOO_patch_dat,off_patch = T,fam="neg_binom")
      fits_drop=fit_ricker_drop(LOO_patch_dat,off_patch = T,fam="neg_binom",patch_col = "patchv")
      fits_EM=fit_ricker_EM(LOO_patch_dat,fam="neg_binom",off_patch = T)
      # still need to work on this
      #saveRDS(LOO_patch_dat,file="/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/sample_data_set.rds")
      set.seed(1340598)
      #LOO_patch_dat=readRDS(file="sample_data_set.rds")
      fits_DA=tryCatch({
        withTimeout({fit_ricker_DA(LOO_patch_dat,fam="neg_binom",off_patch = T,
                                      priors_list = list(
                                        m_r = 1,
                                        sd_r = 0.3,
                                        m_lalpha = -3,
                                        sd_lalpha = 0.8,
                                        m_lpsi = 5,
                                        sd_lpsi = 2.5
                                      ))}, timeout=300)
      }, TimeoutException=function(ex){
          message("The fit timed out and was aborted")
        return(list(NA,reason="excessive compute time"))
        },
        error=function(e) {
          message("An error occurred: ", conditionMessage(e))
          return(list(NA, reason="model fitting error"))
        }
      )
      
      

      
      patches2=unique(LOO_patch_dat$patchv)
      # imputed=ricker_MI(LOO_patch_dat[LOO_patch_dat$patchv==patches2[1],],off_patch = T)
      # for(m in 2:length(patches2)){ # treat each patch individually for MI
      #   imputed1=ricker_MI(LOO_patch_dat[LOO_patch_dat$patchv==patches2[m],],off_patch = T)
      #   imputed=Map(rbind,imputed,imputed1)
      # }
      #fits_MI=fit_ricker_MI(imputed,fam="neg_binom",off_patch = T)
      
      
      # fit_ricker_DA
      
      # do predictions
      true_values=data_MAR[[1]]$y$number[which(data_MAR[[1]]$y$patch==patches[l])]
      # start them all at the first value: 20
      fore_cc=tryCatch({forecast_rmse(ralpha=fits_cc$estim,true_values)},
                       error=function(e){
                         message("forecasting issue here")
                         return(NA)})
      fore_drop=tryCatch({forecast_rmse(ralpha=fits_drop$estim,true_values)},
                         error=function(e){
                           message("forecasting issue here")
                           return(NA)})
      #fore_MI=forecast_rmse(ralpha=fits_MI$estim,true_values)
      fore_DA=tryCatch({forecast_rmse(ralpha=fits_DA$estim,true_values)},
                        error=function(e){
                          message("forecasting issue here")
                          return(NA)})
      fore_EM=tryCatch({forecast_rmse(ralpha=fits_EM$estim,true_values)},
                       error=function(e){
                         message("forecasting issue here")
                         return(NA)})
      
      
      # or fit with full data
      true_values2=data_MAR[[1]]$y$number[-which(data_MAR[[1]]$y$patch==patches[l])]
      # compile into sliced dataframe
      dat_true <- data.frame(
        yt = true_values2[2:length(true_values2)],
        ytm1 = true_values2[1:(length(true_values2) - 1)]
      )
      # remove overlap between patches
      dat_true=dat_true[-seq(21,188,by=21),]
      
      fits_noMiss=fit_ricker_cc(dat_true,off_patch = T,fam="neg_binom")
        
      fore_noMiss=forecast_rmse(ralpha=fits_noMiss$estim,true_values)
      # record results as well as actual missingness and autocorrelation
      
      #index=4
      # extract auto cor and i
      print(j)
      print(names(data_MAR)[j])
      match0 <- regexec("_(\\d+)_i(\\d+)$", names(data_MAR)[j])
      matches <- regmatches(names(data_MAR)[j], match0)
      autocor_result=c(autocor_result,matches[[1]][2])
      print(autocor_i)
      print(matches[[1]][3])
      
      autocor_i=c(autocor_i,matches[[1]][3])
      print(autocor_i)
      match0 <- regexec("_(\\d+\\.\\d+)_autoCorr_(\\d+\\.\\d+)$", names(data_MAR[[j]]$y)[k])
      matches <- regmatches(names(data_MAR[[j]]$y)[k], match0)
      autocor_result2=c(autocor_result2,matches[[1]][2])
      missing_result=c(missing_result,matches[[1]][3])
      drop_fits=c(drop_fits,list(fits_drop))
      cc_fits=c(cc_fits,list(fits_cc))
      DA_fits=c(DA_fits,list(fits_DA))
      EM_fits=c(EM_fits,list(fits_EM))
      #MI_fits=c(MI_fits,list(fits_MI))
      noMiss_fits=c(noMiss_fits,list(fits_noMiss))
      rmse_drop=c(rmse_drop,fore_drop)
      rmse_cc=c(rmse_cc,fore_cc)
      rmse_DA=c(rmse_DA,fore_DA)
      rmse_EM=c(rmse_EM,fore_EM)
      #rmse_MI=c(rmse_MI,fore_MI)
      rmse_noMiss=c(rmse_noMiss,fore_noMiss)
      LOO_rep=c(LOO_rep,l)
      print("nice we are finishing a rep here")
      print(paste("j=",j))
      print(paste("k=",k))
      print(paste("l=",l))
    }
   
    
  }
  print("saving what we have so far and incrementing j")
  result=tibble(autocor_result=autocor_result,
                autocor_i=autocor_i,
                autocor_result2=autocor_result2,
                missing_result=missing_result,
                LOO_rep=LOO_rep,
                drop_fits=drop_fits,
                cc_fits=cc_fits,
                DA_fits=DA_fits,
                EM_fits=EM_fits,
                #MI_fits=MI_fits,
                noMiss_fits=noMiss_fits,
                rmse_drop=rmse_drop,
                rmse_cc=rmse_cc,
                rmse_EM=rmse_EM,
                rmse_DA=rmse_DA,
                #rmse_MI=rmse_MI,
                rmse_noMiss=rmse_noMiss
  )
  #View(result)
  saveRDS(result,here("data/model_results/pois_real_randMiss_ccdropEMDA.rds"))
  
}

result=tibble(autocor_result=autocor_result,
              autocor_i=autocor_i,
              autocor_result2=autocor_result2,
              missing_result=missing_result,
              LOO_rep=LOO_rep,
              drop_fits=drop_fits,
              cc_fits=cc_fits,
              DA_fits=DA_fits,
              EM_fits=EM_fits,
              #MI_fits=MI_fits,
              noMiss_fits=noMiss_fits,
              rmse_drop=rmse_drop,
              rmse_cc=rmse_cc,
              rmse_EM=rmse_EM,
              rmse_DA=rmse_DA,
              #rmse_MI=rmse_MI,
              rmse_noMiss=rmse_noMiss
)
View(result)
saveRDS(result,here("data/model_results/pois_real_randMiss_ccdropEMDA.rds"))


##################################################
##################################################
################ MNAR results ####################
##################################################
##################################################

missing_i=numeric(0)
missing_result=numeric(0)
drop_fits=list()
cc_fits=list()
DA_fits=list()
EM_fits=list()
#MI_fits=list()
noMiss_fits=list()
rmse_drop=numeric(0)
rmse_cc=numeric(0)
rmse_DA=numeric(0)
rmse_EM=numeric(0)
#rmse_MI=numeric(0)
rmse_noMiss=numeric(0)
LOO_rep=numeric(0)

offset_ref=data_MNAR$y$patch # this is the patches reference
patches=levels(offset_ref) # this is a vector of the different patches
  for(k in 19:length(data_MNAR$y)){ # iterate over datasets in data_MNAR
    
    all_patch_dat=matrix(nrow=0,ncol=3)
    for(i in 1:length(patches)){ # within each random set separate the 10 patches for alignment
      
      patch_i=which(data_MNAR$y$patch==patches[i])
      
      print(data_MNAR$y[[k]][patch_i])
      print(length(data_MNAR$y[[k]][patch_i]))
      if(all(is.na(data_MNAR$y[[k]][patch_i]))){
        yt=NULL
        ytm1=NULL
      } else {
        # offset data
        yt=c(data_MNAR$y[[k]][patch_i],NA)
        ytm1=c(NA,data_MNAR$y[[k]][patch_i])
        
        # remove leading and trailing NAs
        
        if(is.na(ytm1[1])){
          cut1=min(which(!is.na(ytm1)))
        } else {
          cut1=1
        }
        if(is.na(yt[length(yt)])){
          cut2=max(which(!is.na(yt)))
        } else {
          cut2=(length(ytm1)-1)
        }
        yt=yt[cut1:cut2]
        ytm1=ytm1[cut1:cut2]
      }
      
      
      patchv=rep(patches[i],length(ytm1))
      # combine offset data
      patch_dat=cbind(ytm1,yt,patchv)
      all_patch_dat=rbind(all_patch_dat,patch_dat)
    }
    # offset data for all patches
    all_patch_dat=as.data.frame(all_patch_dat)
    
    for(l in 1:length(unique(all_patch_dat$patchv))){ # iterate over patches for LOO length(patches)
      # create LOO data sets by excluding patches[l]
      LOO_patch_dat=all_patch_dat[-which(all_patch_dat$patchv==unique(all_patch_dat$patchv)[l]),]
      LOO_patch_dat$ytm1=as.numeric(LOO_patch_dat$ytm1)
      LOO_patch_dat$yt=as.numeric(LOO_patch_dat$yt)
      # fit actual model
      fits_cc=fit_ricker_cc(LOO_patch_dat,off_patch = T,fam="neg_binom")
      fits_drop=fit_ricker_drop(LOO_patch_dat,off_patch = T,fam="neg_binom",patch_col = "patchv")
      fits_EM=fit_ricker_EM(LOO_patch_dat,fam="neg_binom",off_patch = T)
      set.seed(1340598)
      #LOO_patch_dat=readRDS(file="sample_data_set.rds")
      fits_DA=tryCatch({fit_ricker_DA(LOO_patch_dat,fam="neg_binom",off_patch = T,
                                      priors_list = list(
                                        m_r = 1,
                                        sd_r = 0.3,
                                        m_lalpha = -3,
                                        sd_lalpha = 0.8,
                                        m_lpsi = 5,
                                        sd_lpsi = 2.5
                                      ))}, 
                       error=function(e) {
                         message("An error occurred: ", conditionMessage(e))
                         return(list(NA, reason="model fitting error"))
                       }
      )
      patches2=unique(LOO_patch_dat$patchv)
      # imputed=ricker_MI(LOO_patch_dat[LOO_patch_dat$patchv==patches2[1],],off_patch = T)
      # for(m in 2:length(patches2)){ # treat each patch individually for MI
      #   imputed1=ricker_MI(LOO_patch_dat[LOO_patch_dat$patchv==patches2[m],],off_patch = T)
      #   imputed=Map(rbind,imputed,imputed1)
      # }
      #fits_MI=fit_ricker_MI(imputed,fam="neg_binom",off_patch = T)
      
      # fit_ricker_EM
      # fit_ricker_DA
      
      # do predictions
      true_values=data_MNAR$y$number[which(data_MNAR$y$patch==patches[l])]
      # start them all at the first value: 20
      fore_cc=tryCatch({forecast_rmse(ralpha=fits_cc$estim,true_values)},
                       error=function(e){
                         message("forecasting issue here")
                         return(NA)})
      fore_drop=tryCatch({forecast_rmse(ralpha=fits_drop$estim,true_values)},
                         error=function(e){
                           message("forecasting issue here")
                           return(NA)})
      #fore_MI=forecast_rmse(ralpha=fits_MI$estim,true_values)
      fore_DA=tryCatch({forecast_rmse(ralpha=fits_DA$estim,true_values)},
                       error=function(e){
                         message("forecasting issue here")
                         return(NA)})
      fore_EM=tryCatch({forecast_rmse(ralpha=fits_EM$estim,true_values)},
                       error=function(e){
                         message("forecasting issue here")
                         return(NA)})
      #fore_MI=forecast_rmse(ralpha=fits_MI$estim,true_values)
      
      
      
      # or fit with full data
      true_values2=data_MNAR$y$number[-which(data_MNAR$y$patch==patches[l])]
      # compile into sliced dataframe
      dat_true <- data.frame(
        yt = true_values2[2:length(true_values2)],
        ytm1 = true_values2[1:(length(true_values2) - 1)]
      )
      # remove overlap between patches
      dat_true=dat_true[-seq(21,188,by=21),]
      fits_noMiss=fit_ricker_cc(dat_true,off_patch = T,fam="neg_binom")
      
      fore_noMiss=forecast_rmse(ralpha=fits_noMiss$estim,true_values)
      
      # record results as well as actual missingness and autocorrelation
      
      #index=4
      # extract auto cor and i
      print(k)
      print(names(data_MNAR$y)[k])
      match0 <- regexec("_(\\d+\\.\\d+)_i(\\d+)$", names(data_MNAR$y)[k])
      matches <- regmatches(names(data_MNAR$y)[k], match0)
      missing_result=c(missing_result,matches[[1]][2])
      missing_i=c(missing_i,matches[[1]][3])
    
      
      drop_fits=c(drop_fits,list(fits_drop))
      cc_fits=c(cc_fits,list(fits_cc))
      DA_fits=c(DA_fits,list(fits_DA))
      EM_fits=c(EM_fits,list(fits_EM))
      #MI_fits=c(MI_fits,list(fits_MI))
      noMiss_fits=c(noMiss_fits,list(fits_noMiss))
      rmse_drop=c(rmse_drop,fore_drop)
      rmse_cc=c(rmse_cc,fore_cc)
      rmse_DA=c(rmse_DA,fore_DA)
      rmse_EM=c(rmse_EM,fore_EM)
      #rmse_MI=c(rmse_MI,fore_MI)
      rmse_noMiss=c(rmse_noMiss,fore_noMiss)
      LOO_rep=c(LOO_rep,l)
    }
    
    
  }
  


result=tibble(missing_result=missing_result,
              missing_i=missing_i,
              LOO_rep=LOO_rep,
              drop_fits=drop_fits,
              cc_fits=cc_fits,
              DA_fits=DA_fits,
              EM_fits=EM_fits,
              #MI_fits=MI_fits,
              noMiss_fits=noMiss_fits,
              rmse_drop=rmse_drop,
              rmse_cc=rmse_cc,
              rmse_DA=rmse_DA,
              rmse_EM=rmse_EM,
              #rmse_MI=rmse_MI,
              rmse_noMiss=rmse_noMiss
)
View(result)
saveRDS(result,here("data/model_results/pois_real_minMaxMiss_dropccEMDA_secondhalf.rds"))

