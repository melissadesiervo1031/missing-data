
library(here)
library(tidyverse)
source(here("Functions/ricker_MI_function.R"))
source(here("Functions/ricker_drop_function.R"))

#in_args=c("1", "4", "1","1","1492")
in_args <- commandArgs(trailingOnly = T)
cat(in_args)

j_cur=as.numeric(in_args[1])
k_cur=as.numeric(in_args[2])
l_cur=as.numeric(in_args[3])
set.seed(in_args[4])

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
MI_fits=list()
rmse_MI=numeric(0)
LOO_rep=numeric(0)

for(j in j_cur){ # iterate over autocor
  offset_ref=data_MAR[[j]]$y$patch # this is the patches reference
  patches=levels(offset_ref) # this is a vector of the different patches
  
  for(k in k_cur){ # iterate over missingness
    
    all_patch_dat=matrix(nrow=0,ncol=3)
    for(i in 1:length(patches)){ # within each random set separate the 10 patches for alignment
      
      patch_i=which(data_MAR[[j]]$y$patch==patches[i])
      
      print(data_MAR[[j]]$y[[k]][patch_i])
      print(length(data_MAR[[j]]$y[[k]][patch_i]))
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
    
    for(l in l_cur){ # iterate over patches for LOO length(patches)
      # create LOO data sets by excluding patches[l]
      LOO_patch_dat=all_patch_dat[-which(all_patch_dat$patchv==patches[l]),]
      LOO_patch_dat$ytm1=as.numeric(LOO_patch_dat$ytm1)
      LOO_patch_dat$yt=as.numeric(LOO_patch_dat$yt)
      # fit actual model

      patches2=unique(LOO_patch_dat$patchv)
      imputed=ricker_MI(LOO_patch_dat[LOO_patch_dat$patchv==patches2[1],],off_patch = T)
      for(m in 2:length(patches2)){ # treat each patch individually for MI
        if(length(which(is.na(LOO_patch_dat[LOO_patch_dat$patchv==patches2[m],])))==0){
          imputed_0=LOO_patch_dat[LOO_patch_dat$patchv==patches2[m],1:2]
          imputed1=list()
          for(n_imp in 1:5){
            imputed1[[n_imp]]=imputed_0
          }
          
        } else {
          imputed1=ricker_MI(LOO_patch_dat[LOO_patch_dat$patchv==patches2[m],],off_patch = T)
        }
        print(paste("we are at m=", m))
        print('check this one')
        print(imputed1)
        imputed=Map(rbind,imputed,imputed1)
      }
      fits_MI=fit_ricker_MI(imputed,fam="neg_binom",off_patch = T)
      
      # fit_ricker_EM
      # fit_ricker_DA
      
      # do predictions
      true_values=data_MAR[[1]]$y$number[which(data_MAR[[1]]$y$patch==patches[l])]
      # start them all at the first value: 20

      fore_MI=forecast_rmse(ralpha=fits_MI$estim,true_values)

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


      MI_fits=c(MI_fits,list(fits_MI))


      rmse_MI=c(rmse_MI,fore_MI)

      LOO_rep=c(LOO_rep,l)
    }
    
    
  }
  
}

result=tibble(autocor_result=autocor_result,
              autocor_i=autocor_i,
              autocor_result2=autocor_result2,
              LOO_rep=LOO_rep,
              MI_fits=MI_fits,
              rmse_MI=rmse_MI,
)

saveRDS(result,paste0("/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/res_MI_MAR_real/result_",j,"_",k,"_",l,".rds"))



