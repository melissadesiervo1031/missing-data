
library(here)
library(tidyverse)
source(here("Functions/ricker_MI_function.R"))
source(here("Functions/ricker_drop_function.R"))

#in_args=c("4","1","1492")
in_args <- commandArgs(trailingOnly = T)
cat(in_args)


k_cur=as.numeric(in_args[1])
l_cur=as.numeric(in_args[2])
set.seed(as.numeric(in_args[3]))

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

missing_i=numeric(0)
missing_result=numeric(0)
MI_fits=list()
rmse_MI=numeric(0)
LOO_rep=numeric(0)

offset_ref=data_MNAR$y$patch # this is the patches reference
patches=levels(offset_ref) # this is a vector of the different patches
for(k in k_cur){ # iterate over datasets in data_MNAR
  
  all_patch_dat=matrix(nrow=0,ncol=3)
  for(i in 1:length(patches)){ # within each random set separate the 10 patches for alignment
    
    patch_i=which(data_MNAR$y$patch==patches[i])
    
    print(data_MNAR$y[[k]][patch_i])
    print(length(data_MNAR$y[[k]][patch_i]))
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
    if(length(which(is.na(LOO_patch_dat[LOO_patch_dat$patchv==patches2[1],])))==0){
      imputed_0=LOO_patch_dat[LOO_patch_dat$patchv==patches2[1],1:2]
      imputed=list()
      for(n_imp in 1:5){
        imputed[[n_imp]]=imputed_0
      }
      
    } else {
      imputed=ricker_MI(LOO_patch_dat[LOO_patch_dat$patchv==patches2[1],],off_patch = T)
    }
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
      imputed=Map(rbind,imputed,imputed1)
    }
    fits_MI=fit_ricker_MI(imputed,fam="neg_binom",off_patch = T)
    
    
    # do predictions
    true_values=data_MNAR$y$number[which(data_MNAR$y$patch==patches[l])]
    # start them all at the first value: 20

    fore_MI=forecast_rmse(ralpha=fits_MI$estim,true_values)

    # record results as well as actual missingness and autocorrelation
    
    #index=4
    # extract auto cor and i
    print(k)
    print(names(data_MNAR$y)[k])
    match0 <- regexec("_(\\d+\\.\\d+)_i(\\d+)$", names(data_MNAR$y)[k])
    matches <- regmatches(names(data_MNAR$y)[k], match0)
    missing_result=c(missing_result,matches[[1]][2])
    missing_i=c(missing_i,matches[[1]][3])
    
    MI_fits=c(MI_fits,list(fits_MI))
    rmse_MI=c(rmse_MI,fore_MI)
    LOO_rep=c(LOO_rep,l)
  }
  
  
}



result=tibble(missing_result=missing_result,
              missing_i=missing_i,
              LOO_rep=LOO_rep,
              MI_fits=MI_fits,
              rmse_MI=rmse_MI
)

View(result)
saveRDS(result,paste0("/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/res_MI_MNAR_real/result_",k,"_",l,".rds"))




