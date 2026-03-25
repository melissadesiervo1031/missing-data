# Code to run the DA extinctions locally
library(here)
library(tidyverse)


dat_flat=readRDS(here("data/model_results/ricker_Sim_reruns/extinctSets.rds"))

# trim off those terminal zeros that trouble the DA function
dat_flat_trim <- lapply(dat_flat, function(x) x[-length(x)]) 


# do the results in chunks of 1000, in case anything goes wrong, so as not to waste all the compute time
for(i in 1:12){
  if(i<12){
    dat_flat_part=dat_flat_trim[(1000*(i-1)+1):(1000*i)]
  } else {
    dat_flat_part=dat_flat_trim[(1000*(i-1)+1):11223]
  }
  
  DA_fits <- lapply(dat_flat_part, fit_ricker_DA)
  saveRDS(DA_fits,here(paste0("data/model_results/ricker_Sim_reruns/DA_reruns/extinctDA",i,".rds")))
}


