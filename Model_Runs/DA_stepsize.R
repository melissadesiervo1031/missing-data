################################################################
# This script tries a few different values of the stepsize parameter for the DA method
# 
################################################################

# load libraries
library(tidyverse)
library(here)
library(patchwork)
func_list <- list.files(here("Functions/"), pattern = ".R", full.names = T)
lapply(func_list, source)

in_args <- commandArgs(trailingOnly = T)
cat(in_args)

stepSizeSeq=seq(0.01,0.15,length.out=20)


# establish some global params
set.seed(9621)
n_params <- 100
uns <- c(50, 100, 1000, 5000, 10000)
nns <- length(uns)

# create data
ns <- rep(uns, each = n_params)
r <- rep(runif(n_params, min = 0.5, max = log(2)), nns) 
alpha <- rep(runif(n_params, min = 0.01, max = 0.05), nns)  

y <- mapply(
  ricker_sim,
  n = ns, r = r, alpha = alpha, N0 = 10
)

# add these to a dataframe
dat <- tibble(
  y = y,
  r = r,
  alpha = alpha,
  n = ns
)

# find and remove simulations in which pop went extinct
dat <- dat %>% mutate(
  ext = map_dbl(
    dat$y,
    ~ as.numeric(.x[length(.x)] == 0)
  )
)
dat <- dat %>% filter(
  ext == 0
)

# create datasets with 15% missing
dat <- dat %>% mutate(
  y_miss = map(
    y,
    ~ makeMissing(.x, "random", propMiss = 0.2)
  ) %>% flatten()
)

cat(stepSizeSeq[in_args])

# fit the models DA
dat <- dat %>% mutate(
  estims_full_DA = lapply(
    y,
    FUN = fit_ricker_DA,
    stepsize=stepSizeSeq[in_args]
  ),
  estims_miss_DA = lapply(
    y_miss,
    FUN = fit_ricker_DA
  )
)

# fit the models EM
dat <- dat %>% mutate(
  estims_full_EM = lapply(
    y,
    FUN = fit_ricker_EM
  ),
  estims_miss_EM = lapply(
    y_miss,
    FUN = fit_ricker_EM
  )
)

# fit the models MI
dat <- dat %>% mutate(
  estims_full_MI = lapply(
    y,
    FUN = fit_ricker_MI
  ),
  estims_miss_MI = lapply(
    y_miss,
    FUN = fit_ricker_MI
  )
)

# fit the models drop simple
dat <- dat %>% mutate(
  estims_full_drop = lapply(
    y,
    FUN = fit_ricker_drop
  ),
  estims_miss_drop = lapply(
    y_miss,
    FUN = fit_ricker_drop
  )
)

# fit the models drop cc
dat <- dat %>% mutate(
  estims_full_cc = lapply(
    y,
    FUN = fit_ricker_cc
  ),
  estims_miss_cc = lapply(
    y_miss,
    FUN = fit_ricker_cc
  )
)



rhat_full_EM=numeric(nrow(dat))
ahat_full_EM=numeric(nrow(dat))
rhat_miss_EM=numeric(nrow(dat))
ahat_miss_EM=numeric(nrow(dat))

# rhat_full_MI=numeric(nrow(dat))
# ahat_full_MI=numeric(nrow(dat))
rhat_miss_MI=numeric(nrow(dat))
ahat_miss_MI=numeric(nrow(dat))

rhat_full_DA=numeric(nrow(dat))
ahat_full_DA=numeric(nrow(dat))
rhat_miss_DA=numeric(nrow(dat))
ahat_miss_DA=numeric(nrow(dat))

rhat_full_drop=numeric(nrow(dat))
ahat_full_drop=numeric(nrow(dat))
rhat_miss_drop=numeric(nrow(dat))
ahat_miss_drop=numeric(nrow(dat))

rhat_full_cc=numeric(nrow(dat))
ahat_full_cc=numeric(nrow(dat))
rhat_miss_cc=numeric(nrow(dat))
ahat_miss_cc=numeric(nrow(dat))

for(i in 1:nrow(dat)){
  
  rhat_full_EM[i]=dat$estims_full_EM[[i]]$estim[1]
  ahat_full_EM[i]=dat$estims_full_EM[[i]]$estim[2]
  rhat_miss_EM[i]=dat$estims_miss_EM[[i]]$estim[1]
  ahat_miss_EM[i]=dat$estims_miss_EM[[i]]$estim[2]
  
  
  rhat_miss_MI[i]=dat$estims_miss_MI[[i]]$estim[1]
  ahat_miss_MI[i]=dat$estims_miss_MI[[i]]$estim[2]
  
  rhat_full_DA[i]=dat$estims_full_DA[[i]]$estim[1]
  ahat_full_DA[i]=dat$estims_full_DA[[i]]$estim[2]
  rhat_miss_DA[i]=dat$estims_miss_DA[[i]]$estim[1]
  ahat_miss_DA[i]=dat$estims_miss_DA[[i]]$estim[2]
  
  rhat_full_cc[i]=dat$estims_full_cc[[i]]$estim[1]
  ahat_full_cc[i]=dat$estims_full_cc[[i]]$estim[2]
  rhat_miss_cc[i]=dat$estims_miss_cc[[i]]$estim[1]
  ahat_miss_cc[i]=dat$estims_miss_cc[[i]]$estim[2]
  
  rhat_full_drop[i]=dat$estims_full_drop[[i]]$estim[1]
  ahat_full_drop[i]=dat$estims_full_drop[[i]]$estim[2]
  rhat_miss_drop[i]=dat$estims_miss_drop[[i]]$estim[1]
  ahat_miss_drop[i]=dat$estims_miss_drop[[i]]$estim[2]
}


res=cbind(dat$r,dat$alpha,dat$n,rhat_full_EM,ahat_full_EM,rhat_miss_EM,ahat_miss_EM,rhat_miss_MI,ahat_miss_MI,rhat_full_DA,ahat_full_DA,rhat_miss_DA,ahat_miss_DA,rhat_full_cc,ahat_full_cc,rhat_miss_cc,ahat_miss_cc,rhat_full_drop,ahat_full_drop,rhat_miss_drop,ahat_miss_drop)

colnames(res)[1:3]=c("r","alpha","n")

newnames=paste0("rel_bias_",colnames(res))

for(i in 4:21){
  if(i%%2==0){
    res=cbind(res,(res[,i]-res[,1])/res[,1])
    colnames(res)[ncol(res)]=newnames[i]
  } else {
    res=cbind(res,(res[,i]-res[,2])/res[,2])
    colnames(res)[ncol(res)]=newnames[i]
  }
  
}


res_sum=c(50,100,1000,5000,10000)

for(i in 22:39){
  res_sum=cbind(res_sum,tapply(res[,i],res[,3],mean))
  colnames(res_sum)[ncol(res_sum)]=paste0(colnames(res)[i],"mean")
  res_sum=cbind(res_sum,tapply(res[,i],res[,3],sd))
  colnames(res_sum)[ncol(res_sum)]=paste0(colnames(res)[i],"sd")
  
}


n_adj = log(res_sum[,1]) + log(0.2) * (1 - 1:nrow(.) %% 2)

res_sum=cbind(res_sum, n_adj_miss=log(res_sum[,1]) + log(0.2))
res_sum=cbind(res_sum, n_adj_full=log(res_sum[,1]))

res_sum=as.data.frame(res_sum)


# create plots
r <- ggplot() +
  geom_line(data = res_sum, aes(x = n_adj_miss, y = rel_bias_rhat_miss_EMmean, color = "EM"))+
  geom_point(data = res_sum, aes(x = n_adj_miss, y = rel_bias_rhat_miss_EMmean, color = "EM"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_rhat_miss_EMmean - rel_bias_rhat_miss_EMsd, ymax = rel_bias_rhat_miss_EMmean + rel_bias_rhat_miss_EMsd,color="EM"),
                width = 0.1
  ) + geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_MImean,color="MI"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_MImean,color="MI"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_rhat_miss_MImean - rel_bias_rhat_miss_MIsd, ymax = rel_bias_rhat_miss_MImean + rel_bias_rhat_miss_MIsd,color="MI"),
                width = 0.1
  ) + geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_DAmean,color="DA"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_DAmean,color="DA"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_rhat_miss_DAmean - rel_bias_rhat_miss_DAsd, ymax = rel_bias_rhat_miss_DAmean + rel_bias_rhat_miss_DAsd,color="DA"),
                width = 0.1
  ) +geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_ccmean,color="cc"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_ccmean,color="cc"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_rhat_miss_ccmean - rel_bias_rhat_miss_ccsd, ymax = rel_bias_rhat_miss_ccmean + rel_bias_rhat_miss_ccsd,color="cc"),
                width = 0.1
  ) +geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_dropmean,color="drop"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_rhat_miss_dropmean,color="drop"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_rhat_miss_dropmean - rel_bias_rhat_miss_dropsd, ymax = rel_bias_rhat_miss_dropmean + rel_bias_rhat_miss_dropsd,color="drop"),
                width = 0.1
  ) +geom_line(data=res_sum,aes(x=n_adj_full,y=rel_bias_rhat_full_dropmean,color="full"))+
  geom_point(data=res_sum,aes(x=n_adj_full,y=rel_bias_rhat_full_dropmean,color="full"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_full,ymin = rel_bias_rhat_full_dropmean - rel_bias_rhat_full_dropsd, ymax = rel_bias_rhat_full_dropmean + rel_bias_rhat_full_dropsd,color="full"),
                width = 0.1
  ) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  labs(color="Legend text")+
  xlab("ln(Number of observations)") +
  ylab("(estimate - true)/true") +
  ggtitle(expression(r))



# create plots
a <- ggplot() +
  geom_line(data = res_sum, aes(x = n_adj_miss, y = rel_bias_ahat_miss_EMmean, color = "EM"))+
  geom_point(data = res_sum, aes(x = n_adj_miss, y = rel_bias_ahat_miss_EMmean, color = "EM"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_ahat_miss_EMmean - rel_bias_ahat_miss_EMsd, ymax = rel_bias_ahat_miss_EMmean + rel_bias_ahat_miss_EMsd,color="EM"),
                width = 0.1
  ) + geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_MImean,color="MI"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_MImean,color="MI"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_ahat_miss_MImean - rel_bias_ahat_miss_MIsd, ymax = rel_bias_ahat_miss_MImean + rel_bias_ahat_miss_MIsd,color="MI"),
                width = 0.1
  ) + geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_DAmean,color="DA"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_DAmean,color="DA"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_ahat_miss_DAmean - rel_bias_ahat_miss_DAsd, ymax = rel_bias_ahat_miss_DAmean + rel_bias_ahat_miss_DAsd,color="DA"),
                width = 0.1
  ) +geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_ccmean,color="cc"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_ccmean,color="cc"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_ahat_miss_ccmean - rel_bias_ahat_miss_ccsd, ymax = rel_bias_ahat_miss_ccmean + rel_bias_ahat_miss_ccsd,color="cc"),
                width = 0.1
  ) +geom_line(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_dropmean,color="drop"))+
  geom_point(data=res_sum,aes(x=n_adj_miss,y=rel_bias_ahat_miss_dropmean,color="drop"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_miss,ymin = rel_bias_ahat_miss_dropmean - rel_bias_ahat_miss_dropsd, ymax = rel_bias_ahat_miss_dropmean + rel_bias_ahat_miss_dropsd,color="drop"),
                width = 0.1
  ) +geom_line(data=res_sum,aes(x=n_adj_full,y=rel_bias_ahat_full_dropmean,color="full"))+
  geom_point(data=res_sum,aes(x=n_adj_full,y=rel_bias_ahat_full_dropmean,color="full"))+
  geom_errorbar(data = res_sum,
                aes(x=n_adj_full,ymin = rel_bias_ahat_full_dropmean - rel_bias_ahat_full_dropsd, ymax = rel_bias_ahat_full_dropmean + rel_bias_ahat_full_dropsd,color="full"),
                width = 0.1
  ) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  labs(color="Legend text")+
  xlab("ln(Number of observations)") +
  ylab("(estimate - true)/true") +
  ggtitle(expression(alpha))




ggsave(
  filename = here(paste0("Population sim and real/Figures/bias_checks_expanded",in_args[1],".png")),
  plot = r + a,
  device = "png",
  width = 8, height = 3,
  units = "in", dpi = 300
)




