# Evaluating Amelia!

# Load in a few ricker data sets
library(Amelia)

gauss_sim_randMiss <- readRDS("data/missingDatasets/gauss_sim_randMiss_A.rds")
randMissA=readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")
missParams=readRDS("data/missingDatasets/pois_sim_params.rds")


# Data set 1
y_full <- gauss_sim_randMiss[[5]]$y[[1]]
y=gauss_sim_randMiss[[5]]$y[[4]]

# we want to use tscsplot in order to make comparable plots, 
# this means no consecutive missing values since that would require multiple rounds of Amelia
# this is why we are adding in numbers for any consecutive missing values
which(is.na(y_full))
which(is.na(y))
dat <- data.frame(time = seq(1:100),
                   y = y[1:100],
                   light = gauss_sim_randMiss[[5]]$sim_params$X[1:100,2],
                   discharge = gauss_sim_randMiss[[5]]$sim_params$X[1:100,3],
                   cs_dummy = 1
)


# Data set 2
y_full <- gauss_sim_randMiss[[6]]$y[[1]]
y=gauss_sim_randMiss[[6]]$y[[4]]
which(is.na(y_full))
which(is.na(y))
y[3]=y_full[3]
y[6]=y_full[6]
y[27]=y_full[27]
y[35]=y_full[35]
y[41]=y_full[41]
y[43]=y_full[43]
y[78]=y_full[78]
dat <- data.frame(time = seq(1:100),
                   y = y[1:100],
                   light = gauss_sim_randMiss[[6]]$sim_params$X[1:100,2],
                   discharge = gauss_sim_randMiss[[6]]$sim_params$X[1:100,3],
                   cs_dummy = 1
)


# Data set 3
y_full <- gauss_sim_randMiss[[7]]$y[[1]]
y=gauss_sim_randMiss[[7]]$y[[4]]
which(is.na(y_full))
which(is.na(y))
y[34]=y_full[34]
y[44]=y_full[44]
y[82]=y_full[82]
dat <- data.frame(time = seq(1:100),
                   y = y[1:100],
                   light = gauss_sim_randMiss[[7]]$sim_params$X[1:100,2],
                   discharge = gauss_sim_randMiss[[7]]$sim_params$X[1:100,3],
                   cs_dummy = 1
           )

# do the 4 Amelia methods on each

# make data frame with pop at time t and time t-1
# add dummy cross section variable in order to use tscsPlot function
imputationsnum=5

# basic screen printout for keeping track
p2samelia=1 

par(mfrow=c(3,2))

# original (none)
amelia1sim<-amelia(dat, m=imputationsnum, ts="time", p2s=p2samelia, cs="cs_dummy")
tscsPlot(amelia1sim,2,1,main="none")
points(y_full[1:100], pch = 1)

# 1 lags and leads
amelia1sim<-amelia(dat, m=imputationsnum, ts="time", p2s=p2samelia, cs="cs_dummy",
                   lags="y", leads="y")
tscsPlot(amelia1sim,2,1,main="lags and leads")
points(y_full[1:100], pch = 1)

# 2 polytime=1
amelia1sim<-amelia(dat, m=imputationsnum, ts="time", p2s=p2samelia, cs="cs_dummy", polytime = 1)
tscsPlot(amelia1sim,2,1,main="polytime 1")
points(y_full[1:100], pch = 1)

# 3 polytime=2
amelia1sim<-amelia(dat, m=imputationsnum, ts="time", p2s=p2samelia, cs="cs_dummy", polytime = 2)
tscsPlot(amelia1sim,2,1,main="polytime 2")
points(y_full[1:100], pch = 1)

# 4 polytime=3
amelia1sim<-amelia(dat, m=imputationsnum, ts="time",  p2s=p2samelia, cs="cs_dummy", polytime = 3)
tscsPlot(amelia1sim,2,1,main="polytime 3")
points(y_full[1:100], pch = 1)



# try to add comparison to true values, to see if there are any improvements?
