################################################################
# This script double checks the Ricker data augmentation function 
# some simulated datasets
################################################################



# Test
set.seed(5526)
library(tidyverse)

flist <- list.files("Functions/", pattern = "\\.R$", full.names = T)
lapply(flist, source)


# ---- one-off Neg Binom testing ----

## ---- No missing ----

y <- ricker_sim(50, 0.8, 0.01, N0 = 50, err_fam = "neg_binom", psi = 5)

( test_1 <- fit_ricker_DA(y, fam = "neg_binom") )


## ---- With missing obs ----

y2 <- y
y2[sample(1:length(y), 6)] <- NA

( test2 <- fit_ricker_DA(y2, fam = "neg_binom") )


## ---- Creating offset dataframe ahead of time ----

# remove starting NAs
if(is.na(y2[1])){
  start <- min(which(!is.na(y2)))
  y2 <- y2[start:length(y2)]
}

y3 <- data.frame(
  yt = y2[2:length(y2)],
  ytm1 = y2[1:(length(y2) - 1)]
)

(test_3 <- fit_ricker_DA(y3, fam = "neg_binom", off_patch = TRUE))


# ---- one-off Poisson testing ----

## ---- No missing ----

y <- ricker_sim(50, 0.8, 0.01, N0 = 50, err_fam = "poisson")

( test_4 <- fit_ricker_DA(y, fam = "poisson") )


## ---- With missing obs ----

y2 <- y
y2[sample(1:length(y), 6)] <- NA

( test_5 <- fit_ricker_DA(y2, fam = "poisson") )


## ---- Creating offset dataframe ahead of time ----

# remove starting NAs
if(is.na(y2[1])){
  start <- min(which(!is.na(y2)))
  y2 <- y2[start:length(y2)]
}

y3 <- data.frame(
  yt = y2[2:length(y2)],
  ytm1 = y2[1:(length(y2) - 1)]
)

( test_6 <- fit_ricker_DA(y3, fam = "poisson", off_patch = TRUE) )

  