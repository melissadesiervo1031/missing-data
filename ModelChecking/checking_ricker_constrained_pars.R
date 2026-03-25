
# test

library(tidyverse)

flist <- list.files("Functions/", pattern = "\\.R$", full.names = T)
lapply(flist, source)

# ---- one-off Poisson testing ----

## ---- No missing ----

y <- ricker_sim(50, 0.8, 0.01, N0 = 20, err_fam = "poisson")

init <- c(r=0.5, lalpha=log(0.02))

fit <- optim(init, ricker_count_neg_ll_cnstr, y = y, fam = "poisson", hessian = T)
