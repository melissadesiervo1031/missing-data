
# test

library(tidyverse)

flist <- list.files("Functions/", pattern = "\\.R$", full.names = T)
lapply(flist, source)

# ---- one-off Poisson testing ----

## ---- No missing ----

y <- ricker_sim(50, 0.8, 0.01, N0 = 20, err_fam = "poisson")

fit <- fit_ricker_cc(y)
fit_d <- fit_ricker_drop(y)
fit_em <- fit_ricker_EM(y)
fit_da <- fit_ricker_DA(y)

## ---- Some missing ----

y2 <- y
y2[sample(1:length(y), size = 6)] <- NA

fit2 <- fit_ricker_cc(y2)
fit2_d <- fit_ricker_drop(y2)
fit2_em <- fit_ricker_EM(y2)
fit2_da <- fit_ricker_DA(y2)


## ---- multi-series ----

y3 <- lapply(1:5, \(i) { ricker_sim(50, 0.8, 0.01, N0 = rpois(1, 10), err_fam = "poisson") })
for(i in 1:length(y3)){
  y3[[i]][sample(1:length(y3[[i]]), size = rpois(1, 5))] <- NA
}

y_df <- lapply(1:length(y3), function(i, y){
  data.frame(
    yt = y[[i]][2:length(y[[i]])],
    ytm1 = y[[i]][1:(length(y[[i]]) - 1)],
    patch = as.character(i)
  )
}, y = y3)

y_df <- Reduce(rbind, y_df)

fit3 <- fit_ricker_cc(y_df, off_patch = T)
fit3_d <- fit_ricker_drop(y_df, off_patch = T, patch_col = "patch")
fit3_em <- fit_ricker_EM(y_df, off_patch = T)


# ---- one-off Negative Binomial testing ----

## ---- No missing ----

y_nb <- ricker_sim(50, 0.8, 0.01, N0 = 20, err_fam = "neg_binom", psi = 5)

fit_nb <- fit_ricker_cc(y_nb, fam = "neg_binom")
fit_nb_d <- fit_ricker_drop(y_nb, fam = "neg_binom")

## ---- Some missing ----

y2_nb <- y_nb
y2_nb[sample(1:length(y_nb), size = 6)] <- NA

fit2_nb <- fit_ricker_cc(y2_nb, fam = "neg_binom")
fit2_nb_d <- fit_ricker_drop(y2_nb, fam = "neg_binom")


## ---- multi-series ----

y3_nb <- lapply(1:5, \(i) { ricker_sim(50, 0.8, 0.01, N0 = rpois(1, 10), err_fam = "neg_binom", psi = 5) })
for(i in 1:length(y3_nb)){
  y3_nb[[i]][sample(1:length(y3_nb[[i]]), size = rpois(1, 5))] <- NA
}

y_df_nb <- lapply(1:length(y3_nb), function(i, y){
  data.frame(
    yt = y[[i]][2:length(y[[i]])],
    ytm1 = y[[i]][1:(length(y[[i]]) - 1)],
    patch = as.character(i)
  )
}, y = y3_nb)

y_df_nb <- Reduce(rbind, y_df_nb)

fit3_nb <- fit_ricker_cc(y_df_nb, off_patch = T, fam = "neg_binom")
fit3_nb_d <- fit_ricker_drop(y_df_nb, off_patch = T, patch_col = "patch", fam = "neg_binom")









