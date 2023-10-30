

# Test

flist <- list.files("Functions/", full.names = T)
lapply(flist, source)

y <- readRDS(
  here::here("data/missingDatasets/pois_sim_randMiss_A.rds")
)[[1]]$y[[5]]

fit <- fit_ricker_DA(
  y,
  burnin = 5000,
  priors_list = list(
    m_r = 0,
    sd_r = 2.5,
    m_lalpha = -3,
    sd_lalpha = 1
  ),
  nthin = 5,
  return_y = T
)


df <- tibble(
  y = fit$estim[grep("y", names(fit$estim))],
  se = fit$se[grep("y", names(fit$se))],
  low = fit$lower[grep("y", names(fit$lower))],
  high = fit$upper[grep("y", names(fit$upper))],
  t = 1:length(y)
)

df$miss <- 0
df$miss[which(df$se != 0)] <- 1
df$miss <- as.factor(df$miss)

ggplot(data = df, aes(x = t, y = y)) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
  geom_point(aes(shape = miss)) +
  theme_classic() +
  geom_line(linewidth = 0.1) +
  scale_shape_manual(values = c(19, 1)) +
  ylab("Population size")

ggsave(
  here("figures/DA_example.png"),
  device = "png",
  height = 3,
  width = 6,
  units = "in"
)
  