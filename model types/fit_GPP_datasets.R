# Run stan model with fixed observation error on actual GPP data
# Author: AM Carter

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

dat <- read_csv('data/NWIS_MissingTS_subset.csv')
mdat <- read_csv('data/NWIS_MissingTSinfo_subset.csv')

plot_gpp <- function(df, pp_fit = FALSE){
    ratio_QL <- max(df$light)/max(df$Q)
    GPP_plot <- ggplot(df, aes(date, GPP))+
      #geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="chartreuse4")+
      geom_point(size = 2, color="chartreuse4") + 
      geom_line(size = 1, color="chartreuse4")+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$site_name[1])+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12))
    if(pp_fit){
        GPP_plot <- df %>%
          rename('posterior estimate' = GPP_rep)%>%
          pivot_longer(cols = c('GPP', 'posterior estimate'), 
                       names_to = 'GPP', values_to = 'value') %>%
          ggplot(aes(date, value, col = GPP))+
          geom_errorbar(aes(ymin = GPP_rep.lower, ymax = GPP_rep.upper), 
                        width=0.2,color="grey")+
          geom_point(size=2) + geom_line(size=1)+
          scale_color_manual(values = c('chartreuse4', 'grey35'))+
          labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$site_name[1])+
          theme(legend.position = 'top', legend.title = element_blank(),
                panel.background = element_rect(color = "black", fill=NA, size=1),
                axis.title.x = element_blank(), axis.text.x = element_blank(),
                axis.text.y = element_text(size=12),
                axis.title.y = element_text(size=12))
    }
    data_plot <- ggplot(df, aes(date, Q*ratio_QL))+
      geom_point(data=df, aes(date, light), size=1.5, color="darkgoldenrod3")+
      geom_line(size=1, color="deepskyblue4")+
      scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
      labs(y=expression('Daily PPFD'))+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text = element_text(size=12),
            axis.title.y.left = element_text(size=12, color="darkgoldenrod3"),
            axis.title.y.right = element_text(size=12, color="deepskyblue4"),
            axis.text.x = element_text(angle=25, hjust = 1),
            strip.background = element_rect(fill="white", color="black"),
            strip.text = element_text(size=15))
    
    pp <- ggpubr::ggarrange(GPP_plot, data_plot, ncol = 1, align = 'v')
    
    return(pp)
    
}

# start with pine river

id <- mdat$site_name[4]
pr <- filter(dat, site_name == id)

# calculate observation error based on 95% CIs

sigma_obs <- (pr$GPP.upper - pr$GPP.lower)/3.92
pr <- pr %>%
  select(date, GPP, light, Q) %>%
  mutate(across(-date, ~zoo::na.approx(.)))%>%
  mutate(light.rel = light/max(light))


pr %>% mutate(site_name = 'Pine River') %>% plot_gpp()

#Prep models ####
# model file: "model types/fixed_oi_light_centered.stan"
model_l <- stan_model("model types/Stan/AR1_light_centered.stan")
model_lq <- stan_model("model types/Stan/AR1_light_Q_centered.stan")

#Create data object
data <- list(N = nrow(pr), P_obs = pr$GPP,
             sdo = sigma_obs, light = pr$light.rel)

#Run Stan
fit_l <- rstan::sampling(object=model_l, data = data,  
                       iter = 4000, chains = 4)

#Create data object
data <- list(N = nrow(pr), P_obs = pr$GPP,
             sdo = sigma_obs, light = pr$light.rel, Q = pr$Q,
             miss_vec = rep(1, nrow(pr)))

#Run Stan
fit_lq <- rstan::sampling(object=model_lq, data = data,  
                       iter = 4000, chains = 4)
# examine model outputs
traceplot(fit_l, pars=c("phi", "sdp", "beta"))
pairs(fit_l, pars=c("phi", "sdp","beta","lp__"))
stan_dens(fit_l, pars=c("phi", "sdp", "beta"))

traceplot(fit_lq, pars=c("phi", "sdp", "beta"))
pairs(fit_lq, pars=c("phi", "sdp","beta","lp__"))
plot(fit_lq, pars=c("phi", "sdp", "beta"))


# look at posterior predictions:
get_pp <- function(fit){
    fit_extract <- rstan::extract(fit)
    pp <- t(apply(fit_extract$GPP_rep, 2, 
            function(x) quantile(x, probs = c(0.025, 0.5, 0.975), names = F))) %>%
      data.frame() 
    names(pp) <- c('GPP_rep.lower', 'GPP_rep', 'GPP_rep.upper')
    
    return(pp)
    
}

pr_l <- bind_cols(pr, get_pp(fit_l))
plot_gpp(pr_l, pp_fit = TRUE)

pr_q <- bind_cols(pr, get_pp(fit_lq))
plot_gpp(pr_q, pp_fit = TRUE)
