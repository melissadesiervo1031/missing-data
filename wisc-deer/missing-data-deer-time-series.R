# Workspace set up -------------------------------------------------------------
library(raster)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(gridExtra)

# Read in data -----------------------------------------------------------------

#dat <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/DeerAbundance.csv")
dat <- read.csv("DeerAbundance.csv", header=T, check.names=FALSE)
# Reference layers -------------------------------------------------------------

landuse <- raster::raster("./wlciw930.tif")
raster::crs(landuse) <- sp::CRS("+init=EPSG:6610")

wisc <- raster::getData(country="USA", level=2)
wisc <- wisc[wisc$NAME_1=="Wisconsin",]
wisc <- sp::spTransform(wisc, raster::crs(landuse))

raster::plot(landuse)
raster::plot(wisc, add=T)
raster::text(wisc, wisc$NAME_2, cex=1, font=2)

# Not worth processing the land cover data to get the counties
# The distinction between the northern forested counties and the rest of the
# state is pretty clear without reclassifying all the land cover types

# Northern counties

no_cty <- c("Douglas", "Bayfield", "Ashland", "Iron", "Vilas", "Forest", 
            "Florence", "Burnett", "Washburn", "Sawyer", "Rusk", "Price", 
            "Taylor", "Oneida", "Marinette", "Lincoln", "Langlade")

# Menominee should be included, but has no deer data, possibly because it
# is tribal land and managed differently

# Clean data, covariates -------------------------------------------------------

deerdat <- dat[dat$County%in%no_cty,-1]
#snow <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/SnowDepth.csv")
#mei <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/WinterMEI.csv")[1,-1]
pdo <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/WinterPDO.csv")[1,-1]
pdo<-(t(unname(pdo)))
pdo<-data.frame(pdo)
snowdepth <- read.csv("SnowDepth.csv", header=T, check.names=FALSE) ## average snow depth from December- March. Data were the average ofall available weather stations for each of 71 counties (rows) in Wisconsin across 36 years(columns, 1981-2016). Some missing values (NA), resulting in some counties beingdropped from final analyses. Data from https://www.ncdc.noaa.gov/.
mei <- read.csv("WinterMEI.csv", header=T, check.names=FALSE) ##Average Multivariate El Niño/Southern Oscillation Index (MEI) valuefrom December-March across 36 years (columns, 1981-2016). Values of the matrix arerepeated across rows to match the spatial dimensions of the DeerAbundance data matrix.Data from https://www.esrl.noaa.gov/psd/enso/mei
snow <- snowdepth[snowdepth$County%in%no_cty,-1] # only covar that differs by county



deer_long_county <- gather(deerdat, year, abundance, '1981':'2016', factor_key=TRUE)
deer_long<-deer_long_county %>% group_by(year) %>% summarise(totalpop = sum(abundance))


snow_long_county <- gather(snow, year, snowdepth, '1981':'2016', factor_key=TRUE)
snow_long<-snow_long_county %>% group_by(year) %>% summarise(avgsnow = mean(snowdepth, na.rm=TRUE))


# Model specification ----------------------------------------------------------

sink("deer_model.stan")

cat("
  
  /*----------------------- Data --------------------------*/
    data {
      int N; // Number of data
      int K; //Number of predictors
      vector[N] y; // Deer population observations
      vector[N] snow; // snow data
      vector[N] pdo; // PDO data
      real y0; // Initial data prior
    }
  
  /*----------------------- Parameters --------------------------*/
    
    parameters {
      real<lower=0> sdp; // Standard deviation of the process equation
      vector[K] beta;
      real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    }
  
  
  /*----------------------- Model --------------------------*/
    model {
      // Prior distributions
      sdp ~ normal(0, 1);
      phi ~ beta(1,1);
      beta[1] ~ normal(0,5);
      beta[2] ~ normal(0,5);
      beta[2] ~ normal(0,5);
      
      // Distribution for the first data
      y[1] ~ normal(y0, sdp);
      
      // Distributions for all other data
      for(t in 2:N){
        y[t] ~ normal(beta[1]+phi*y[t-1]+ beta[2]*snow+ beta[3]*pdo, sdp);// process model with error
      }
    }
   

  "
,fill=TRUE)
sink()

#Prep model
model<-"deer_model.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(deer_long$totalpop),y=log(deer_long$totalpop), y0=log(deer_long$totalpop[1]),
             snow=snow_long$avgsnow, pdo=pdo$pdo, K=3)

#Run Stan
fit<-rstan::sampling(object=model, data = data,  iter = 2000, chains = 4)

print(fit)
traceplot(fit, pars=c("sdp", "phi", "beta"))
fit_extract<-rstan::extract(fit)


##Plot density plots
plot_sdp <- stan_dens(fit, pars="sdp")+xlab("Process Error")
plot_phi <- stan_dens(fit, pars="phi")+xlab("Auto-correlation")
plot_b1 <- stan_dens(fit, pars="beta[1]")+xlab("Intercept")
plot_b2 <- stan_dens(fit, pars="beta[2]")+xlab("Snow-depth")
plot_b3 <- stan_dens(fit, pars="beta[3]")+xlab("PDO")

grid.arrange(plot_sdp,plot_phi, plot_b1,plot_b2,plot_b3,nrow=2)


y<-summary(fit, pars = c("y_rep"), probs = c(0.1,0.5,0.9))$summary
est<-exp(y[,1])
plot(deer_long$totalpop~deer_long$year)
points(est)


# JAGS is declarative, so autoregressive models aren't really possible 
# because including previous timesteps redefines the node


