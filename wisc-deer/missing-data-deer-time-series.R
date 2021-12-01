# Workspace set up -------------------------------------------------------------

# Read in data -----------------------------------------------------------------

dat <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/DeerAbundance.csv")

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
snow <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/SnowDepth.csv")
mei <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/WinterMEI.csv")[1,-1]
pdo <- read.csv("https://raw.githubusercontent.com/melissadesiervo1031/timeseries/main/WinterPDO.csv")[1,-1]

snow <- snow[snow$County%in%no_cty,-1] # only covar that differs by county

# Model specification ----------------------------------------------------------

# placeholder for Matt to add Stan autoregressive model

# JAGS is declarative, so autoregressive models aren't really possible 
# because including previous timesteps redefines the node


