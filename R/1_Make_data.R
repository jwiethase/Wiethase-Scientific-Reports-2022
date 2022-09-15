library(raster)
library(rgdal)
library(tidyverse)

fileprefix = "BG200m"
deg_name_pattern = '200m'

# Load the data  ---------------------------
deg <- do.call(stack, lapply(list.files(path = "data", pattern = '_200m', full.names = T), raster))

CAs <- readOGR("data/CAs.shp")
CAs <- spTransform(CAs, proj4string(deg))
CAs$Type <- as.factor(CAs$Type)

studyArea <- readOGR("data/NTRI_outline.shp")
studyArea <- spTransform(studyArea, proj4string(deg))

# Prepare data ---------------------------------
CA_rast <- rasterize(CAs, deg, "Type")   ### 1 = CCRO, 2 = NP / NCAA, 3 = WMA
CA_rast[is.na(CA_rast) & !is.na(deg[[1]])] <- 0  ### 0 = no protection, 1 = CCRO, 2 = NP / NCAA, 3 = WMA

# Most and least degraded at end of time, median of the 3 most recent years:
last_3_raster <- subset(deg, (dim(deg)[3]-2):dim(deg)[3])
last_3_median <- calc(last_3_raster, median)

deg.quants <- quantile(values(last_3_median), p = c(0, 0.25, .75, 1), na.rm = TRUE)   ### higher values are most degraded
deg_q <- last_3_median
deg_q[!is.na(deg_q)] <- 0
deg_q[last_3_median < deg.quants[2]] <- -1
deg_q[last_3_median > deg.quants[3]] <- 1  ## -1 is least degraded areas, 0 is middle 50%, 1 is upper 25%

# By rainfall - actual and quantiles - chirps
overall_annual_median <- raster("/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/GEE_files/GEE_export_images/NTRI_annual_median_rain_99-20.tif")  ## from here: https://code.earthengine.google.com/11749041205711977212e1483321579c
rain <- projectRaster(overall_annual_median, deg)
rain.quants <- quantile(values(rain), p = c(0, 0.25, .75, 1), na.rm = TRUE)
rain_q <- rain
rain_q[!is.na(rain_q)] <- 0
rain_q[rain < rain.quants[2]] <- -1
rain_q[rain > rain.quants[3]] <- 1  ## -1 is driest areas, 0 is middle 50%, 1 is wettest 25%
# writeRaster(rain_q, paste0('data/', fileprefix, '_rain_q.tif'), overwrite = TRUE)

 # Annual rainfall
annRain <- read.csv("/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/GEE_files/GEE_export_csv/NTRI_PrecipitationYearly_99-20.csv")  ## from here: https://code.earthengine.google.com/11749041205711977212e1483321579c
annRain$year <- as.factor(1999:2020)
annRain$deg <- NA

# By tlu - actual and quantiles
tlu <- raster('data/NTRI_TLU_2010.tif')     # From https://code.earthengine.google.com/3cca3e20e7be731704c7fc60846a47c9
tlu <- projectRaster(tlu, deg)
tlu.quants <- quantile(values(tlu), p = c(0, 0.25, .75, 1), na.rm = TRUE)
tlu_q <- tlu
tlu_q[!is.na(tlu_q)] <- 0
tlu_q[tlu < tlu.quants[2]] <- -1
tlu_q[tlu > tlu.quants[3]] <- 1  ## -1 is lowest cattle areas, 0 is middle 50%, 1 is highest cattle area 25%

# By human population - actual and quantiles  (landscan data)
pop <- raster("data/population.tiff")
pop <- projectRaster(pop, deg)
pop.quants <- quantile(values(pop), p = c(0, 0.25, .75, 1), na.rm = TRUE)
pop_q <- pop
pop_q[!is.na(pop_q)] <- 0
pop_q[pop < pop.quants[2]] <- -1
pop_q[pop > pop.quants[3]] <- 1  ## -1 is lowest people areas, 0 is middle 50%, 1 is highest people area 25%

# Degradation change year on year
deg_change <- subset(deg, 2:nlayers(deg)) - subset(deg, 1:(nlayers(deg)-1))
deg_change_stack <- stack(deg_change)

# Stack all layers together
raw_stack <- stack(deg, CA_rast, pop, pop_q, rain, rain_q, tlu, tlu_q, deg_q, deg_change)
raw_spdf <- as(raw_stack, "SpatialPointsDataFrame")
raw_spdf <- raw_spdf[!is.na(rowSums(raw_spdf@data)), ]
names(raw_spdf) <- c(paste0("DegInd", 1:21), 
                     "CA", "pop", "pop_q", "rain", "rain_q", "tlu", "tlu_q", "deg_q", paste0("DegChange", 1:20))

save(raw_spdf, file =  paste0('data/', fileprefix, '_raw_spdf.RData'))




