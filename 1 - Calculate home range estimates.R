# Title: Biweekly home ranges
# Subtitle: 1 - Calculate biweekly home range estimates
# Author: Nathan D. Hooven
# Email: nathan.hooven@uky.edu
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 11 May 2022
# Date completed: 12 May 2022
# Date modified: 
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(sp)                # create spatial objects
library(raster)            # use DEM
library(amt)               # create track and resample
library(rgdal)             # write shapefiles
library(lubridate)         # deal with dates
library(adehabitatHR)      # create home ranges
library(adehabitatLT)      # create trajectories
library(stringr)           # string manipulation
library(move)              # create "telemetry" objects
library(ctmm)              # continuous time movement models and AKDEs

#_____________________________________________________________________________________________________________
# 2. Read in cleaned telemetry data ----
#_____________________________________________________________________________________________________________

# read in DEM raster for surface area calculation
DEM <- raster(paste0("G:/Elk project/Elk Zone rasters (7-20-21)", "/DEM_UTM.tif"))

DEM.crs <- crs(DEM)

# read in vectronic relocation data
vectronic.data <- read.csv("G:/Elk project/Data analysis/Raw data processing/Relocations_vectronic_1.csv")

# define UTM projection
projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")

# make sure t is a POSIXct object (even though the time is in Eastern, we'll let it be UTC for as.Date)
vectronic.data$t <- as.POSIXct(mdy_hm(vectronic.data$t))

# select only columns we need
vectronic.data <- vectronic.data %>% dplyr::select(x, y, t, burst, DOP, Animal, Year, Season)

# filter out relocations from after Jan 2021 for 2020 collars (retained originally for movement project)
vectronic.data.2020 <- vectronic.data %>% filter(Year == 2020) %>%
                                          filter(t < as.POSIXct("2021-02-01"))

vectronic.data.2021 <- vectronic.data %>% filter(Year == 2021)

# bind together
vectronic.data.all <- rbind(vectronic.data.2020, vectronic.data.2021)

# load error model
load("G:/Elk project/Data analysis/Raw data processing/vectronic_uere.RData")

#_____________________________________________________________________________________________________________
# 3. Define animals per year (removing confirmed brainers) ----
#_____________________________________________________________________________________________________________

collars.2020 <- c(37703, 37704, 37705, 37706, 37707, 37708, 37709, 37710, 37711, 37712,
                  37714, 37715, 37716, 37717, 37718, 37719, 37720, 37721, 37722,
                  37723, 37724, 37725, 37726, 37727)

collars.2021 <- c(45469, 45470, 45492, 45493, 45494, 45496, 45497, 45498, 45499,
                  45500, 45501, 45502, 45503, 45504, 45505, 45506, 45507, 45508, 
                  45510, 45511, 46391, 46392, 46393, 46394, 46396, 46397, 46399,
                  46400, 46401, 46402)

#_____________________________________________________________________________________________________________
# 4. Define analysis periods ----
#_____________________________________________________________________________________________________________
# 4a. 2020 ----
#_____________________________________________________________________________________________________________

periods.2020 <- data.frame("Period" = 1:26,
                           "Start.init" = as.Date("2020-02-01"))

periods.2020 <- periods.2020 %>% mutate(Start = Start.init + (Period - 1)*14,
                                        End = Start.init + (Period - 1)*14 + 14) %>%
                                 dplyr::select(-Start.init)

# change final date to include last days of Dec
periods.2020[26, 3] <- as.Date("2021-02-01")

#_____________________________________________________________________________________________________________
# 4b. 2021 ----
#_____________________________________________________________________________________________________________

periods.2021 <- data.frame("Period" = 1:26,
                           "Start.init" = as.Date("2021-02-01"))

periods.2021 <- periods.2021 %>% mutate(Start = Start.init + (Period - 1)*14,
                                        End = Start.init + (Period - 1)*14 + 14) %>%
                                 dplyr::select(-Start.init)

# change final date to include last days of Dec
periods.2021[26, 3] <- as.Date("2022-02-01")

#_____________________________________________________________________________________________________________
# 5. Calculate HR estimates ----
#_____________________________________________________________________________________________________________
# 5a. 2020 ----
#_____________________________________________________________________________________________________________

all.data.2020 <- data.frame()

for (x in collars.2020[11:24]) {
  
  # subset collar data
  CollarID <- x
  
  collar.data <- vectronic.data.all %>% filter(Animal == CollarID)
  
  # loop to subset all biweekly periods
  collar.hrs <- data.frame("Animal" = CollarID,
                           "Estimator" = c("AKDE", "AKDE", "BBMM", "BBMM"),
                           "Area" = c("plan", "surf", "plan", "surf"))
  
  for (y in periods.2020$Period) {
    
    # subset each period
    start.date <- periods.2020$Start[periods.2020$Period == y]
    end.date <- periods.2020$End[periods.2020$Period == y]
    
    focal.period <- collar.data %>% filter(t >= start.date & t < end.date)
    
    if (nrow(focal.period) >= 15) {
    
    # AKDE
    # create a track in amt
    focal.track <- make_track(focal.period, x, y, t, burst, Animal, DOP,
                              crs = 32617)
    
    # convert to telemetry object
    focal.telem <- as_telemetry(focal.track,
                                timeformat = "auto",
                                timezone = "UTC",
                                projection = 32617,
                                keep = TRUE)
    
    # fit error-informed movement models
    # assign error using the error model
    uere(focal.telem) <- all.animals.uere
    
    # guesstimated model parameters from the variogram
    guess.param <- variogram.fit(variogram(focal.telem), 
                                 name = "guess.param", 
                                 interactive = FALSE)
  
    # model selection
    fitted.mods <- ctmm.select(focal.telem, 
                               CTMM = guess.param, 
                               verbose = TRUE)
    
    # fit autocorrelated kernel density estimator and extract areas
    # fit an AKDE
    akde.1 <- akde(focal.telem, CTMM = fitted.mods[[1]])
    
    # extract 95% contour area
    akde.area.plan <- summary(akde.1, level.UD = 0.95)$CI[2]
    
    # extract 95% contour as a polygon
    akde.HR.all <- SpatialPolygonsDataFrame.UD(akde.1, level.UD = 0.95)
    
    akde.HR <- akde.HR.all[akde.HR.all$name == "1 95% est", ]
    
    # convert to UTM
    akde.HR <- spTransform(akde.HR, projection)
    
    # BBMM
    # create a trajectory
    focal.ltraj <- as.ltraj(xy = focal.period[,c("x", "y")],
                           date = focal.period$t,
                           id = 1,
                           burst = focal.period$burst,
                           proj4string = projection)
    
    # fit Brownian bridge movement model UD
    # Sig1 using tested SD of location error (16.57726 m)
    lik <- liker(focal.ltraj, sig2 = 16.57726, rangesig1 = c(0, 50))
  
    # estimate BBMM
    bbmm <- kernelbb(focal.ltraj, 
                     sig1 = lik$`1`$sig1, 
                     sig2 = 16.57726, 
                     grid = 200,
                     extent = 1)
    
    # vectorize 95% contour
    bbmm.HR <- getverticeshr(bbmm, percent = 95)
    
    # add CRS
    bbmm.HR@proj4string <- projection
    
    # Estimate area
    bbmm.area.plan <- kernel.area(bbmm,
                                  percent = c(95),
                                  unin = "m",
                                  unout = "km2")
    
    # calculate topographic home ranges (surface area)
    # crop DEM to each polygon
    DEM.crop.akde <- crop(DEM, akde.HR)
    DEM.crop.bbmm <- crop(DEM, bbmm.HR)
    
    # mask DEM to each polygon
    DEM.mask.akde <- mask(DEM.crop.akde, akde.HR)
    DEM.mask.bbmm <- mask(DEM.crop.bbmm, bbmm.HR)
    
    # convert to SpatialGridDataFrames
    DEM.sgdf.akde <- as(DEM.mask.akde, "SpatialGridDataFrame")
    DEM.sgdf.bbmm <- as(DEM.mask.bbmm, "SpatialGridDataFrame")
    
    # compute surface areas
    akde.area.surf <- surfaceArea(DEM.sgdf.akde) / 1000000
    bbmm.area.surf <- surfaceArea(DEM.sgdf.bbmm) / 1000000
      
    
    # assign areas to a df
    collar.areas <- data.frame("col1" = c(akde.area.plan,
                                          akde.area.surf,
                                          bbmm.area.plan,
                                          bbmm.area.surf))
    
    colnames(collar.areas) <- paste0("Period.", y)
    
    # bind to all individual home ranges
    collar.hrs <- cbind(collar.hrs, collar.areas)
    
    } else {
      
    collar.areas <- data.frame("col1" = c(NA,
                                          NA,
                                          NA,
                                          NA))  
    
    colnames(collar.areas) <- paste0("Period.", y)
      
    collar.hrs <- cbind(collar.hrs, collar.areas)
      
    }
    
  }
  
  # bind to master df
  all.data.2020 <- rbind(all.data.2020, collar.hrs)
  
}

#_____________________________________________________________________________________________________________
# 5b. 2021 ----
#_____________________________________________________________________________________________________________

all.data.2021 <- data.frame()

for (x in collars.2021) {
  
  # subset collar data
  CollarID <- x
  
  collar.data <- vectronic.data.all %>% filter(Animal == CollarID)
  
  # loop to subset all biweekly periods
  collar.hrs <- data.frame("Animal" = CollarID,
                           "Estimator" = c("AKDE", "AKDE", "BBMM", "BBMM"),
                           "Area" = c("plan", "surf", "plan", "surf"))
  
  for (y in periods.2021$Period) {
    
    # subset each period
    start.date <- periods.2021$Start[periods.2021$Period == y]
    end.date <- periods.2021$End[periods.2021$Period == y]
    
    focal.period <- collar.data %>% filter(t >= start.date & t < end.date)
    
    if (nrow(focal.period) >= 15) {
    
    # AKDE
    # create a track in amt
    focal.track <- make_track(focal.period, x, y, t, burst, Animal, DOP,
                              crs = 32617)
    
    # convert to telemetry object
    focal.telem <- as_telemetry(focal.track,
                                timeformat = "auto",
                                timezone = "UTC",
                                projection = 32617,
                                keep = TRUE)
    
    # fit error-informed movement models
    # assign error using the error model
    uere(focal.telem) <- all.animals.uere
    
    # guesstimated model parameters from the variogram
    guess.param <- variogram.fit(variogram(focal.telem), 
                                 name = "guess.param", 
                                 interactive = FALSE)
  
    # model selection
    fitted.mods <- ctmm.select(focal.telem, 
                               CTMM = guess.param, 
                               verbose = TRUE)
    
    # fit autocorrelated kernel density estimator and extract areas
    # fit an AKDE
    akde.1 <- akde(focal.telem, CTMM = fitted.mods[[1]])
    
    # extract 95% contour area
    akde.area.plan <- summary(akde.1, level.UD = 0.95)$CI[2]
    
    # extract 95% contour as a polygon
    akde.HR.all <- SpatialPolygonsDataFrame.UD(akde.1, level.UD = 0.95)
    
    akde.HR <- akde.HR.all[akde.HR.all$name == "1 95% est", ]
    
    # convert to UTM
    akde.HR <- spTransform(akde.HR, projection)
    
    # BBMM
    # create a trajectory
    focal.ltraj <- as.ltraj(xy = focal.period[,c("x", "y")],
                           date = focal.period$t,
                           id = 1,
                           burst = focal.period$burst,
                           proj4string = projection)
    
    # fit Brownian bridge movement model UD
    # Sig1 using tested SD of location error (16.57726 m)
    lik <- liker(focal.ltraj, sig2 = 16.57726, rangesig1 = c(0, 50))
  
    # estimate BBMM
    bbmm <- kernelbb(focal.ltraj, 
                     sig1 = lik$`1`$sig1, 
                     sig2 = 16.57726, 
                     grid = 200,
                     extent = 1)
    
    # vectorize 95% contour
    bbmm.HR <- getverticeshr(bbmm, percent = 95)
    
    # add CRS
    bbmm.HR@proj4string <- projection
    
    # Estimate area
    bbmm.area.plan <- kernel.area(bbmm,
                                  percent = c(95),
                                  unin = "m",
                                  unout = "km2")
    
    # calculate topographic home ranges (surface area)
    # crop DEM to each polygon
    DEM.crop.akde <- crop(DEM, akde.HR)
    DEM.crop.bbmm <- crop(DEM, bbmm.HR)
    
    # mask DEM to each polygon
    DEM.mask.akde <- mask(DEM.crop.akde, akde.HR)
    DEM.mask.bbmm <- mask(DEM.crop.bbmm, bbmm.HR)
    
    # convert to SpatialGridDataFrames
    DEM.sgdf.akde <- as(DEM.mask.akde, "SpatialGridDataFrame")
    DEM.sgdf.bbmm <- as(DEM.mask.bbmm, "SpatialGridDataFrame")
    
    # compute surface areas
    akde.area.surf <- surfaceArea(DEM.sgdf.akde) / 1000000
    bbmm.area.surf <- surfaceArea(DEM.sgdf.bbmm) / 1000000
      
    
    # assign areas to a df
    collar.areas <- data.frame("col1" = c(akde.area.plan,
                                          akde.area.surf,
                                          bbmm.area.plan,
                                          bbmm.area.surf))
    
    colnames(collar.areas) <- paste0("Period.", y)
    
    # bind to all individual home ranges
    collar.hrs <- cbind(collar.hrs, collar.areas)
    
    } else {
      
    collar.areas <- data.frame("col1" = c(NA,
                                          NA,
                                          NA,
                                          NA))  
    
    colnames(collar.areas) <- paste0("Period.", y)
      
    collar.hrs <- cbind(collar.hrs, collar.areas)
      
    }
    
  }
  
  # bind to master df
  all.data.2021 <- rbind(all.data.2021, collar.hrs)
  
}

#_____________________________________________________________________________________________________________
# 6. Write to .csv ----
#_____________________________________________________________________________________________________________

all.data.2020$Year <- 2020
all.data.2021$Year <- 2021

all.data <- rbind(all.data.2020, all.data.2021)

write.csv(all.data, "all_data.csv")
