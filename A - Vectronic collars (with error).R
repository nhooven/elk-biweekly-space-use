# Title: Raw data processing
# Subtitle: 1 - Vectronic collars (with error estimates)
# Author: Nathan D. Hooven
# Email: nathan.hooven@uky.edu
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 16 Feb 2021
# Date completed: 16 Feb 2021
# Date modified: 8 Feb 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(sp)                # create spatial objects
library(amt)               # create track and resample
library(rgdal)             # write shapefiles
library(lubridate)         # deal with dates
library(adehabitatHR)      # create home ranges
library(adehabitatLT)      # create trajectories
library(stringr)           # string manipulation
library(ctmm)              # apply error models

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

# all elk names to include (I'll include all of them so we can model all home ranges and assess autocorrelation)
# 2020 collars
elk.names.2020 <- as.vector(c("37703", "37704", "37705", "37706", "37707", "37708", "37709", "37710", "37711",
                              "37712", "37713", "37714", "37715", "37716", "37717", "37718", "37719", "37720",
                              "37721", "37722", "37723", "37724", "37725", "37726", "37727"))

# 2021 collars
elk.names.2021 <- as.vector(c("45469", "45470", "45492", "45493", "45494", "45495", "45496", "45497", "45498",
                              "45499", "45500", "45501", "45502", "45503", "45504", "45505", "45506", "45507",
                              "45508", "45509", "45510", "45511", "46391", "46392", "46393", "46394", "46395", 
                              "46396", "46397", "46399", "46400", "46401", "46402"))

# folder with csvs
exp.dir <- "G:/Elk project/Vectronic data/" 

# load error model
load("vectronic_uere.RData")

#_____________________________________________________________________________________________________________
# 3. 2020 - For loop to run through all elk ----
#_____________________________________________________________________________________________________________

# Create blank data frames to hold all data
elk.relocations.2020 <- data.frame()

for (x in elk.names.2020){
  
  #_____________________________________________________________________________________________________________
  # 3a. Select data ----
  #_____________________________________________________________________________________________________________
  
  CollarID <- x
  
  Animal  <- read.csv(paste0(exp.dir, "Collar", CollarID, "_GPS_Default Storage.csv"))
  
  # Timeframe
  Start.Date <- as.POSIXct("2/2/2020", tz = "America/New_York", tryFormats = "%m/%d/%Y")
  End.Date <- as.POSIXct("8/1/2021", tz = "America/New_York", tryFormats = "%m/%d/%Y")
  
  #_____________________________________________________________________________________________________________
  # 3b. Filter relocations by date ----
  #_____________________________________________________________________________________________________________
  
  # coerce local date to Date format
  Animal$LMT_Date <- as.Date(Animal$LMT_Date, tryFormats = "%m/%d/%Y")
  
  # filter dates of interest
  Animal.1 <- Animal %>% filter(LMT_Date >= Start.Date & LMT_Date <= End.Date)
  
  #_____________________________________________________________________________________________________________
  # 3c. Resample to ~ 13 hours and filter "bad fixes" using collar error model ----
  #_____________________________________________________________________________________________________________
  
  # filter "no fix" relocations out
  Animal.1 <- Animal.1[Animal.1$FixType == "val. GPS-3D",]
  
  # define projection
  projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
  
  # delete "17" zone from beginning of Easting string
  Animal.1$Easting <- str_sub(Animal.1$Easting, 3, 11)
  Animal.1$Easting <- as.numeric(Animal.1$Easting)
  
  # Add date and time variables together convert time variable to 24-hour time
  # Create timestamp variable
  Animal.1$datetime <- paste(Animal.1$LMT_Date, Animal.1$LMT_Time)
  Animal.1$Timestamp <- as.POSIXct(ymd_hms(Animal.1$datetime), tz = "America/New_York")
  
  # create track in amt
  Animal.track <- make_track(Animal.1, 
                             .x = Easting,
                             .y = Northing,
                             .t = Timestamp,
                             crs = 32617,
                             DOP = DOP)
  
  # resample
  Animal.track.1 <- track_resample(Animal.track,
                                   rate = hours(13),
                                   tolerance = hours(2),
                                   start = 1)
  
  # create telemetry object
  Animal.telem <- as_telemetry(Animal.track.1,
                               timeformat = "auto",
                               timezone = "America/New_York",
                               projection = 32617,
                               keep = TRUE)
  
  # assign error using calibration-data error model
  uere(Animal.telem) <- all.animals.uere
  
  # remove outliers with speed > 0.5 m/s and distance > 10 km
  Animal.telem.out <- outlie(Animal.telem)
  
  out.rules <- Animal.telem.out$distance > 8000 | Animal.telem.out$speed > 0.08
  
  # remove outliers from track
  Animal.track.2 <- Animal.track.1[!out.rules,]
  
  #_____________________________________________________________________________________________________________
  # 3d. Create a trajectory and bind to master df ----
  #_____________________________________________________________________________________________________________
  
  Animal.ltraj <- as.ltraj(xy = Animal.track.2[,c("x_", "y_")],
                           date = Animal.track.2$t_,
                           id = 1,
                           proj4string = projection)
  
  # Name columns "x" and "y" in trajectory
  Animal.ltraj[[1]]$x <- Animal.ltraj[[1]]$x_
  Animal.ltraj[[1]]$y <- Animal.ltraj[[1]]$y_
  
  # Create a SpatialPoints object too
  Animal.sp <- SpatialPoints(coords = Animal.track.2[,c("x_", "y_")],
                             proj4string = projection)
  
  # Create data frame and bind to master relocation data frame
  Animal.relocations <- data.frame(x = Animal.track.2$x_,
                                   y = Animal.track.2$y_,
                                   t = Animal.track.2$t_,
                                   burst = Animal.track.2$burst_,
                                   DOP = Animal.track.2$DOP,
                                   Animal = CollarID,
                                   Year = "2020")
  
  # write to master data frame
  elk.relocations.2020 <- rbind(elk.relocations.2020, Animal.relocations)
  
}

#_____________________________________________________________________________________________________________
# 4. 2021 - For loop to run through all elk ----
#_____________________________________________________________________________________________________________

# Create blank data frames to hold all data
elk.relocations.2021 <- data.frame()

for (x in elk.names.2021){
  
  #_____________________________________________________________________________________________________________
  # 3a. Select data ----
  #_____________________________________________________________________________________________________________
  
  CollarID <- x
  
  Animal  <- read.csv(paste0(exp.dir, "Collar", CollarID, "_GPS_Default Storage.csv"))
  
  # Timeframe
  Start.Date <- as.POSIXct("2/2/2021", tz = "America/New_York", tryFormats = "%m/%d/%Y")
  End.Date <- as.POSIXct("1/31/2022", tz = "America/New_York", tryFormats = "%m/%d/%Y")
  
  #_____________________________________________________________________________________________________________
  # 3b. Filter relocations by date ----
  #_____________________________________________________________________________________________________________
  
  # coerce local date to Date format
  Animal$LMT_Date <- as.Date(Animal$LMT_Date, tryFormats = "%m/%d/%Y")
  
  # filter dates of interest
  Animal.1 <- Animal %>% filter(LMT_Date >= Start.Date & LMT_Date <= End.Date)
  
  # filter "no fix" relocations out
  Animal.1 <- Animal.1[Animal.1$FixType == "val. GPS-3D",]
  Animal.1 <- Animal.1[Animal.1$DOP < 10,]
  
  #_____________________________________________________________________________________________________________
  # 3c. Resample to ~ 13 hours and filter "bad fixes" using collar error model ----
  #_____________________________________________________________________________________________________________
  
  # filter "no fix" relocations out
  Animal.1 <- Animal.1[Animal.1$FixType == "val. GPS-3D",]
  
  # define projection
  projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
  
  # delete "17" zone from beginning of Easting string
  Animal.1$Easting <- str_sub(Animal.1$Easting, 3, 11)
  Animal.1$Easting <- as.numeric(Animal.1$Easting)
  
  # Add date and time variables together convert time variable to 24-hour time
  # Create timestamp variable
  Animal.1$datetime <- paste(Animal.1$LMT_Date, Animal.1$LMT_Time)
  Animal.1$Timestamp <- as.POSIXct(ymd_hms(Animal.1$datetime), tz = "America/New_York")
  
  # create track in amt
  Animal.track <- make_track(Animal.1, 
                             .x = Easting,
                             .y = Northing,
                             .t = Timestamp,
                             crs = 32617,
                             DOP = DOP)
  
  # resample
  Animal.track.1 <- track_resample(Animal.track,
                                   rate = hours(13),
                                   tolerance = hours(2),
                                   start = 1)
  
  # create telemetry object
  Animal.telem <- as_telemetry(Animal.track.1,
                               timeformat = "auto",
                               timezone = "America/New_York",
                               projection = 32617,
                               keep = TRUE)
  
  # assign error using calibration-data error model
  uere(Animal.telem) <- all.animals.uere
  
  # remove outliers with speed > 0.8 m/s and distance > 8 km
  Animal.telem.out <- outlie(Animal.telem)
  
  out.rules <- Animal.telem.out$distance > 8000 | Animal.telem.out$speed > 0.08
  
  # remove outliers from track
  Animal.track.2 <- Animal.track.1[!out.rules,]
  
  #_____________________________________________________________________________________________________________
  # 3d. Create a trajectory and bind to master df ----
  #_____________________________________________________________________________________________________________
  
  Animal.ltraj <- as.ltraj(xy = Animal.track.2[,c("x_", "y_")],
                           date = Animal.track.2$t_,
                           id = 1,
                           proj4string = projection)
  
  # Name columns "x" and "y" in trajectory
  Animal.ltraj[[1]]$x <- Animal.ltraj[[1]]$x_
  Animal.ltraj[[1]]$y <- Animal.ltraj[[1]]$y_
  
  # Create a SpatialPoints object too
  Animal.sp <- SpatialPoints(coords = Animal.track.2[,c("x_", "y_")],
                             proj4string = projection)
  
  # Create data frame and bind to master relocation data frame
  Animal.relocations <- data.frame(x = Animal.track.2$x_,
                                   y = Animal.track.2$y_,
                                   t = Animal.track.2$t_,
                                   burst = Animal.track.2$burst_,
                                   DOP = Animal.track.2$DOP,
                                   Animal = CollarID,
                                   Year = "2021")
  
  # write to master data frame
  elk.relocations.2021 <- rbind(elk.relocations.2021, Animal.relocations)
  
}

#_____________________________________________________________________________________________________________
# 5. Bind together and create season variables ----
#_____________________________________________________________________________________________________________

# bind together
elk.relocations <- rbind(elk.relocations.2020, elk.relocations.2021)

# create season variables (in a smart way)
elk.relocations <- elk.relocations %>% mutate(Season = dplyr::case_when(substr(elk.relocations$t, 6, 7) %in%  c("02", "03", "04") ~ "WS",  
                                                                        substr(elk.relocations$t, 6, 7) %in%  c("05", "06", "07") ~ "SU", 
                                                                        substr(elk.relocations$t, 6, 7) %in%  c("08", "09", "10") ~ "UA",
                                                                        substr(elk.relocations$t, 6, 7) %in%  c("11", "12", "01") ~ "AW"))

#_____________________________________________________________________________________________________________
# 6. Bind and write to final csv ----
#_____________________________________________________________________________________________________________

write.csv(elk.relocations, "Relocations_vectronic.csv")
