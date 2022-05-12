# Title: Raw data processing
# Subtitle: 8a - Vectronic collar error model
# Author: Nathan D. Hooven
# Email: nathan.hooven@uky.edu
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 29 Dec 2021
# Date completed: 29 Dec 2021
# Date modified: 
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(sp)                # create spatial objects
library(lubridate)         # deal with dates
library(adehabitatHR)      # create home ranges
library(adehabitatLT)      # create trajectories
library(stringr)           # string manipulation
library(amt)               # make a track
library(ctmm)              # error modeling

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

# collars to include
elk.names.2020 <- as.vector(c(37703, 37704, 37705, 37706, 37707, 37708, 37709, 37710, 37711,
                              37712, 37713, 37714, 37715, 37716, 37717, 37718, 37719, 37720,
                              37721, 37722, 37723, 37724, 37725, 37726, 37727))

# Vectronic data folder
exp.dir <- "G:/Elk project/Vectronic data/"

start.date <- as.Date("11/18/2019", tryFormats = "%m/%d/%Y", tz = "America/New_York")
end.date <- as.Date("11/22/2019", tryFormats = "%m/%d/%Y", tz = "America/New_York")

#_____________________________________________________________________________________________________________
# 3. Fit error models per individual ----
#_____________________________________________________________________________________________________________

all.animals.list <- list()

all.animals <- data.frame()

for (x in elk.names.2020){
  
  CollarID <- x
  
  Animal  <- read.csv(paste0(exp.dir, "Collar", CollarID, "_GPS_Default Storage.csv"))
  
  Animal <- Animal %>% dplyr::select(No, CollarID, LMT_Date, LMT_Time, 
                                     Latitude...., Longitude...., 
                                     Height..m., DOP, FixType)
  
  # filter locations by date
  # coerce local date to Date format
  Animal$LMT_Date <- as.Date(Animal$LMT_Date, tryFormats = "%m/%d/%Y", tz = "America/New_York")
  
  # filter dates of interest
  Animal.1 <- Animal %>% filter(LMT_Date >= start.date & LMT_Date <= end.date)
  
  # create a timestamp
  Animal.1$datetime <- paste(Animal.1$LMT_Date, Animal.1$LMT_Time)
  Animal.1$timestamp <- as.POSIXct(ymd_hms(Animal.1$datetime), tz = "America/New_York")
  
  # keep only columns we need
  Animal.1 <- Animal.1 %>% dplyr::select(CollarID, timestamp, Latitude...., Longitude...., 
                                         Height..m., DOP, FixType)
  
  # create a track
  animal.track <- make_track(Animal.1, 
                             Longitude...., 
                             Latitude...., 
                             timestamp,
                             crs = 4326,
                             all_cols = TRUE)
  
  # remove rows with NAs 
  animal.track <- animal.track[complete.cases(animal.track),]
  
  # create a telemetry object
  animal.telem <- as_telemetry(animal.track,
                               timeformat = "auto",
                               timezone = "America/New_York",
                               projection = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
                               keep = TRUE)
  
  # pack into a list
  all.animals.list[[x - 37702]] <- animal.telem
  
}

#_____________________________________________________________________________________________________________
# 4. Error modeling ----
#_____________________________________________________________________________________________________________

# fit error model based upon DOP
all.animals.uere <- uere.fit(all.animals.list)

summary(all.animals.uere)

# save for later
save(all.animals.uere, file = "vectronic_uere.RData")
