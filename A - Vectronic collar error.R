# Title: Raw data processing
# Subtitle: 4 - Vectronic collar location error
# Author: Nathan D. Hooven
# Email: nathan.hooven@uky.edu
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 16 Feb 2021
# Date completed: 16 Feb 2021
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

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

# collars to include
elk.names.2020 <- as.vector(c("37703", "37704", "37705", "37706", "37707", "37708", "37709", "37710", "37711",
                              "37712", "37713", "37714", "37715", "37716", "37717", "37718", "37719", "37720",
                              "37721", "37722", "37723", "37724", "37725", "37726", "37727"))

# 2021 collars

#input your export directory as per GPS Plus X
exp.dir <- "C:/Users/Public/Documents/VAS/GPS Plus X/AutoExport/" 

# Create blank data frames to hold all data
elk.dist <- data.frame()

#_____________________________________________________________________________________________________________
# 3. 2020 - Robinson forest tests ----
#_____________________________________________________________________________________________________________

# Robinson forest test location
RF.test <- data.frame(x = 309154.06,
                      y = 4148138.30)

for (x in elk.names.2020){
  
  #_____________________________________________________________________________________________________________
  # 3a. Select data ----
  #_____________________________________________________________________________________________________________
  CollarID <- x
  
  Animal  <- read.csv(paste0(exp.dir, "Collar", CollarID, "_GPS_Default Storage.csv"))
  
  # Timeframe
  Start.Date <- as.Date("11/18/2019", tryFormats = "%m/%d/%Y")
  End.Date <- as.Date("11/22/2019", tryFormats = "%m/%d/%Y")
  
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
  # 3c. Determine distance (in m) from locations ----
  #_____________________________________________________________________________________________________________
  
  # define projection
  projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
  
  # delete "17" zone from beginning of Easting string
  Animal.1$Easting <- str_sub(Animal.1$Easting, 3, 11)
  Animal.1$Easting <- as.numeric(Animal.1$Easting)
  
  # Create a SpatialPointsDataFrame to hold all distances
  Animal.spdf <- SpatialPointsDataFrame(coords = Animal.1[ ,c("Easting", "Northing")],
                                        proj4string = projection,
                                        data = Animal.1[ ,c("Easting", "Northing")])
  
  Animal.spdf@data$Dist <- NULL
  
  # calculate straight line distances from each point to the actual location
  for(i in 1:nrow(Animal.spdf)) {
    
    loc <- Animal.spdf@data[i,]
    
    dist <- sqrt((loc[1] - RF.test[1])^2 + (loc[2] - RF.test[2])^2)
    
    Animal.spdf@data$Dist[i] <- as.numeric(dist)
    
  }
  
  # create animal ID column
  Animal.spdf@data$Animal <- CollarID
  
  # write to master data frame
  elk.dist <- rbind(elk.dist, Animal.spdf@data)
  
}

#_____________________________________________________________________________________________________________
# 4. Summarize error ----
#_____________________________________________________________________________________________________________

ggplot(data = elk.dist, aes(x = Dist)) +
  geom_density(fill = "blue", alpha = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = mean(elk.dist$Dist), size = 1.25, color = "blue") +
  geom_vline(xintercept = mean(elk.dist$Dist) + sd(elk.dist$Dist), size = 1.25, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = mean(elk.dist$Dist) - sd(elk.dist$Dist), size = 1.25, color = "blue", linetype = "dashed")

mean(elk.dist$Dist, na.rm = TRUE)
sd(elk.dist$Dist, na.rm = TRUE)