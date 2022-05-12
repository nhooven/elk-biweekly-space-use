# Title: Raw data processing
# Subtitle: 3 - Remove mort locations
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

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

vectronic.relocations <- read.csv("Relocations_vectronic.csv")
lotek.relocations <- read.csv("Relocations_lotek.csv")

time.zone <- "America/New_York"

# convert timestamps to POSIXct
vectronic.relocations$t <- as.POSIXct(vectronic.relocations$t, tz = time.zone)
lotek.relocations$t <- as.POSIXct(lotek.relocations$t, tz = time.zone)

#_____________________________________________________________________________________________________________
# 3. Remove mortality/post-mortality relocations ----
#_____________________________________________________________________________________________________________

# relocations
vectronic.relocations <- vectronic.relocations %>% filter(!((Animal == "37717" & t > as.POSIXct("2020-10-31 00:00:00", tz = time.zone)) | 
                                                           ((Animal == "37721" & t > as.POSIXct("2020-04-02 00:00:00", tz = time.zone)) |
                                                           ((Animal == "37727" & t > as.POSIXct("2020-12-29 00:00:00", tz = time.zone)) |
                                                           ((Animal == "46395" & t > as.POSIXct("2021-05-19 00:00:00", tz = time.zone)) |
                                                           ((Animal == "45495" & t > as.POSIXct("2021-09-07 00:00:00", tz = time.zone)) |
                                                           ((Animal == "45499" & t > as.POSIXct("2021-09-23 00:00:00", tz = time.zone)) |
                                                           ((Animal == "45509" & t > as.POSIXct("2021-11-05 00:00:00", tz = time.zone)) |
                                                            (Animal == "46401" & t > as.POSIXct("2021-12-27 00:00:00", tz = time.zone))))))))))

lotek.relocations <- lotek.relocations %>% filter(!((Animal == "101959" & t > as.POSIXct("2021-11-01 00:00:00", tz = time.zone)) |
                                                   ((Animal == "101972" & t > as.POSIXct("2020-06-05 00:00:00", tz = time.zone)) | 
                                                   ((Animal == "102489" & t > as.POSIXct("2020-11-28 00:00:00", tz = time.zone)) |
                                                   ((Animal == "102492" & t > as.POSIXct("2020-09-19 00:00:00", tz = time.zone)) |
                                                   ((Animal == "102494" & t > as.POSIXct("2020-03-28 00:00:00", tz = time.zone)) |
                                                   ((Animal == "102495" & t > as.POSIXct("2020-11-29 00:00:00", tz = time.zone)) |
                                                   ((Animal == "103183" & t > as.POSIXct("2021-01-27 00:00:00", tz = time.zone)) |
                                                    (Animal == "103184" & t > as.POSIXct("2021-12-01 00:00:00", tz = time.zone))))))))))

# write to csvs
write.csv(vectronic.relocations, "Relocations_vectronic_1.csv")
write.csv(lotek.relocations, "Relocations_lotek_1.csv")
