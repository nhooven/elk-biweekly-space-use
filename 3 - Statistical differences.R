# Title: Biweekly home ranges
# Subtitle: 3 - Statistical differences
# Author: Nathan D. Hooven
# Email: nathan.hooven@uky.edu
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 12 May 2022
# Date completed: 12 May 2022
# Date modified: 
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

all.data <- read.csv("all_data.csv")

#_____________________________________________________________________________________________________________
# 3. Pivot ----
#_____________________________________________________________________________________________________________

all.data.longer <- all.data %>% pivot_longer(5:30)

# remove NAs
all.data.longer <- all.data.longer %>% drop_na()

#_____________________________________________________________________________________________________________
# 4. Look at distributions ----
#_____________________________________________________________________________________________________________

ggplot(data = all.data.longer, aes(x = log(value))) +
       theme_bw() +
       geom_density()

#_____________________________________________________________________________________________________________
# 5. Repeated measures ANOVA ----
#_____________________________________________________________________________________________________________
# 5a. AKDE ----
#_____________________________________________________________________________________________________________

all.data.longer.akde <- all.data.longer %>% filter(Estimator == "AKDE")

akde.aov <- aov(log(value) ~ factor(Area) + Error(factor(Animal)), data = all.data.longer.akde)

summary(akde.aov)

#_____________________________________________________________________________________________________________
# 5b. BBMM ----
#_____________________________________________________________________________________________________________

all.data.longer.bbmm <- all.data.longer %>% filter(Estimator == "BBMM")

bbmm.aov <- aov(log(value) ~ factor(Area) + Error(factor(Animal)), data = all.data.longer.bbmm)

summary(bbmm.aov)
