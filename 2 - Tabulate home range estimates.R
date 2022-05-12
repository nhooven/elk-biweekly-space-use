# Title: Biweekly home ranges
# Subtitle: 2 - Tabulate home range estimates
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
# 3. Group and summarize ----
#_____________________________________________________________________________________________________________

all.data.longer <- all.data %>% pivot_longer(5:30)

# remove NAs
all.data.longer <- all.data.longer %>% drop_na()

all.data.summary <- all.data.longer %>% group_by(Estimator, Area, name) %>%
                                        summarize(mean = mean(value, na.rm = TRUE),
                                                  sd = sd(value, na.rm = TRUE),
                                                  n = n()) %>%
                                        mutate(Period.num = as.numeric(substr(name, 8, 9))) %>%
                                        arrange(Period.num, by_group = TRUE) %>%
                                        dplyr::select(-name)

write.table(all.data.summary, "clipboard", sep = "\t")                                                

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
# 4c. Merge into summary dataset ----
#_____________________________________________________________________________________________________________

colnames(periods.2020)[1] <- "Period.num"

all.data.summary <- all.data.summary %>% left_join(periods.2020, by = "Period.num")

#_____________________________________________________________________________________________________________
# 5. Create plot ----
#_____________________________________________________________________________________________________________

# create a labeling vector
week.label <- substr(periods.2020$Start, 6, 10)

all.data.summary$Period.fact <- factor(all.data.summary$Period.num,
                                       labels = week.label)

# compute 95% CIs
all.data.summary <- all.data.summary %>% mutate(SE = sd/sqrt(n)) %>% 
                                         mutate(lwr = mean - 1.96*SE,
                                                upr = mean + 1.96*SE)

ggplot(data = all.data.summary, aes(x = Period.fact, y = mean)) +
       facet_wrap(~Estimator, nrow = 2, scales = "free_y") +
       theme_bw() +
       geom_line(aes(color = Estimator, 
                     linetype = Area,
                     group = Estimator:Area),
                  position = position_dodge(width = 1),
                 size = 1.05) +
       scale_color_manual(values = c("darkorange1", "dodgerblue")) +
       scale_fill_manual(values = c("darkorange1", "dodgerblue")) +
       geom_point(aes(fill = Estimator,
                      shape = Area),
                  position = position_dodge(width = 1),
                  size = 2) +
       scale_shape_manual(values = c(21, 24)) +
       geom_errorbar(aes(ymin = lwr,
                         ymax = upr,
                         color = Estimator, 
                         linetype = Area),
                      position = position_dodge(width = 1),
                     width = 0) +
       theme(legend.position = "none",
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             axis.text.x = element_text(angle = 90, vjust = 0.5),
             axis.title.x = element_blank()) +
       ylab(expression("Mean area " (km^2)))
