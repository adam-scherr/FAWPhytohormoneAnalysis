###SA and IAA Sequestration Data Analysis, version 2
###09/06/2023, Adam Scherr

rm(list=ls())

#upload the csv from Adam Scherr's github
seq_data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/IAA_SA_Seq_Data%20-%20Tidy%20SA%26IAA-Seq.csv")

str(seq_data)
  #there's the problem! the numeric variables are listed as characters in R

#make better names for some columns, using colnames() function
colnames(seq_data)[c(4:12)] <- c("ba_avg", "ba_high", "ba_low",
                                 "sa_avg", "sa_high", "sa_low",
                                 "iaa_avg", "iaa_high", "iaa_low")

#make your numeric variables properly numeric, using mutate_at() in dplyr
library(dplyr)
seq_data <- seq_data %>%
  mutate_at(c("ba_avg", "ba_high", "ba_low",
              "sa_avg", "sa_high", "sa_low",
              "iaa_avg", "iaa_high", "iaa_low"), as.numeric) %>%
  mutate_at(c("microbe_level", "diet_type"), as.factor)
str(seq_data)
  #SCORE! lets see if it translated well into making good plots


#remove the rows containing labeled compounds and the benzoic acid-spiked diet
trimmed_data <- seq_data %>%
filter(id != "",
       diet_type != "benzoic")
View(trimmed_data)

#filter out the "diet" from the rest of the samples
sample_data <- trimmed_data %>%
  filter(microbe_level != "diet")

View(sample_data)

#make a preliminary box and whisker plot to see the data spread
library(ggplot2)
BA_box <- ggplot(sample_data, aes(x = microbe_level, y= ba_avg, 
                                   color = diet_type))+
  geom_boxplot()+
  geom_jitter(size = 3, width = 0.1)
BA_box

ba2 <- ggplot(sample_data, aes(x = id, y = ba_adjusted_conc_sample_ngperg,
                               color = diet_type,
                               shape = microbe_level))+
  geom_point(size = 3)
ba2

SA_box <- ggplot(sample_data, aes(x = microbe_level, y= sa_avg, 
                                  color = diet_type))+
  geom_boxplot()+
  geom_jitter(size = 3, width = 0.1)

IAA_box <- ggplot(sample_data, aes(x = microbe_level, y= iaa_avg, 
                                  color = diet_type))+
  geom_boxplot()+
  geom_jitter(size = 3, width = 0.1)

#arrange all the dot plots in one window
library(gridExtra)
grid.arrange(BA_box, SA_box, IAA_box, ncol=3)

##lets make the individual graphs
library(ggdark) #https://r-charts.com/ggplot2/themes/
library(ggthemes)

#IAA GRAPHS
IAA_point <- ggplot(sample_data, aes(x= id, y = iaa_avg, color = diet_type))+
  geom_point(size = 6)+
  geom_errorbar(aes(ymin = iaa_low,
                    ymax = iaa_high),
                position = position_dodge(0.9), width=0.5, size=1,
                color = "blue")+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan"), name="Diet Type",
                    labels = c("Control", "Indole Acetic Acid", "Salicylic Acid"))+
  labs(x= "", y= "IAA Concentration (ng/g)")+
  theme_clean(base_size=14)
IAA_point


IAA_point2 <- ggplot(sample_data, aes(x= microbe_level, y = iaa_avg, color = diet_type))+
  geom_jitter(size = 6, width=0.2)+
  geom_errorbar(aes(ymin = iaa_low,
                    ymax = iaa_high),
                position = position_dodge(0.9), width=0.3, size=1,
                color = "blue")+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan"), name="Diet Type",
                     labels = c("Control", "Indole Acetic Acid", "Salicylic Acid"))+
  labs(x= "", y= "IAA Concentration (ng/g)")+
  theme_clean(base_size=14)
IAA_point2

IAA_point3 <- ggplot(sample_data, aes(x= microbe_level, y = iaa_avg, color = diet_type))+
  geom_pointrange(aes(ymin = iaa_low, ymax = iaa_high),
                  position=position_jitter(width=0.2), size =1.5)+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan"), name="Diet Type",
                     labels = c("Control", "Indole Acetic Acid", "Salicylic Acid"))+
  labs(x= "", y= "IAA Concentration (ng/g)")+
  theme_clean(base_size=14)
IAA_point3


#BA GRAPHS
BA_point3 <- ggplot(sample_data, aes(x= microbe_level, y = ba_avg, color = diet_type))+
  geom_pointrange(aes(ymin = ba_low, ymax = ba_high),
                  position=position_jitter(width=0.2), size =1.5)+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan"), name="Diet Type",
                     labels = c("Control", "Indole Acetic Acid", "Salicylic Acid"))+
  labs(x= "", y= "BA Concentration (ng/g)")+
  theme_clean(base_size=14)
BA_point3


#SA GRAPHS
SA_point <- ggplot(sample_data, aes(x= id, y = sa_avg, color = diet_type))+
  geom_point(size = 6)+
  geom_errorbar(aes(ymin = sa_low,
                    ymax = sa_high),
                position = position_dodge(0.9), width=0.5, size=1,
                color = "blue")+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan"), name="Diet Type",
                     labels = c("Control", "Indole Acetic Acid", "Salicylic Acid"))+
  labs(x= "", y= "SA Concentration (ng/g)")+
  theme_clean(base_size=14)
SA_point

SA_point2 <- ggplot(sample_data, aes(x= microbe_level, y = sa_avg, color = diet_type))+
  geom_pointrange(aes(ymin = sa_low, ymax = sa_high),
                  position=position_jitter(width=0.2), size =1.5)+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan"), name="Diet Type",
                     labels = c("Control", "Indole Acetic Acid", "Salicylic Acid"))+
  labs(x= "", y= "SA Concentration (ng/g)")+
  theme_clean(base_size=14)
SA_point2
