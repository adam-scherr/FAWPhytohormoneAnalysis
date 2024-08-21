####Relative Growth Rate for larvae from Sequestration Experiment 3
### Adam Scherr, August 20th, 2024
### Version 1

#import dataset
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_08_14_DIetWtData/Relative%20Growth%20Rate%20Calculations%20-%20Larval%20Weights.csv")

#load packages
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(ggthemes)

#visualize the data with a quick graph
#larval weight data first
ggplot(data, aes(x = microbe_level, y = final_larval_wt, color = microbe_level))+
  geom_boxplot(width=0.8)+
  geom_point() +
  scale_color_manual(values = c("#9b26b7", "#bf2b6c", "#ce382f"), name = "") 
 # labs(x = "Microbe Treatment", y = "Avg Diet Eaten (g) per Gram Larval Wt") +
  # theme(axis.text = element_text(size =15),
  #       axis.title = element_text(size=15),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size=10)) +
  # scale_y_continuous(expand = c(0,0), limits = c(-9, 10)) +
  # geom_hline(yintercept = 0, color = "limegreen", linetype = "solid", size = 1)

#now larvalWt divided by daysOfGrowth
ggplot(data, aes(x = microbe_level, y = larvalWt_div_daysOfGrowth, color = microbe_level))+
  geom_boxplot(width=0.8)+
  geom_point()+
  scale_color_manual(values=c("#9b26b7", "#bf2b6c", "#ce382f"), name = "")



####################################everything below this point is just copied from
#DietWtStats_2_2024_08_15; it is NOT yet applicable to Relative Growth Rate data


##find outliers in the data
outliers <- data %>%
  filter(wt > 6 | wt < 0)

#make a new dataset that removes the outlier of -8. That's wild
data.1 <- data %>%
  filter(wt > -7)

#graph the data now that the -8 outlier has been removed
ggplot(data.1, aes(x = microbe_level, y = wt, color = microbe_level))+
  geom_boxplot(width=0.8)+
  geom_point() +
  geom_hline(yintercept = 0, color = "limegreen", linetype = "solid", size = 1)+
  scale_color_manual(values = c("#9b26b7", "#bf2b6c", "#ce382f"), name = "") +
  labs(x = "Microbe Treatment", y = "Avg Diet Eaten (g) per Gram Larval Wt") +
  theme(axis.text = element_text(size =15),
        axis.title = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

##make a new dataset that removes all negative numbers
data.2 <- data %>%
  filter(wt > 0)

#graph the data that has all negatives removed
ggplot(data.2, aes(x = microbe_level, y = wt, color = microbe_level))+
  geom_boxplot(width=0.8)+
  geom_point() +
  geom_hline(yintercept = 0, color = "limegreen", linetype = "solid", size = 1)+
  scale_color_manual(values = c("#9b26b7", "#bf2b6c", "#ce382f"), name = "") +
  labs(x = "Microbe Treatment", y = "Avg Diet Eaten (g) per Gram Larval Wt") +
  theme(axis.text = element_text(size =15),
        axis.title = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

##make a new dataset that replaces the negatives with zeros
data.3 <- data
data.3$wt[data.3$wt < 0] <- 0

#graph the data that has all negatives replaced with zeros
ggplot(data.3, aes(x = microbe_level, y = wt, color = microbe_level))+
  geom_boxplot(width=0.8)+
  geom_point() +
  geom_hline(yintercept = 0, color = "limegreen", linetype = "solid", size = 1)+
  scale_color_manual(values = c("#9b26b7", "#bf2b6c", "#ce382f"), name = "") +
  labs(x = "Microbe Treatment", y = "Avg Diet Eaten (g) per Gram Larval Wt") +
  theme(axis.text = element_text(size =15),
        axis.title = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

####Running Two-Way ANOVA for each microbe_level, using the data WITH -8 OUTLIER REMOVED####
##assumptions of a two-way ANOVA:
#1. homogeneity of variance (all variance is roughly equal); homoscedasticity
#2. Independence of observations (this one is true because each point is an independent observation)
#3. Normally-distributed dependent variable
#The second assumption is true, and the first and third assumptions will be tested
#by plotting the residuals of our anova model

library(lme4)
library(emmeans)
library(multcomp)
#two-way ANOVA for the weights, looking at effect of both microbe_level and diet_type
wt.model.1 <- aov(wt ~ microbe_level*diet_type, data = data.1)

#check assumptions
par(mfrow = c(2,2))
plot(wt.model.1)
  #the residuals look evenly above and below the fit line, which shows homoscedasticity
  #the normality is okay, not a perfectly straight line, but there are very few outliers

#running a TukeyHSD posthoc test
TukeyHSD(wt.model.1, conf.level = 0.95)
cld(emmeans(wt.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #statistically significant differences for all possible interactions between
  #the three microbe levels. So when only plotting microbe levels, we have SUPER low
  #p values(p = 0 or p < 0.0001).

  #statistical significance is present in the differences between all diets EXCEPT
  #for the difference between SA and control diet

  #and then the interactions accounting for both microbe_level and diet_type are a whole
  #mess, and I think we don't need to report that data. The real actually important
  #data was between microbe levels. We would only care about the effect of diet type
  #on diet consumption if we were trying to compare how well FAW sequestered one 
  #phytohormone compared to how well they sequestered a different phytohormone, given
  #the same gut microbe conditions. This is not a comparison we are looking at.

####general finding: removing the outlier DID NOT change the ANOVA results

####Running ANOVA again, but with all negatives removed####
wt.model.2 <- aov(wt ~ microbe_level*diet_type, data = data.2)

#check assumptions
par(mfrow = c(2,2))
plot(wt.model.2)
#the residuals look evenly above and below the fit line, which shows homoscedasticity
#the normality is not the best. Some very large upper outliers

#running a TukeyHSD posthoc test
TukeyHSD(wt.model.2, conf.level = 0.95)
#statistically significant differences for all possible interactions between
#the three microbe levels. So when only plotting microbe levels, we have SUPER low
#p values(p = 0 or p < 0.0001).

#statistical significance is present in the differences between all diets EXCEPT
#for the difference between SA and control diet

####Running ANOVA again, but with all negatives replaced with zeros####
wt.model.3 <- aov(wt ~ microbe_level*diet_type, data = data.3)

#check assumptions
plot(wt.model.3)
#variance looks homoscedastic
#normality could be better, because of those upper outliers

#running TukeyHSD posthoc test
TukeyHSD(wt.model.3, conf.level = 0.95)
#statistically significant differences for all possible interactions between
#the three microbe levels. So when only plotting microbe levels, we have SUPER low
#p values(p = 0 or p < 0.0001).

#statistical significance is present in the differences between all diets EXCEPT
#for the difference between SA and control diet