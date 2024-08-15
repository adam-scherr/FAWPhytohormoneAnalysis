####Running stats for diet weight data
### Adam Scherr, August 13th, 2024
### Version 1

#import dataset
data <- All_Weights_Combined_tidy

#load packages
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(ggthemes)

#visualize the data with a quick graph
ggplot(data, aes(x = microbe_level, y = wt, color = microbe_level))+
  geom_boxplot()+
  geom_point()

####Running Two-Way ANOVA for each microbe_level####
##assumptions of a two-way ANOVA:
#1. homogeneity of variance (all variance is roughly equal); homoscedasticity
#2. Independence of observations (this one is true because each point is an independent observation)
#3. Normally-distributed dependent variable
#The second assumption is true, and the first and third assumptions will be tested
#by plotting the residuals of our anova model

#two-way ANOVA for the weights, looking at effect of both microbe_level and diet_type
wt.model <- aov(wt ~ microbe_level*diet_type, data = data)

#check assumptions
par(mfrow = c(2,2))
plot(wt.model)
  #the residuals look evenly above and below the fit line, which shows homoscedasticity
  #the normality is okay, not a perfectly straight line, but there are very few outliers

#running a TukeyHSD posthoc test
TukeyHSD(wt.model, conf.level = 0.95)
  #statistically significant differences for all possible interactions between
  #the three microbe levels. So when only plotting microbe levels, we have SUPER low
  #p values(p = 0 or p < 0.0001).

  #statistical significance is present in the differences between most diets, but not
  #for the difference between 1)IAA and control diet, and 2)SA and control diet

  #and then the interactions accounting for both microbe_level and diet_type are a whole
  #mess, and I think we don't need to report that data. The real actually important
  #data was between microbe levels. We would only care about the effect of diet type
  #on diet consumption if we were trying to compare how well FAW sequestered one 
  #phytohormone compared to how well they sequestered a different phytohormone, given
  #the same gut microbe conditions. This is not a comparison we are looking at.

#two-way ANOVA on SA data
sa.model <- aov(sa_avg ~ microbe_level*diet_type, data = sa.sal)

#check assumptions
par(mfrow = c(2,2))
plot(sa.model)
#also not great equal variance, and not great normality, but once again the small
#number of points messes with these assumptions easily

#SA post hoc
TukeyHSD(sa.model, conf.level=0.95)
#all interactions between SA and control/BA diets are significant (p<0.05)
