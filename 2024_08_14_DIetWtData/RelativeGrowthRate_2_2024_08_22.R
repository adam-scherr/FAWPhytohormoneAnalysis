####Relative Growth Rate for larvae from Sequestration Experiment 3
### Adam Scherr, August 20th, 2024
### Version 1

#import dataset
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_08_14_DIetWtData/Relative%20Growth%20Rate%20Calculations%20-%20Larval%20Weights.csv")

str(data)
data <- data %>%
  mutate_at(c("microbe_level", "diet_type", "tray_number"), as.factor)
str(data)

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
  scale_color_manual(values = c("#9b26b7", "#bf2b6c", "#ce382f"), name = "",
                     labels = c("Axenic", "E. mundtii", "Xenic"))+ 
 labs(x = "Microbe Treatment", y = "Maximum Larval Weight (g)") +
  theme(axis.text = element_text(size =15),
  axis.title = element_text(size=15),
         legend.title = element_blank(),
         legend.text = element_text(size=10))
  # scale_y_continuous(expand = c(0,0), limits = c(-9, 10)) +
  # geom_hline(yintercept = 0, color = "limegreen", linetype = "solid", size = 1)

#days of growth per microbe treatment
ggplot(data, aes(x = microbe_level, y = days_of_growth, color = microbe_level))+
  #geom_boxplot(width = 0.8)+
  geom_jitter(height=0)+
  scale_color_manual(values = c("#9b26b7", "#bf2b6c", "#ce382f"), name = "",
                     labels = c("Axenic", "E. mundtii", "Xenic"))+
  labs(x = "Microbe Treatment", y = "Days to grow to maximum size") +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), limits = c(0, 10))+
  theme(axis.text = element_text(size =15),
        axis.title = element_text(size=15),
        legend.title = element_blank())

#now larvalWt divided by daysOfGrowth
ggplot(data, aes(x = microbe_level, y = larvalWt_div_daysOfGrowth, color = microbe_level))+
  geom_boxplot(width=0.8)+
  geom_point()+
  scale_color_manual(values=c("#9b26b7", "#bf2b6c", "#ce382f"), name = "",
                     labels = c("Axenic", "E. mundtii", "Xenic"))+
  labs(x = "Microbe Treatment", y = "Max. Larval Wt (g)/ Days of Growth")+
theme(axis.text = element_text(size =15),
      axis.title = element_text(size=15),
      legend.title = element_blank(),
      legend.text = element_text(size=10))



####Running Two-Way ANOVA for each microbe_level, based on final larval weight values####
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
larval.wt.model <- aov(final_larval_wt ~ microbe_level*diet_type, data = data)

#check assumptions
par(mfrow = c(2,2))
plot(larval.wt.model)
  #looks gorgeous, variance is similar and residuals are normal

#running a TukeyHSD posthoc test
TukeyHSD(larval.wt.model, conf.level = 0.95)
cld(emmeans(larval.wt.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #statistically significant differences for all possible interactions between
  #the three microbe levels, although when you introduce diet_type, you get more overlap


#####Two-Way ANOVA for each microbe_level, based on larval weight divided by days of growth
growth.rate.model <- aov(larvalWt_div_daysOfGrowth ~ microbe_level*diet_type, data = data)

#check assumptions
plot(growth.rate.model)
  #variance looks similar and residuals are normal

#running a TukeyHSD posthoc test
TukeyHSD(growth.rate.model, conf.level = 0.95)
cld(emmeans(growth.rate.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #statistically significant differences for all possible interactions between 
  #the three microbe levels, slight overlap once you introduce diet_type into the model