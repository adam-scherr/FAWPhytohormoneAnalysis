####Running stats for CFU counts (FAW Gut Microbe Quantification)
### Adam Scherr, August 29th, 2024
### Version 1

#upload dataset
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_08_29_CFU%20Count%20Data/Gut%20Microbe%20Quantification%20-%20Tidy%20CFU%20Counts.csv")

#rename columns for easier coding
colnames(data)[c(3:5)] <- c("cfu_count", "midgut_wt", "cfu_conc")

#load packages
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(ggthemes)
# library(emmeans)
# library(multcomp)

#trim extraneous data and set variable types
str(data)
data.1 <- data %>%
  mutate_at(c("microbe_level"), as.factor) %>%
  filter(microbe_level != "xenic, yes antibiotics") #%>%
 # select(1:5)

#visualize the data
ggplot(data.1, aes(x = microbe_level, y = cfu_conc, color = microbe_level))+
  geom_boxplot()+
  geom_point()

#####Running a one-way ANOVA on these data####
#ANOVA Assumptions:
#1. homogeneity of variance (all variance is roughly equal); homoscedasticity
#2. Independence of observations (this one is true because each point is an independent observation)
#3. Normally-distributed dependent variable
#The second assumption is true, and the first and third assumptions will be tested
#by plotting the residuals of our anova model

#one-way ANOVA for the cfu concentration, as it differs by microbe_level grouping
cfu.model <- aov(cfu_conc ~ microbe_level, data = data.1)

#check assumptions
par(mfrow = c(2,2))
plot(cfu.model)
  #the variance looks similar, but the normality could be improved. Lets try transforming

cfu.model.1 <- aov(cfu_conc^0.25 ~ microbe_level, data = data.1)
plot(cfu.model.1)
  #normality looks better on this plot. We transformed by taking the 4th root of cfu_conc

#see the summary
summary(cfu.model.1)
  #we have significance for one of the interactions

##running post hoc using Tukey test
TukeyHSD(cfu.model.1, conf.level = 0.95)
cld(emmeans(cfu.model.1, pairwise ~ microbe_level, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level         emmean   SE df lower.CL upper.CL .group
# xenic, no antibiotics   5.68 1.18 21     3.23     8.13  a    
# OD 0.65                 7.51 1.18 21     5.06     9.96  ab   
# OD 0.32                11.05 1.18 21     8.60    13.50   b 

#the results:
#OD 0.32 has a higher CFU level than the antibiotic-free xenic FAW guts
#OD 0.65 has a statistically similar CFU level as both OD 0.32 AND antibiotic-free
  #xenic FAW guts

#just because I am curious, what would the result have been if we did not transform?
TukeyHSD(cfu.model, conf.level = 0.95)
cld(emmeans(cfu.model, pairwise ~ microbe_level, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #the result is the same in the untransformed data

