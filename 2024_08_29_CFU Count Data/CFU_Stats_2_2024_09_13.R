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
library(emmeans)
library(multcomp)


#trim extraneous data and set variable types
str(data)
data.1 <- data %>%
  mutate_at(c("microbe_level"), as.factor)

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

cfu.model.1 <- aov(cfu_conc^0.3 ~ microbe_level, data = data.1)
plot(cfu.model.1)
  #normality looks better on this plot. We transformed by taking the 4th root of cfu_conc

#see the summary
summary(cfu.model.1)
  #we have significance for microbe_level
  #output:
#               Df   Sum Sq   Mean Sq   F value   Pr(>F)    
# microbe_level  3   1334.0   444.7     17.6      1.29e-06 ***
# Residuals     28    707.6   25.3                     


##running post hoc using Tukey test
TukeyHSD(cfu.model.1, conf.level = 0.95)
cld(emmeans(cfu.model.1, pairwise ~ microbe_level, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level          emmean   SE df lower.CL upper.CL .group
# xenic, yes antibiotics   0.00 1.78 28    -3.64     3.64  a    
# xenic, no antibiotics    8.62 1.78 28     4.98    12.26   b   
# OD 0.65                 11.42 1.78 28     7.78    15.06   bc  
# OD 0.32                 17.99 1.78 28    14.35    21.63    c  

#the results:
#OD 0.32 has a higher CFU level than the two xenic treatments
#xenic, yes antibiotics has lower CFU level than all three other treatments
#OD 0.65 has a statistically similar CFU level as both OD 0.32 AND antibiotic-free
  #xenic FAW guts

#just because I am curious, what would the result have been if we did not transform?
TukeyHSD(cfu.model, conf.level = 0.95)
cld(emmeans(cfu.model, pairwise ~ microbe_level, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #the data IS  different. 
    #output:
# microbe_level          emmean   SE df lower.CL upper.CL .group
# xenic, yes antibiotics      0 3704 28    -7587     7587  a    
# xenic, no antibiotics    4007 3704 28    -3580    11593  a    
# OD 0.65                  6355 3704 28    -1231    13942  ab   
# OD 0.32                 19824 3704 28    12238    27411   b   

