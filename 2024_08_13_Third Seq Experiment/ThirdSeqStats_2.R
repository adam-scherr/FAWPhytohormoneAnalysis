### Making Graphics for Third (Final) Sequestration Experiment
### Adam Scherr, August 12th, 2024
### Version 1

#upload the dataset
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_08_13_Third%20Seq%20Experiment/1_2024_07_25%20GCMS%20Data%20-%20tidy_data.csv")

#load in packages
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(ggthemes)

#### Setting up datasets ####
#improve value types and add factor levels
str(data)
data.1 <- data %>%
  select(!3) %>%
  mutate_at(c("microbe_level", "diet_type"), as.factor) %>%
  filter(microbe_level != "special")

data.1$diet_type <- factor(data.1$diet_type, 
                           levels = c("control", "benzoic", "salicylic", "indole"))
str(data.1)

#make the negative values into zeros since they are effectively zero
data.1$ba_avg[data.1$ba_avg < 0] <- 0
data.1$sa_avg[data.1$sa_avg < 0] <- 0
data.1$iaa_avg[data.1$iaa_avg < 0] <- 0

#make datasets for each phytohormone to make graphing easier
ba.data <- data.1 %>%
  select(1:6) %>%
  filter(diet_type == "benzoic" | diet_type == "control")

ba.sa.data <- data.1 %>%
  filter(diet_type != "indole")

sa.data <- data.1 %>%
  select(1:3, 7:9) %>%
  filter(diet_type == "salicylic" | diet_type == "control")

iaa.data <- data.1 %>%
  select(1:3, 10:12) %>%
  filter(diet_type=="indole" | diet_type == "control")

#### Preliminary graphs ####
ba_box <- ggplot(ba.data, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"), name = "")+
  scale_color_manual(values = c("#9ecae1", "#2171b5", "#08306b"), name ="")
ba_box

sa_box <-  ggplot(sa.data, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"), name = "")+
  scale_color_manual(values = c("#a1d99b", "#41ab5d", "#00441b"), name ="")
sa_box

iaa_box <-  ggplot(iaa.data, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"), name = "")+
  scale_color_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"), name ="")
iaa_box

####Two-Way ANOVA for combined unspiked and spiked diet data####
#this is based on the file entitled "inoc&seq_5_negsReplaced"

##start by testing the assumptions for a two-way ANOVA
#we already know the observations are independent, so you need to see if
  #there is homogeneous variance and there is normality

#testing assumptions for BA-spiked and unspiked diet, looking at BA levels in salivary glands
ba.model <- aov(ba_avg^.3 ~ microbe_level*diet_type, data = ba.data)
par(mfrow=c(2,2))
plot(ba.model)
  #the variance is pretty well distributed, and the normality is not perfect but I think still sufficient
summary(ba.model)

####whacky fun time with Jared Adam ####
lm1 <- lm(ba_avg^0.3 ~ microbe_level*diet_type, data = ba.data)
summary(lm1)
anova(lm1)
summary(ba.model)
lm2 <- lm(ba_avg^0.3 ~ microbe_level*diet_type+id, data = ba.data)
summary(lm2)
install.packages("effects")
library(effects)
plot(allEffects(lm1))
plot(allEffects(lm2))

#now for BA-spiked, SA-spiked and unspiked diet effect on SA levels
ba.sa.model <- aov(sa_avg^0.3 ~ microbe_level*diet_type, data = ba.sa.data)
plot(ba.sa.model) #looks good now that it's transformed
    #transformed by raising sa_avg to power of 3/10


#now for SA-spiked and unspiked diet, looking at SA levels in salivary glands
sa.model <- aov(sa_avg ~ microbe_level*diet_type, data = sa.data)
plot(sa.model)

sa.model.1 <- aov(sa_avg^.3 ~ microbe_level*diet_type, data = sa.data)
plot(sa.model.1)  #winner winner chicken dinner! We shall transform by raising to the power of 3/10,
                    #which is the same as cubing it, then finding the 10th root (or finding the 10th root then cubing it)
summary(sa.model.1)
# sa.model.2 <- aov(sa_avg^.2 ~ microbe_level*diet_type, data = sa.data)
# plot(sa.model.2)

  #variance is pretty well distributed, but the normality is sort of rocky. May need to transform for normality

#now for IAA-spiked and unspiked diet, looking at IAA levels in salivary glands
iaa.model <- aov(iaa_avg ~ microbe_level*diet_type, data = iaa.data)
plot(iaa.model)
  #variance is good, normality is about the same as it was for BA. Workable

#here's an attempt to imrpove the normality
iaa.model.1 <- aov(iaa_avg^0.5 ~ microbe_level*diet_type, data = iaa.data)
plot(iaa.model.1) #I think that's as good as we're going to get, not a real improvement
summary(iaa.model.1)

##posthoc for our two-way ANOVAs using TukeyHSD, followed by compact letter display
#first load the packages needed
library(lme4)
library(emmeans)
library(multcomp)

#post hoc for the BA data
TukeyHSD(ba.model, conf.level=0.95)
cld(emmeans(ba.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
ba.model
  #output:
# microbe_level diet_type emmean  SE df lower.CL upper.CL .group
# xenic         control     0.00 1.4 18  -2.9320     2.93  a    
# axenic        control     3.02 1.4 18   0.0885     5.95  a    
# E. mundtii    control     5.11 1.4 18   2.1796     8.04  a    
# E. mundtii    benzoic    12.34 1.4 18   9.4130    15.28   b   
# axenic        benzoic    13.14 1.4 18  10.2059    16.07   b   
# xenic         benzoic    13.93 1.4 18  10.9959    16.86   b 


#post hoc for SA levels in BA-spiked, SA-spiked, and unspiked diet, transformed version
TukeyHSD(ba.sa.model, conf.level = 0.95)
cld(emmeans(ba.sa.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean   SE df lower.CL upper.CL .group
# axenic        control    0.226 1.11 27   -2.061     2.51  a    
# E. mundtii    benzoic    0.749 1.11 27   -1.538     3.04  a    
# xenic         control    1.233 1.11 27   -1.053     3.52  a    
# xenic         benzoic    3.025 1.11 27    0.738     5.31  a    
# axenic        benzoic    3.075 1.11 27    0.788     5.36  a    
# E. mundtii    control    3.349 1.11 27    1.063     5.64  a    
# axenic        salicylic 15.005 1.11 27   12.718    17.29   b   
# E. mundtii    salicylic 17.476 1.11 27   15.189    19.76   b   
# xenic         salicylic 19.843 1.11 27   17.556    22.13   b 

#post hoc for SA data, using the transformed version of the data (sa.model.1)
TukeyHSD(sa.model.1, conf.level=0.95)
cld(emmeans(sa.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha=0.05)
  #output:
# microbe_level diet_type emmean   SE df lower.CL upper.CL .group
# axenic        control    0.226 1.21 18   -2.316     2.77  a    
# xenic         control    1.233 1.21 18   -1.309     3.78  a    
# E. mundtii    control    3.349 1.21 18    0.807     5.89  a    
# axenic        salicylic 15.005 1.21 18   12.463    17.55   b   
# E. mundtii    salicylic 17.476 1.21 18   14.934    20.02   b   
# xenic         salicylic 19.843 1.21 18   17.301    22.38   b   

#just curious, are the results any different if I don't transform the data (sa.model)
TukeyHSD(sa.model, conf.level=0.95)
cld(emmeans(sa.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha=0.5)
  #output: there actually is a difference here
# microbe_level diet_type   emmean   SE df lower.CL upper.CL .group
# axenic        control   1.77e-01 3265 18    -6859     6859  a    
# xenic         control   1.03e+01 3265 18    -6848     6869  a    
# E. mundtii    control   8.81e+01 3265 18    -6771     6947  a    
# axenic        salicylic 9.52e+03 3265 18     2660    16377   b   
# E. mundtii    salicylic 1.53e+04 3265 18     8484    22201   bc  
# xenic         salicylic 2.26e+04 3265 18    15709    29426    c  

#post hoc for IAA, using untransformed data (iaa.model)
TukeyHSD(iaa.model, conf.level=0.95)
cld(emmeans(iaa.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha=0.5)
  #output:
# microbe_level diet_type   emmean   SE df lower.CL upper.CL .group
# axenic        control       0.00 1985 18    -4170     4170  a    
# E. mundtii    control       9.51 1985 18    -4160     4179  a    
# xenic         control     102.52 1985 18    -4067     4272  a    
# E. mundtii    indole     6716.37 1985 18     2547    10886   b   
# axenic        indole    13414.21 1985 18     9245    17584    c  
# xenic         indole    21176.73 1985 18    17007    25346     d 

#let's see if it is any different with the slightly more normal, transformed data (iaa.model.1)
TukeyHSD(iaa.model.1, conf.level=0.95)
cld(emmeans(iaa.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha=0.5)
#output: the significance levels are the same as the untransformed data
# microbe_level diet_type emmean   SE df lower.CL upper.CL .group
# axenic        control     0.00 6.73 18   -14.14     14.1  a    
# E. mundtii    control     1.54 6.73 18   -12.60     15.7  a    
# xenic         control     6.29 6.73 18    -7.85     20.4  a    
# E. mundtii    indole     81.85 6.73 18    67.71     96.0   b   
# axenic        indole    115.74 6.73 18   101.60    129.9    c  
# xenic         indole    143.06 6.73 18   128.92    157.2     d 

####One-Way ANOVA for unspiked diet only data####
#start by making the unspiked diet datasets
unspiked.ba.data <- ba.data %>%
  filter(diet_type == "control")

unspiked.sa.data <- sa.data %>%
  filter(diet_type == "control")

unspiked.iaa.data <- iaa.data %>%
  filter(diet_type == "control")

##make ANOVA models and check for significance

#first for BA
unspiked.ba.model <- aov(ba_avg ~ microbe_level, data= unspiked.ba.data)
plot(unspiked.ba.model)
  #normality and homoscedasticity look pretty good
summary(unspiked.ba.model)
  #nearly significant, p = 0.0586

#now for SA
unspiked.sa.model <- aov(sa_avg ~ microbe_level, data = unspiked.sa.data)
plot(unspiked.sa.model) #not the best normality, but workable
summary(unspiked.sa.model)
  #not significant, p = 0.104

#trying for a better normal distribution transformation in SA data
unspiked.sa.model.1 <- aov(sa_avg^.4 ~ microbe_level, data = unspiked.sa.data)
plot(unspiked.sa.model.1)
  #transforming by putting it to the power of 0.4 makes for better normality
summary(unspiked.sa.model.1)
  #AND this model DOES give us significance. Interesting

#now for IAA
unspiked.iaa.model <- aov(iaa_avg ~ microbe_level, data = unspiked.iaa.data)
plot(unspiked.iaa.model)
  #not the best homoscedasticity or normality
summary(unspiked.iaa.model) #not significant, p = 0.365

#trying for better assumptions in IAA data
unspiked.iaa.model.1 <- aov(iaa_avg^0.2 ~ microbe_level, data = unspiked.iaa.data)
plot(unspiked.iaa.model.1)  
  #this has slightly better normality, but the homoscedasticity is about the same
summary(unspiked.iaa.model.1)

##post hoc tests for these models
#BA posthoc
TukeyHSD(unspiked.ba.model, level=0.95)
cld(emmeans(unspiked.ba.model, pairwise ~ microbe_level, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level emmean  SE df lower.CL upper.CL .group
# xenic              0 114  9   -256.8      257  a    
# axenic           205 114  9    -51.6      462  ab   
# E. mundtii       451 114  9    194.0      708   b   

#SA posthoc
TukeyHSD(unspiked.sa.model.1, level=0.95)
cld(emmeans(unspiked.sa.model.1, pairwise~microbe_level, adjust="tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level emmean   SE df lower.CL upper.CL .group
# axenic         0.218 1.02  9   -2.084     2.52  a    
# xenic          1.668 1.02  9   -0.634     3.97  ab   
# E. mundtii     5.176 1.02  9    2.874     7.48   b 

#just curious, is the output the same for the untransformed SA data?
TukeyHSD(unspiked.sa.model, level=0.95)
cld(emmeans(unspiked.sa.model, pairwise~microbe_level, adjust="tukey"),
    Letters = letters, alpha = 0.05)
#NO, in the untransformed data, all 3 groups are statistically the same

#IAA posthoc
TukeyHSD(unspiked.iaa.model.1, level=0.95)
cld(emmeans(unspiked.iaa.model.1, pairwise~microbe_level, adjust="tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level emmean    SE df lower.CL upper.CL .group
# axenic         0.000 0.553  9  -1.2507     1.25  a    
# E. mundtii     0.518 0.553  9  -0.7331     1.77  a    
# xenic          1.322 0.553  9   0.0716     2.57  a  

#just curious, is the output the same for the untransformed IAA data?
TukeyHSD(unspiked.iaa.model, level=0.95)
cld(emmeans(unspiked.iaa.model, pairwise~microbe_level, adjust="tukey"),
    Letters = letters, alpha = 0.05)
#YES, the output is the same for the transformed versus the untransformed IAA data
