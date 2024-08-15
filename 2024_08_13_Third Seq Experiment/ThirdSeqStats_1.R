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

#now for SA-spiked and unspiked diet, looking at SA levels in salivary glands
sa.model <- aov(sa_avg ~ microbe_level*diet_type, data = sa.data)
plot(sa.model)

sa.model.1 <- aov(sa_avg^.3 ~ microbe_level*diet_type, data = sa.data)
plot(sa.model.1)  #winner winner chicken dinner! We shall transform by raising to the power of 3/10,
                    #which is the same as cubing it, then finding the 10th root (or finding the 10th root then cubing it)

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


##posthoc for our two-way ANOVAs using TukeyHSD, followed by compact letter display
#first load the packages needed
library(lme4)
library(emmeans)
library(multcomp)

#post hoc for the BA data
TukeyHSD(ba.model, conf.level=0.95)
cld(emmeans(ba.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean  SE df lower.CL upper.CL .group
# xenic         control     0.00 1.4 18  -2.9320     2.93  a    
# axenic        control     3.02 1.4 18   0.0885     5.95  a    
# E. mundtii    control     5.11 1.4 18   2.1796     8.04  a    
# E. mundtii    benzoic    12.34 1.4 18   9.4130    15.28   b   
# axenic        benzoic    13.14 1.4 18  10.2059    16.07   b   
# xenic         benzoic    13.93 1.4 18  10.9959    16.86   b 

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