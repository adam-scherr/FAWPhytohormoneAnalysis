### Running stats for combined third and second sequestration experiment data
### Adam Scherr, August 21th, 2024
### Version 1

#upload the dataset
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_08_21_Combined%202nd%263rd%20Seq%20Exp/1_2nd%20and%203rd%20seq%20combined%20-%20combined_data_1.csv")

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
  mutate_at(c("microbe_level", "diet_type", "tray_number", "seq_experiment"), as.factor)

data.1$diet_type <- factor(data.1$diet_type, 
                           levels = c("control", "benzoic", "salicylic", "indole"))
str(data.1)

#make the negative values into zeros since they are effectively zero
data.1$ba_avg[data.1$ba_avg < 0] <- 0
data.1$sa_avg[data.1$sa_avg < 0] <- 0
data.1$iaa_avg[data.1$iaa_avg < 0] <- 0

#make datasets for each phytohormone to make graphing easier
ba.data <- data.1 %>%
  select(1:8) %>%
  filter(diet_type == "benzoic" | diet_type == "control")

ba.only.data <- ba.data%>%
  filter(diet_type == "benzoic")

sa.data <- data.1 %>%
  select(1:5, 9:11) %>%
  filter(diet_type == "salicylic" | diet_type == "control")

sa.only.data <- sa.data %>%
  filter(diet_type == "salicylic")

iaa.data <- data.1 %>%
  select(1:5, 12:14) %>%
  filter(diet_type=="indole" | diet_type == "control")

iaa.only.data <- iaa.data %>%
  filter(diet_type == "indole")

#### Graphs comparing phytohormone levels for spiked diet-fed FAW, between the two seq experiments####
ba_box <- ggplot(ba.only.data, aes(x = seq_experiment, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"))+
  scale_color_manual(values = c("#9ecae1", "#2171b5", "#08306b")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,17000))
ba_box

sa_box <-  ggplot(sa.only.data, aes(x = seq_experiment, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_color_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,45000))
sa_box

iaa_box <-  ggplot(iaa.only.data, aes(x = seq_experiment, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_color_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,50000))
iaa_box


#####Graphs of just the Unspiked Diet Data####
#first make the unspiked diet data sets
unspiked.ba.data <- ba.data %>%
  filter(diet_type == "control")

unspiked.sa.data <- sa.data %>%
  filter(diet_type == "control")

unspiked.iaa.data <- iaa.data %>%
  filter(diet_type == "control")


#####2 diets, 3 microbe treatments, and 2 sequestration experiments, graphics combined####
all_ba_box <- ggplot(ba.data, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black")+
  geom_point(size = 4, alpha = 1, position=position_dodge(width = 0.5),aes(color = seq_experiment))+
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"))+
  scale_color_manual(values = c("#C594D2", "#A01C1C")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,17000))
# labs(x = "Diet Treatment", y = "BA Concentration (ng/g)") +
#theme_clean() +
# scale_x_discrete(labels = c("Control", "BA-Spiked")) +
# theme(axis.text = element_text(size =15),
#       panel.grid.major = element_line(size = 2),
#       axis.title = element_text(size=20),
#       legend.title = element_blank(),
#       legend.text = element_text(size=15))

all_ba_box

all_sa_box <-  ggplot(sa.data, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(0.8),aes(color = seq_experiment))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_color_manual(values = c("#C594D2", "#A01C1C"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,45000))
all_sa_box

all_iaa_box <-  ggplot(iaa.data, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8),aes(color = seq_experiment))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_color_manual(values = c("#C594D2", "#A01C1C"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,50000))
all_iaa_box

####Running Two-Way ANOVAs on unspiked and spiked data, with both seq. experiments' data####
#this is based on the file entitled "inoc&seq_5_negsReplaced" and on "ThirdSeqStats_2"

##start by testing the assumptions for a two-way ANOVA
#we already know the observations are independent, so you need to see if
#there is homogeneous variance and there is normality

#testing assumptions for BA-spiked and unspiked diet, looking at BA levels in salivary glands
ba.model <- aov(ba_avg ~ microbe_level*diet_type, data = ba.data)
par(mfrow=c(2,2))
plot(ba.model)
  #variance isn't perfect, and the QQ-plot definitely shows poor normality

ba.model.1 <- aov(ba_avg^0.3 ~ microbe_level*diet_type, data = ba.data)
plot(ba.model.1)
  #variance and normality both look good. Transformed by raising to the power of 3/10ths


#teting assujmptions for SA-spiked and unspiked diet, looking at SA levels in salivary glands
sa.model <- aov(sa_avg~microbe_level*diet_type, data = sa.data)
plot(sa.model)
  #variance is okay, normality is quite poor

sa.model.1 <- aov(sa_avg^0.25 ~ microbe_level*diet_type, data = sa.data)
plot(sa.model.1)
  #variance and normality both look good. Transformed by finding the fourth root
  #(or raising to 1/4th power)

#testing assumptions for IAA-spiked and unspiked diet, looking at IAA levels in salivary glands
iaa.model <- aov(iaa_avg ~ microbe_level*diet_type, data= iaa.data)
plot(iaa.model)
  #variance is fine, but normality could be better

iaa.model.1 <- aov(iaa_avg^0.4 ~ microbe_level*diet_type, data= iaa.data)
plot(iaa.model.1)
  #variance could be better, but the normality looks reasonable enough
  #transformed by raising to the power of 4/10ths

##posthoc for our two-way ANOVAs using TukeyHSD, followed by compact letter display
#first load the packages needed
library(lme4)
library(emmeans)
library(multcomp)

#post hoc for the BA data
TukeyHSD(ba.model.1, conf.level=0.95)
cld(emmeans(ba.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
#microbe_level diet_type emmean    SE df lower.CL upper.CL .group
# xenic         control     0.00 1.286 30    -2.63     2.63  a    
# axenic        control     3.41 0.972 30     1.43     5.40  ab   
# E. mundtii    control     5.19 0.972 30     3.20     7.18   b   
# axenic        benzoic    13.64 0.972 30    11.66    15.63    c  
# xenic         benzoic    13.93 1.286 30    11.30    16.55    c  
# E. mundtii    benzoic    14.46 0.972 30    12.47    16.44    c  

#post hoc for SA data
TukeyHSD(sa.model.1, conf.level=0.95)
cld(emmeans(sa.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean    SE df lower.CL upper.CL .group
# axenic        control    0.131 0.549 30   -0.989     1.25  a    
# xenic         control    1.061 0.726 30   -0.422     2.54  a    
# E. mundtii    control    2.406 0.549 30    1.285     3.53  a    
# axenic        salicylic 10.624 0.549 30    9.503    11.74   b   
# E. mundtii    salicylic 11.572 0.549 30   10.451    12.69   b   
# xenic         salicylic 12.045 0.726 30   10.562    13.53   b  

#post hoc for IAA data
TukeyHSD(iaa.model.1, conf.level = 0.95)
cld(emmeans(iaa.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters = letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean   SE df lower.CL upper.CL .group
# E. mundtii    control     2.18 3.49 30    -4.95     9.32  a    
# axenic        control     3.02 3.49 30    -4.11    10.16  a    
# xenic         control     3.70 4.62 30    -5.75    13.14  a    
# E. mundtii    indole     43.19 3.49 30    36.06    50.33   b   
# xenic         indole     52.88 4.62 30    43.44    62.32   b   
# axenic        indole     55.47 3.49 30    48.33    62.60   b 