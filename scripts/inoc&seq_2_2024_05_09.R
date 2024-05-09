### Analysis of phytohormone contents in salivary glands of bacteria-inoculated FAW
### Adam Scherr, May 9th, 2024
### Version 2


#upload your datasets
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_05_08_Bacterial%20Inoculation%26Sequestration/Inoc_Seq_Data.csv")

#load in packages
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(ggthemes)

####setting up datasets####
#change the object types so the right things are characters/factors
str(data)
data.1 <- data %>%
  select(!3) %>%
  mutate_at(c("microbe_level", "diet_type"), as.factor)
str(data.1)

#make datasets for each phytohormone in spiked diet
ba.data <- data.1 %>%
  select(1:6) %>%
  filter(diet_type == "benzoic")

sa.data <- data.1 %>%
  select(1:3, 7:9) %>%
  filter(diet_type == "salicylic")

iaa.data <- data.1 %>%
  select(1:3, 10:12) %>%
  filter(diet_type=="indole")

#datasets for each phytohormone in unspiked (control) diet
unspiked.ba.data <- data.1 %>%
  select(1:6) %>%
  filter(diet_type == "control")

unspiked.sa.data <- data.1 %>%
  select(1:3, 7:9)%>%
  filter(diet_type == "control")

unspiked.iaa.data <- data.1 %>%
  select(1:3, 10:12)%>%
  filter(diet_type == "control")

####graphs of the data####
##lets see what these datasets look like graphically
  #still need to have axis titles changes,legends removed, and theme added
ba_box <- ggplot(ba.data, aes(x = microbe_level, y=ba_avg, fill=microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#3182bd", "#08306b"))+
  scale_color_manual(values = c("#2171b5", "#2171b5", "#2171b5"))
ba_box

sa_box <- ggplot(sa.data, aes(x = microbe_level, y = sa_avg, fill = microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b", "#238b45", "#00441b"))+
  scale_color_manual(values = c("#64bc6e", "#64bc6e", "#64bc6e"))
sa_box

iaa_box <- ggplot(iaa.data, aes(x=microbe_level, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b", "#d94801", "#7f2704"))+
  scale_color_manual(values = c("#fd8d3c", "#fd8d3c", "#fd8d3c"))
iaa_box


##graphs of the unspiked diet data
unspiked_ba_box <- ggplot(unspiked.ba.data, aes(x = microbe_level, y=ba_avg, fill=microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#3182bd", "#08306b"))+
  scale_color_manual(values = c("#2171b5", "#2171b5", "#2171b5"))

unspiked_sa_box <- ggplot(unspiked.sa.data, aes(x = microbe_level, y = sa_avg, fill = microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b", "#238b45", "#00441b"))+
  scale_color_manual(values = c("#64bc6e", "#64bc6e", "#64bc6e"))

unspiked_iaa_box <- ggplot(unspiked.iaa.data, aes(x=microbe_level, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b", "#d94801", "#7f2704"))+
  scale_color_manual(values = c("#fd8d3c", "#fd8d3c", "#fd8d3c"))

unspiked_ba_box
unspiked_sa_box
unspiked_iaa_box
####statistical analysis_ANOVA####
#Use ANOVA because we have one categorical explanatory variable (microbe_level) and want to 
  #compare continuous response variable (phytohormone concentration) between 3 
#groups (axenic, Enterococcus mundtii -inoculated, and Frass Community-inoculated)

##testing ANOVA assumptions
#resource: https://www.statology.org/one-way-anova-r/

##checking assumptions for BA data
ba.model <- aov(ba_avg ~ microbe_level, data = ba.data)

#QQ plot for BA
par(mfrow = c(2, 2))
plot(ba.model)
 #QQ-plot looks like it fits the diagonal line very well, data is normal
 #the distribution looks roughly equal, but a Levene's test will show for certain

#Levene's test for equal variances for BA
library(car)
leveneTest(ba_avg ~ microbe_level, data = ba.data)
  #p>0.05, so we DO have equal variance, we passed both our assumptions for BA

##running ANOVA for BA data
summary(ba.model)
  #p = 0.0802, nearly significant for p<0.05 significance
  #I think I should run the Tukey-HSD posthoc test anyway to see where the significance is


##checking assumptions for SA data
sa.model <- aov(sa_avg ~ microbe_level, data = sa.data)

#QQ-plot for SA
par(mfrow = c(2, 2))
plot(sa.model)
  #QQ-plot is not perfectly fitting along diagonal line, let's double check with Shapiro-Wilk test

#Shapiro-Wilk test for SA
shapiro.test(sa.data$sa_avg)
  #p=0.8847, far above 0.05. We have normality

#Levene's test for equal variances for SA
leveneTest(sa_avg~microbe_level, data = sa.data)
  #p=0.7389, far above 0.05. Variances are equal

##running ANOVA for SA data
summary(sa.model)
  #p=0.748. There is no significance of microbe_level on SA concentration


##checking assumptions for IAA data
iaa.model <- aov(iaa_avg~microbe_level, data = iaa.data)

##QQ-plot for IAA
par(mfrow = c(2, 2))
plot(iaa.model)

#Shapiro-Wilk test for IAA
shapiro.test(iaa.data$iaa_avg)
  #p = 0.2333, above 0.05 so we have normality

#Levene's test for IAA
leveneTest(iaa_avg ~ microbe_level, data = iaa.data)
  #p=0.4306, above 0.05 so variances are equal

##running ANOVA for IAA data
summary(iaa.model)
  #p=0.012, below 0.05, microbe_level is SIGNIFICANT for IAA concentration

####statistical analysis_Post Hoc Tukey's Test####
##post-hoc for BA data (nearly significant p-value)
TukeyHSD(ba.model, conf.level=0.95)
  #p>0.05 for all comparisions.NO SIGNIFIANCE. The closest to significance is p=0.097 
  #for interaction between E. mundtii and axenic AND p=0.13 for interaction between 
  #E. mundtii and Frass Community

##post-hoc for SA data
TukeyHSD(sa.model, conf.level=0.95)
  #p is far greater than 0.05. NO SIGNIFICANCE

#post-hoc for IAA data
TukeyHSD(iaa.model, conf.level=0.95)
  #p=0.0099 (SIGNIFICANT) for interaction between axenic and Frass Community
  #p=0.11 (close to significant) for interaction between axenic and E. mundtii
