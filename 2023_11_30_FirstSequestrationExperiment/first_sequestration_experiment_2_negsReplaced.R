###Analysis of BA, SA, and IAA Sequestration data, controlled study
####adapted from Second Committee Meeting Code sessions
###Adam Scherr, 11/27/2023
###Version 1, 05/09/2023

####Introduction and Data Sets####
#upload the dataset
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2023_11_30_CommitteeMeeting/EntSoc_SeqData_2%20-%20tidy_data.csv")

#rename the columns so they are easier to type
colnames(data)[c(9:17)] <- c("ba_avg", "ba_high", "ba_low",
                             "sa_avg", "sa_high", "sa_low",
                             "iaa_avg", "iaa_high", "iaa_low")

#check the format of each variable
str(data) #a lot of characters, lets make them numeric

#lets make things the correct object and tidy the data at the same time
#we should have numeric and factor objects, and eliminate the blanks and labeled compounds
library(dplyr)
data.1 <- data %>%
  mutate_at(c("conc_sample_ngperg", "ba_avg", "ba_high", "ba_low",
              "sa_avg", "sa_high", "sa_low",
              "iaa_avg", "iaa_high", "iaa_low"), as.numeric) %>%
  mutate_at(c("microbe_level", "diet_type"), as.factor) %>%
  filter(microbe_level != "blank", nat_lab != "lab") %>%
  filter(id != "NC2H_1")
str(data.1)

data.1$diet_type <- factor(data.1$diet_type, 
                           levels = c("control", "benzoic", "salicylic", "indole"))
data.1$tissue <- factor(data.1$tissue,
                        levels = c("sal", "fat", "hem"))

#now make the negative values into zeros
data.1$ba_avg[data.1$ba_avg < 0] <- 0
data.1$sa_avg[data.1$sa_avg < 0] <- 0
data.1$iaa_avg[data.1$iaa_avg < 0] <- 0

data.2 <- data.1 %>%
  select(1:4, 9, 12, 15)

data.sal <- data.2 %>%
  filter(tissue == "sal")



##set up datasets for each phytohormone-control diet comparison
#BA and control
ba.sal <- data.sal %>%
  filter(diet_type != "indole") %>%
  filter(diet_type != "salicylic")

#SA and control
sa.sal <- data.sal %>%
  filter(diet_type != "indole")
 # filter(diet_type != "benzoic")

#IAA and control
iaa.sal <- data.sal %>%
  filter(diet_type != "benzoic") %>%
  filter(diet_type != "salicylic")

#example graph
library(ggplot2)
ba.box <- ggplot(ba.sal, aes(x = diet_type, y=ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#3182bd", "#08306b"))+
  scale_color_manual(values = c("#2171b5", "#2171b5", "#2171b5")) + 
  theme_bw()
ba.box

####lets try running ANOVA on this####
#resource: http://www.sthda.com/english/wiki/two-way-anova-test-in-r#google_vignette
#goal: use a two-way anova to show how microbe_level and diet_type affect phytohormone concentration
library(lme4)
library(emmeans)
library(multcomp)
##using BA as first two-way ANOVA example
ba.model <- aov(ba_avg ~ microbe_level*diet_type, data = ba.sal)
summary(ba.model)
plot(ba.model)

#trying to transform the data into a better form
ba.sal$scaled.ba_avg <- as.numeric(scale(ba.sal$ba_avg)) ## scaling response (subtract by mean, and diving by SD)
ba.model.scale <- aov(exp(scaled.ba_avg) ~ microbe_level*diet_type, data = ba.sal)
summary(ba.model.scale)

ba.model.lm <- glm(ba_avg ~ microbe_level*diet_type, data = ba.sal)
summary(ba.model.lm)
hist(residuals(ba.model.lm))
cld(emmeans(ba.model.lm,~microbe_level*diet_type))

ba.model.square <- aov(ba_avg^2 ~ microbe_level*diet_type, data = ba.sal)
plot(ba.model.square)
summary(ba.model.square)
TukeyHSD(ba.model.square, conf.level=0.95)

ba.glm.squared <- glm(ba_avg^2 ~ microbe_level * diet_type, data = ba.sal)
summary(ba.glm.squared)
plot(ba.glm.squared)
cld(emmeans(ba.glm.squared,~microbe_level*diet_type, data = ba.sal))
#check our assumptions
par(mfrow = c(2, 2))
plot(ba.model.scale)

#QQ-plot and variance don't look great
#We have so few data points that it's easy to skew the normality and variance
#So let's go through with it since our small number of points is skewing normality and variance so much
summary(ba.model)
  #p=0.0261 for diet_type. Only diet_type has SIGNIFICANT effect on phytohormone concentration

#two-way ANOVA on SA data
sa.model <- aov(sa_avg ~ microbe_level*diet_type, data = sa.sal)

#check assumptions
par(mfrow = c(2,2))
plot(sa.model)
  #also not great equal variance, and not great normality, but once again the small
  #number of points messes with these assumptions easily

summary(sa.model)
  #SIGNIFICANT effect of diet_type on phytohormone concentration. p<0.001

#two-way ANOVA on IAA data
iaa.model <- aov(iaa_avg ~microbe_level*diet_type, data = iaa.sal)

#check assumptions
par(mfrow(2,2))
plot(iaa.model)
  #variance is better, but normality is not optimal
summary(iaa.model)
  #we have significnace for ALL values!! The things that significantly affect
  #IAA concentrations are microbe_level, diet_type, and the interaction of the two
  #p=0.0040, p<0.001, p=0.0042 respectively

###now do post hoc analysis for these two-way ANOVAs
#BA post hoc
TukeyHSD(ba.model, conf.level=0.95) #somehow no significance anywhere. Like. What.
#plot(TukeyHSD(ba.model, conf.level=.95), las = 2)
t.test(ba_avg ~ diet_type, data = ba.sal, var.equal = TRUE) #significance, p=0.024

#SA post hoc
TukeyHSD(sa.model, conf.level=0.95)
  #all interactions between SA and control/BA diets are significant (p<0.05)

#IAA post hoc
TukeyHSD(iaa.model, conf.level = 0.95)
  #all interactions significant EXCEPT nonaxenic control and axenic control are the same (p=0.999)

#####unspiked data phytohormone comparisons####
unspiked.sal <- data.sal %>%
  filter(diet_type == "control")

##graphs of this data
unspiked.ba.point <- ggplot(unspiked.sal, aes(x = microbe_level, y =ba_avg))+
  geom_point(size=4, aes(color = microbe_level))+
  scale_color_manual(values = c("blue", "darkblue"))
unspiked.ba.point

unspiked.sa.point <- ggplot(unspiked.sal, aes(x = microbe_level, y =sa_avg))+
  geom_point(size=4, aes(color = microbe_level))+
  scale_color_manual(values = c("green", "darkgreen"))
unspiked.sa.point

unspiked.iaa.point <- ggplot(unspiked.sal, aes(x = microbe_level, y=iaa_avg))+
  geom_point(size = 4, aes(color = microbe_level))+
  scale_color_manual(values = c("pink", "maroon"))
unspiked.iaa.point

##running one-way ANOVAs, just to see
#BA
unspiked.ba.anova <- aov(ba_avg~microbe_level, data = unspiked.sal)
par(mfrow=c(2,2))
plot(unspiked.ba.anova)
summary(unspiked.ba.anova) #no significance, p=0.394
TukeyHSD(unspiked.ba.anova, conf.level=0.95) #no sig, p=0.394

unspiked.sa.anova <- aov(sa_avg~microbe_level, data = unspiked.sal)
plot(unspiked.sa.anova)
summary(unspiked.sa.anova) #no sig, p=0.201
TukeyHSD(unspiked.sa.anova, conf.level=0.95) #no sig, p=0.201

unspiked.iaa.anova <- aov(iaa_avg~microbe_level, data = unspiked.sal)
plot(unspiked.iaa.anova)
summary(unspiked.iaa.anova) #no sig, p=0.288
TukeyHSD(unspiked.iaa.anova, conf.level=0.95) #no sig, p=0.288
