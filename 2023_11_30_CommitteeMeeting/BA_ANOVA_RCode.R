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
  filter(microbe_level != "blank", nat_lab != "lab")
str(data.1)

data.1$diet_type <- factor(data.1$diet_type, 
                           levels = c("control", "benzoic", "salicylic", "indole"))
data.1$tissue <- factor(data.1$tissue,
                        levels = c("sal", "fat", "hem"))

data.2 <- data.1 %>%
  select(1:4, 9, 12, 15)

data.sal <- data.2 %>%
  filter(tissue == "sal")

##set up datasets for each phytohormone-control diet comparison
#BA and control
ba.sal <- data.sal %>%
  filter(diet_type != "indole") %>%
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

##using BA as first two-way ANOVA example
ba.model <- aov(ba_avg ~ microbe_level*diet_type, data = ba.sal)
summary(ba.model)


# ba.sal$scaled.ba_avg <- as.numeric(scale(ba.sal$ba_avg)) ## scaling response (subtract by mean, and diving by SD)
# ba.model.scale <- aov(exp(scaled.ba_avg) ~ microbe_level*diet_type, data = ba.sal)
# summary(ba.model.scale)
# 
# ba.model.lm <- lm(ba_avg ~ microbe_level*diet_type, data = ba.sal)
# summary(ba.model.lm)

#check our assumptions
par(mfrow = c(2, 2))
plot(ba.model.scale)

#QQ-plot and variance don't look great
#We have so few data points that it's easy to skew the normality and variance
#So let's go through with it since our small number of points is skewing normality and variance so much
summary(ba.model)
  #p=0.0261 for diet_type. Only diet_type has SIGNIFICANT effect on phytohormone concentration

###now do post hoc analysis for these two-way ANOVAs
#BA post hoc
TukeyHSD(ba.model, conf.level=0.95) #somehow no significance anywhere. Like. What.

#t-test shows that, if you group by diet_type alone (ignore microbe_level), you get significance
t.test(ba_avg ~ diet_type, data = ba.sal, var.equal = TRUE) #significance, p=0.024

