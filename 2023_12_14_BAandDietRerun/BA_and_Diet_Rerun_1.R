###BA Rerun and Diet Run R code
###Adam Scherr, 12/14/2023

#step 1, upload your data file
#first time making this, I uploaded it from the local computer data
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2023_12_14_BAandDietRerun/BA%26Diet%20Rerun_1%20-%20tidy_data.2.csv")

#simplify the dataset
library(dplyr)
data.1 <- data %>%
  filter(rep_avg > 1) %>%
  select(1:3, 18) %>%
  mutate_at(c("microbe_level", "color"), as.factor)
str(data.1)

#lets make a graph of this puppy
library(ggplot2)
library(ggthemes)            
BA_rerun <- ggplot(data.1, aes(x = microbe_level, y = rep_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"))+
  scale_color_manual(values = c("#74a9cf", "#023858"))+
  labs(x = "Microbe Level", y = "BA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_rerun

#what if we made the new and old points stand out more?
BA_rerun.2 <- ggplot(data.1, aes(x = microbe_level, y = rep_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= color)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"))+
  scale_color_manual(values = c("purple", "yellow"), labels = c("New Data", "Prelim Data"))+
  labs(x = "Microbe Level", y = "BA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_rerun.2

# BA_sal <- ggplot(data.sal, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
#   geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
#   geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
#   scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
#                     labels = c("Axenic", "Nonaxenic"))+
#   scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
#                      labels = c("Axenic", "Nonaxenic"))+
#   labs(x = "Diet Type", y = "BA Concentration (ng/g)")+
#   theme_clean() +
#   scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
#   theme(axis.text = element_text(size =15),
#         panel.grid.major = element_line(size = 2),
#         axis.title = element_text(size=20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=15))

BA_point <- ggplot(data.1, aes(x = microbe_level, y = rep_avg))+
  geom_point()
BA_point

data.2 <- data.1 %>%
  filter (microbe_level == "axenic")
BA_axenic_point <- ggplot (data.2, aes(x = microbe_level, y=rep_avg)) +
  geom_point()
BA_axenic_point

#for Kelli, one graph that shows the techincal replicate data only
data.2 <- data.1 %>%
  filter(color == "blue")

BA_rerun.1 <- ggplot(data.2, aes(x = microbe_level, y = rep_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"))+
  scale_color_manual(values = c("#74a9cf", "#023858"))+
  labs(x = "Microbe Level", y = "BA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_rerun.1

#two pairs of points from the two axenic runs are VERY close (within 50 ng/g of each other)
#but there is an outlier on the higher end, double the value of one of the pairs


##lets see if these are all statistically significant from each other

#start by seeing if it looks normal
plot(data.1$rep_avg) #eh, maybe normal
qqnorm(data.1$rep_avg) #eh, kinda a straight line (but also curved)
shapiro.test(data.1$rep_avg) #p=0.0098, NOT normally distributed

##lets check our three assumptions discussed in this resource:
#https://www.datanovia.com/en/lessons/t-test-in-r/

#are there outliers?
library(rstatix)
data.1 %>%
  identify_outliers(rep_avg) #it might be saying that blue NB1 is an outlier?

##is it normaly?
#no, it is not, we already tested this

#are variances equal?
library(car)
leveneTest(rep_avg ~ microbe_level, data = data.1)  
  #p = 0.046, this means that there is a sig difference in the variances and they
      #are NOT EQUAL. I repeat, variances are NOT equal. And we failed all over our other tests...

##our nonparametric test: The Mann-Whitney U Test (also called Wilcoxon Rank-Sum Test)
wilcox.test(rep_avg ~ microbe_level, data = data.1)
  #we're getting a p-value of 0.4127, which is NOT significant. Welp, there you have it

#lets try with the goal of seeing if axenic is less than nonaxenic helps much
wilcox.test(rep_avg ~ microbe_level, data = data.1, alternative = "less")
  #p = 0.206, still a good bit above 0.05

#despite all my best efforts, we're still getting a p-value well above 0.05
            
            