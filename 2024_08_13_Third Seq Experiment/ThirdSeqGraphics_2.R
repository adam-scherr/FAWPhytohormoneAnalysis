### Making Graphics for Third (Final) Sequestration Experiment
### Adam Scherr, August 12th, 2024
### Version 2

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
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"), name = "",
                    labels = c("Axenic", "E. mundtii", "Xenic"))+
  scale_color_manual(values = c("#9ecae1", "#2171b5", "#08306b"), name ="",
                     labels = c("Axenic", "E. mundtii", "Xenic")) +
  labs(x = "Diet Treatment", y = "BA Concentration (ng/g)") +
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA-Spiked")) +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
  
ba_box

sa_box <-  ggplot(sa.data, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"), name = "",
                    labels = c("Axenic", "E. mundtii", "Xenic"))+
  scale_color_manual(values = c("#a1d99b", "#41ab5d", "#00441b"), name ="",
                     labels = c("Axenic", "E. mundtii", "Xenic"))+
  labs(x = "Diet Treatment", y = "SA Concentration (ng/g)") +
  theme_clean() +
  scale_x_discrete(labels = c("Control", "SA-Spiked")) +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
sa_box

iaa_box <-  ggplot(iaa.data, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"), name = "",
                    labels = c("Axenic", "E. mundtii", "Xenic"))+
  scale_color_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"), name ="",
                     labels = c("Axenic", "E. mundtii", "Xenic"))+
  labs(x = "Diet Treatment", y = "IAA Concentration (ng/g)") +
  theme_clean() +
  scale_x_discrete(labels = c("Control", "IAA-Spiked")) +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
iaa_box

# unspikedandba_box <- ggplot(unspikedandba.data, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
#   geom_boxplot(width=0.8, color = "black", alpha=0.82)+
#   geom_point(size = 4, alpha =1, position=position_dodge(width=0.8), aes(color=microbe_level))+
#   scale_fill_manual(values = c("#9ecae1",  "#08306b"))+
#   scale_color_manual(values =c("#9ecae1",  "#08306b"))+
#   labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
#   scale_y_log10(n.breaks = 6)+
#   theme_clean() +
#   scale_x_discrete(labels = c("Control", "BA-Spiked"))+
#   theme(axis.text = element_text(size =15),
#         panel.grid.major = element_line(size = 2),
#         axis.title = element_text(size=20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=15))

# ba_box <- ggplot(ba.data, aes(x = microbe_level, y=ba_avg, fill=microbe_level)) +
#   geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
#   geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
#   scale_fill_manual(values = c("#9ecae1", "#08306b"), name = "")+
#   scale_color_manual(values = c("#9ecae1", "#08306b"), name="") +
#   labs(x= "Microbe Treatment", y= "BA Concentration (ng/g)")+
#   theme_clean() +
#   theme(axis.text = element_text(size =15),
#         panel.grid.major = element_line(size = 2),
#         axis.title = element_text(size=20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=15))