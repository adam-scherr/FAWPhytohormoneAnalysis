### Making Graphics for the combined third and second sequestration experiment data
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
 # labs(x = "Diet Treatment", y = "BA Concentration (ng/g)") +
  #theme_clean() +
  # scale_x_discrete(labels = c("Control", "BA-Spiked")) +
  # theme(axis.text = element_text(size =15),
  #       panel.grid.major = element_line(size = 2),
  #       axis.title = element_text(size=20),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size=15))
  
ba_box

sa_box <-  ggplot(sa.only.data, aes(x = seq_experiment, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_color_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,45000))
 # labs(x = "Diet Treatment", y = "SA Concentration (ng/g)") 
  # theme_clean() +
  # scale_x_discrete(labels = c("Control", "SA-Spiked")) +
  # theme(axis.text = element_text(size =15),
  #       panel.grid.major = element_line(size = 2),
  #       axis.title = element_text(size=20),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size=15))
sa_box

iaa_box <-  ggplot(iaa.only.data, aes(x = seq_experiment, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_color_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,50000))
  # labs(x = "Diet Treatment", y = "IAA Concentration (ng/g)") +
  # theme_clean() +
  # scale_x_discrete(labels = c("Control", "IAA-Spiked")) +
  # theme(axis.text = element_text(size =15),
  #       panel.grid.major = element_line(size = 2),
  #       axis.title = element_text(size=20),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size=15))
iaa_box


#####Graphs of just the Unspiked Diet Data####
#first make the unspiked diet data sets
unspiked.ba.data <- ba.data %>%
  filter(diet_type == "control")

unspiked.sa.data <- sa.data %>%
  filter(diet_type == "control")

unspiked.iaa.data <- iaa.data %>%
  filter(diet_type == "control")

#now make the graphs
unspiked.ba.box <- ggplot(unspiked.ba.data, aes(x = seq_experiment, y = ba_avg, fill=microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"))+
  scale_color_manual(values = c("#9ecae1", "#2171b5", "#08306b"))
# labs(x = "Microbe Treatment", y = "BA Concentration (ng/g)") +
#   theme_clean() +
#  scale_x_discrete(labels = c("Axenic", "E. mundtii", "Xenic")) +
#   theme(axis.text = element_text(size =15),
#         panel.grid.major = element_line(size = 2),
#         axis.title = element_text(size=20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=15))
unspiked.ba.box

unspiked.sa.box <- ggplot(unspiked.sa.data, aes(x = seq_experiment, y = sa_avg, fill=microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_color_manual(values = c("#a1d99b", "#41ab5d", "#00441b")) 
  # labs(x = "Microbe Treatment", y = "SA Concentration (ng/g)") +
  # theme_clean() +
  # scale_x_discrete(labels = c("Axenic", "E. mundtii", "Xenic")) +
  # theme(axis.text = element_text(size =15),
  #       panel.grid.major = element_line(size = 2),
  #       axis.title = element_text(size=20),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size=15))
unspiked.sa.box

unspiked.iaa.box <- ggplot(unspiked.iaa.data, aes(x = seq_experiment, y = iaa_avg, fill=microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_color_manual(values = c("#fdae6b", "#f16913" ,"#7f2704")) 
  # labs(x = "Microbe Treatment", y = "IAA Concentration (ng/g)") +
  # theme_clean() +
  # scale_x_discrete(labels = c("Axenic", "E. mundtii", "Xenic")) +
  # theme(axis.text = element_text(size =15),
  #       panel.grid.major = element_line(size = 2),
  #       axis.title = element_text(size=20),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size=15))
unspiked.iaa.box

#####Combined Data Graphs####
#These graphs will show the combined 2nd and 3rd experiments' data as replicates of the same treatment
combined_ba_box <- ggplot(ba.only.data, aes(x = microbe_level, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, aes(color = seq_experiment))+
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

combined_ba_box

combined_sa_box <-  ggplot(sa.only.data, aes(x = microbe_level, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, aes(color = seq_experiment))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_color_manual(values = c("#C594D2", "#A01C1C"))+
 scale_y_continuous(expand = c(0,0), limits = c(0,45000))
combined_sa_box

combined_iaa_box <-  ggplot(iaa.only.data, aes(x = microbe_level, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, aes(color = seq_experiment))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_color_manual(values = c("#C594D2", "#A01C1C"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,50000))
combined_iaa_box
