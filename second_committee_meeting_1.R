###Analysis of BA, SA, and IAA Sequestration data, controled study
####For Committee Meeting and possibly general publication 
###Adam Scherr, 11/27/2023

#upload the dataset
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/EntSoc_SeqData_2%20-%20tidy_data.csv")

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
View(data.1)

#now change the levels so graphs aren't alphabetically ordered
data.1$diet_type <- factor(data.1$diet_type, 
                                levels = c("control", "benzoic", "salicylic", "indole"))
data.1$tissue <- factor(data.1$tissue,
                        levels = c("sal", "fat", "hem"))
levels(data.1$diet_type) #it worked, lets move one

#make a little preliminary box and whisker, lets see what she looks like
library(ggplot2)
library(ggthemes)
BA_whisker <- ggplot(data.1, aes(x = microbe_level, y = ba_avg,
                                 color = diet_type))+
  geom_boxplot()+
  geom_jitter(size = 3, width = 0.1)
BA_whisker
#this isn't so pretty, we'll make them look pretty in the end

#let's make a bunch of data sets to use for our graphs
#1. data set for salivary gland data only
data.sal <- data.1 %>%
  filter(tissue == "sal")
View(data.sal)

#2. data set for fat body
data.fat <- data.1 %>%
  filter(tissue == "fat")
View(data.fat)

#3. data set for hemolymph
data.hem <- data.1 %>%
  filter(tissue == "hem")
View(data.hem)

#4. data set for benzoic acid spiked diet only
data.ba <- data.1 %>%
  filter (diet_type == "benzoic")
View(data.ba)

#5 data set sa spiked diet only
data.sa <- data.1 %>%
  filter(diet_type == "salicylic")
View(data.sa)

#6 data set for iaa spiked diet only
data.iaa <- data.1 %>%
  filter(diet_type == "indole")
View(data.iaa)

#7. data set for control diet only
data.control <- data.1 %>%
  filter(diet_type == "control")
View(data.control)

##LET THE GRAPHING BEGIN!!
##firstly, the 9 graphs showing phytohormone levels in each tissue

#BA in salivary glands
BA_sal <- ggplot(data.sal, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "BA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_sal

#SA in salivary glands
SA_sal <- ggplot(data.sal, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "SA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_sal

#IAA in salivary glands
IAA_sal <- ggplot(data.sal, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_sal

#BA in fat body
BA_fat <- ggplot(data.fat, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "BA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_fat

#SA in fat body
SA_fat <- ggplot(data.fat, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "SA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_fat

#IAA in fat body
IAA_fat <- ggplot(data.fat, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_fat

#BA in hemolymph
BA_hem <- ggplot(data.hem, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "BA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_hem

#SA in hemolymph
SA_hem <- ggplot(data.hem, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "SA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_hem

#IAA in hemolymph
IAA_hem <- ggplot(data.hem, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_hem

##Now make the 4 graphs comparing phytohormone levels among each tissue
#BA-fed FAW in different tissues
BA_tissue <- ggplot(data.ba, aes(x = tissue, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  theme_clean() +
  labs(x = "Tissue", y = "BA Concentration (ng/g)")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_tissue

#SA-fed FAW in different tissues
SA_tissue <- ggplot(data.sa, aes(x = tissue, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  theme_clean() +
  labs(x = "Tissue", y = "SA Concentration (ng/g)")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_tissue
#IAA-fed FAW in different tissues
IAA_tissue <- ggplot(data.iaa, aes(x = tissue, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  theme_clean() +
  labs(x = "Tissue", y = "IAA Concentration (ng/g)")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_tissue

#control-fed FAW in different tissues
control_tissue_ba <- ggplot(data.control, aes(x = tissue, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  theme_clean() +
  labs(x = "Tissue", y = "BA Concentration (ng/g)")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
control_tissue_ba

control_tissue_sa <- ggplot(data.control, aes(x = tissue, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  theme_clean() +
  labs(x = "Tissue", y = "SA Concentration (ng/g)")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
control_tissue_sa

control_tissue_iaa <- ggplot(data.control, aes(x = tissue, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  theme_clean() +
  labs(x = "Tissue", y = "IAA Concentration (ng/g)")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
control_tissue_iaa
