###Analysis of BA, SA, and IAA Sequestration data
####Salivary Gland and Tissue Graphs of Phytohormone Levels, averages and standard errors only
###Adam Scherr, 11/27/2023
###Version 3, 08/28/2024

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

#now change the levels so graphs aren't alphabetically ordered
data.1$diet_type <- factor(data.1$diet_type, 
                                levels = c("control", "benzoic", "salicylic", "indole"))
data.1$tissue <- factor(data.1$tissue,
                        levels = c("sal", "fat", "hem"))
levels(data.1$diet_type) 
levels(data.1$tissue) #it worked, lets move one

#now make the negative values into zeros
data.1$ba_avg[data.1$ba_avg < 0] <- 0
data.1$sa_avg[data.1$sa_avg < 0] <- 0
data.1$iaa_avg[data.1$iaa_avg < 0] <- 0


#make a little preliminary box and whisker, lets see what she looks like
library(ggplot2)
library(ggthemes)


##Making datasets
#data set for salivary gland data only
data.sal <- data.1 %>%
  filter(tissue == "sal")

#data set for benzoic acid spiked diet only
ba.data <- data.sal %>%
  select(1:4, 9)%>%
  filter (diet_type == "benzoic")


# data set sa spiked diet only
sa.data <- data.sal %>%
  select(1:4, 12)%>%
  filter(diet_type == "salicylic")


#data set for iaa spiked diet only
iaa.data <- data.sal%>%
  select(1:4, 15)%>%
  filter(diet_type == "indole")

##subset the data for our averages and standard errors
#for BA
axenic.data.ba <- ba.data %>%
  filter(microbe_level == "axenic")

xenic.data.ba <- ba.data %>%
  filter(microbe_level == "nonaxenic")


#for SA
axenic.data.sa <- sa.data %>%
  filter(microbe_level == "axenic")

xenic.data.sa <- sa.data %>%
  filter(microbe_level == "nonaxenic")


#for IAA
axenic.data.iaa <- iaa.data %>%
  filter(microbe_level == "axenic")

xenic.data.iaa <- iaa.data %>%
  filter(microbe_level == "nonaxenic")

##now calculate the averages
#for BA
axenic_avg.ba <- mean(axenic.data.ba$ba_avg)
xenic_avg.ba <- mean(xenic.data.ba$ba_avg)

#for SA
axenic_avg.sa <- mean(axenic.data.sa$sa_avg)
xenic_avg.sa <- mean(xenic.data.sa$sa_avg)

#for IAA
axenic_avg.iaa <- mean(axenic.data.iaa$iaa_avg)
xenic_avg.iaa <- mean(xenic.data.iaa$iaa_avg)

##calculate standard error
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

#standard error for BA
axenic_sterr.ba <- st.error(axenic.data.ba$ba_avg)
xenic_sterr.ba <- st.error(xenic.data.ba$ba_avg)

#for SA
axenic_sterr.sa <- st.error(axenic.data.sa$sa_avg)
xenic_sterr.sa <- st.error(xenic.data.sa$sa_avg)

#for IAA
axenic_sterr.iaa <- st.error(axenic.data.iaa$iaa_avg)
xenic_sterr.iaa <- st.error(xenic.data.iaa$iaa_avg)

##lastly, making data frames combining all these data
microbe_treatment <- c("Axenic","Axenic","Axenic", "Xenic","Xenic","Xenic")
diet_treatment <- c("benzoic", "salicylic", "indole",
                    "benzoic", "salicylic", "indole")
spiked_avg <- c(axenic_avg.ba, axenic_avg.sa, axenic_avg.iaa,
                xenic_avg.ba,xenic_avg.sa, xenic_avg.iaa)
spiked_sterr <- c(axenic_sterr.ba,axenic_sterr.sa,axenic_sterr.iaa,
                  xenic_sterr.ba,xenic_sterr.sa,xenic_sterr.iaa)
df <- data.frame(microbe_treatment, diet_treatment, spiked_avg, spiked_sterr)

df.ba <- df %>%
  filter(diet_treatment == "benzoic")

df.sa <- df %>%
  filter(diet_treatment == "salicylic")

df.iaa <- df %>%
  filter(diet_treatment == "indole")

#####Making Graphs of spiked salivary glands using only spiked diet####
ba.bar <- ggplot(df.ba, aes(x = microbe_treatment, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8)+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.6, linewidth=1)+
  scale_fill_manual(values = c("#9ecae1", "#08306b"), name = "")+
  labs(x = "Microbe Treatment", y= "BA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
ba.bar

sa.bar <- ggplot(df.sa, aes(x = microbe_treatment, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8)+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.6, linewidth=1)+
  scale_fill_manual(values = c("#a1d99b",  "#00441b"), name = "") +
  labs(x = "Microbe Treatment", y= "SA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
sa.bar

iaa.bar <- ggplot(df.iaa, aes(x = microbe_treatment, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8)+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.6, linewidth=1)+
  scale_fill_manual(values = c("#fdae6b","#7f2704"), name = "")+
  labs(x = "Microbe Treatment", y= "IAA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
iaa.bar


#### tissue graphs ####
##Now make the 4 graphs comparing phytohormone levels among each tissue
#BA-fed FAW in different tissues
BA_tissue <- ggplot(data.ba, aes(x = tissue, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Xenic"))+
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
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Xenic"))+
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
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Xenic"))+
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
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Xenic"))+
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
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Xenic"))+
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
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Xenic"))+
  theme_clean() +
  labs(x = "Tissue", y = "IAA Concentration (ng/g)")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
control_tissue_iaa

####Salivary Gland Graphs####
#BA and BA-spiked diet dataset
BA_simple <- data.sal %>%
  filter(diet_type != "salicylic", diet_type != "indole") %>%
  select(1:9)

BA_simpler <- BA_simple %>%
  filter(diet_type == "benzoic")

#SA and SA-spiked diet dataset
SA_simple <- data.sal %>%
  filter(diet_type != "indole") %>%
  select(1:8, 12)

#IAA and IAA_spiked diet dataset
IAA_simple <- data.sal %>%
  filter(diet_type != "benzoic",  diet_type != "salicylic") %>%
  select(1:8, 15)

#Let's visualize these graphs
#graph for BA
BA_simple_sal <- ggplot(BA_simple, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "",
                     labels = c("Axenic", "Xenic"))+
  labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
  scale_y_log10(n.breaks = 6)+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_simple_sal

library(scales)
#graph for SA
SA_simple_sal <- ggplot (SA_simple, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "",
                    labels= c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "",
                     labels = c("Axenic", "Xenic"))+
  labs(x = "Diet Treatment", y = "SA Concentration (ng/g)")+
  scale_y_log10(n.breaks = 6)+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA-Spiked", "SA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_simple_sal

#graph for IAA
IAA_simple_sal <- ggplot (IAA_simple, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "",
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "",
                     labels = c("Axenic", "Xenic"))+
  labs(x = "Diet Treatment", y = "IAA Concentration (ng/g)")+
  scale_y_log10(n.breaks = 6)+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "IAA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_simple_sal

#attempts to make the breaks work, but the results were not ideal
scale_y_log10(breaks = seq(0, 11000, by = 100))
scale_y_log10(breaks = pretty(IAA_simple$iaa_avg, n = 100))
library(scales)
scale_y_log10(breaks = scales ::pretty_breaks(n=10)
              
