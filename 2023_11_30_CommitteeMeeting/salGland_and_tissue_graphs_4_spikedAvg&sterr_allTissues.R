###Analysis of BA, SA, and IAA Sequestration data
####Tissue Graphs of Phytohormone Levels, averages and standard errors only
###Adam Scherr, 11/27/2023
###Version  4, 08/28/2024

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
#data set for benzoic acid spiked diet only
ba.data <- data.1 %>%
  select(1:4, 9)%>%
  filter (diet_type == "benzoic")


# data set sa spiked diet only
sa.data <- data.1 %>%
  select(1:4, 12)%>%
  filter(diet_type == "salicylic")


#data set for iaa spiked diet only
iaa.data <- data.1%>%
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

##subset the data by tissue
#all salivary glands
sal.ax.avg.ba <- axenic.data.ba %>%
  filter(tissue == "sal")
sal.ax.avg.sa <- axenic.data.sa%>%
  filter(tissue == "sal")
sal.ax.avg.iaa <- axenic.data.iaa%>%
  filter(tissue == "sal")
sal.x.avg.ba <- xenic.data.ba %>%
  filter(tissue == "sal")
sal.x.avg.sa <- xenic.data.sa%>%
  filter(tissue == "sal")
sal.x.avg.iaa <- xenic.data.iaa%>%
  filter(tissue == "sal")

#all fat body
fat.ax.avg.ba <- axenic.data.ba %>%
  filter(tissue == "fat")
fat.ax.avg.sa <- axenic.data.sa%>%
  filter(tissue == "fat")
fat.ax.avg.iaa <- axenic.data.iaa%>%
  filter(tissue == "fat")
fat.x.avg.ba <- xenic.data.ba %>%
  filter(tissue == "fat")
fat.x.avg.sa <- xenic.data.sa%>%
  filter(tissue == "fat")
fat.x.avg.iaa <- xenic.data.iaa%>%
  filter(tissue == "fat")

#all hemolymph
hem.ax.avg.ba <- axenic.data.ba %>%
  filter(tissue == "hem")
hem.ax.avg.sa <- axenic.data.sa%>%
  filter(tissue == "hem")
hem.ax.avg.iaa <- axenic.data.iaa%>%
  filter(tissue == "hem")
hem.x.avg.ba <- xenic.data.ba %>%
  filter(tissue == "hem")
hem.x.avg.sa <- xenic.data.sa%>%
  filter(tissue == "hem")
hem.x.avg.iaa <- xenic.data.iaa%>%
  filter(tissue == "hem")

##now calculate the averages
#for BA
sal.ax.ba <- mean(sal.ax.avg.ba$ba_avg)
sal.x.ba<- mean(sal.x.avg.ba$ba_avg)
fat.ax.ba <- mean(fat.ax.avg.ba$ba_avg)
fat.x.ba<- mean(fat.x.avg.ba$ba_avg)
hem.ax.ba <- mean(hem.ax.avg.ba$ba_avg)
hem.x.ba<- mean(hem.x.avg.ba$ba_avg)

#for SA
sal.ax.sa <- mean(sal.ax.avg.sa$sa_avg)
sal.x.sa <- mean(sal.x.avg.sa$sa_avg)
fat.ax.sa <- mean(fat.ax.avg.sa$sa_avg)
fat.x.sa <- mean(fat.x.avg.sa$sa_avg)
hem.ax.sa <- mean(hem.ax.avg.sa$sa_avg)
hem.x.sa <- mean(hem.x.avg.sa$sa_avg)

#for IAA
sal.ax.iaa <- mean(sal.ax.avg.iaa$iaa_avg)
sal.x.iaa <- mean(sal.x.avg.iaa$iaa_avg)
fat.ax.iaa <- mean(fat.ax.avg.iaa$iaa_avg)
fat.x.iaa <- mean(fat.x.avg.iaa$iaa_avg)
hem.ax.iaa <- mean(hem.ax.avg.iaa$iaa_avg)
hem.x.iaa <- mean(hem.x.avg.iaa$iaa_avg)

##calculate standard error
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

#standard error for BA
sal.ax.ba.err <- st.error(sal.ax.avg.ba$ba_avg)
sal.x.ba.err <- st.error(sal.x.avg.ba$ba_avg)
fat.ax.ba.err <- st.error(fat.ax.avg.ba$ba_avg)
fat.x.ba.err <- st.error(fat.x.avg.ba$ba_avg)
hem.ax.ba.err <- st.error(hem.ax.avg.ba$ba_avg)
hem.x.ba.err <- st.error(hem.x.avg.ba$ba_avg)
#for SA
sal.ax.sa.err <- st.error(sal.ax.avg.sa$sa_avg)
sal.x.sa.err <- st.error(sal.x.avg.sa$sa_avg)
fat.ax.sa.err <- st.error(fat.ax.avg.sa$sa_avg)
fat.x.sa.err <- st.error(fat.x.avg.sa$sa_avg)
hem.ax.sa.err <- st.error(hem.ax.avg.sa$sa_avg)
hem.x.sa.err <- st.error(hem.x.avg.sa$sa_avg)

#for IAA
sal.ax.iaa.err <- st.error(sal.ax.avg.iaa$iaa_avg)
sal.x.iaa.err <- st.error(sal.x.avg.iaa$iaa_avg)
fat.ax.iaa.err <- st.error(fat.ax.avg.iaa$iaa_avg)
fat.x.iaa.err <- st.error(fat.x.avg.iaa$iaa_avg)
hem.ax.iaa.err <- st.error(hem.ax.avg.iaa$iaa_avg)
hem.x.iaa.err <- st.error(hem.x.avg.iaa$iaa_avg)

##lastly, making data frames combining all these data
microbe_treatment <- c("Axenic","Axenic","Axenic", "Xenic","Xenic","Xenic",
                       "Axenic","Axenic","Axenic", "Xenic","Xenic","Xenic",
                       "Axenic","Axenic","Axenic", "Xenic","Xenic","Xenic")
diet_treatment <- c("benzoic", "salicylic", "indole", "benzoic", "salicylic", "indole",
                    "benzoic", "salicylic", "indole", "benzoic", "salicylic", "indole",
                    "benzoic", "salicylic", "indole", "benzoic", "salicylic", "indole")
tissue <- c("sal","sal","sal","sal","sal","sal",
            "fat","fat","fat","fat","fat","fat",
            "hem","hem","hem","hem","hem","hem")
spiked_avg <- c(sal.ax.ba, sal.ax.sa, sal.ax.iaa, sal.x.ba, sal.x.sa, sal.x.iaa,
                fat.ax.ba, fat.ax.sa, fat.ax.iaa, fat.x.ba, fat.x.sa, fat.x.iaa,
                hem.ax.ba, hem.ax.sa, hem.ax.iaa, hem.x.ba, hem.x.sa, hem.x.iaa)
spiked_sterr <- c(sal.ax.ba.err, sal.ax.sa.err, sal.ax.iaa.err, sal.x.ba.err, sal.x.sa.err, sal.x.iaa.err,
                  fat.ax.ba.err, fat.ax.sa.err, fat.ax.iaa.err, fat.x.ba.err, fat.x.sa.err, fat.x.iaa.err,
                  hem.ax.ba.err, hem.ax.sa.err, hem.ax.iaa.err, hem.x.ba.err, hem.x.sa.err, hem.x.iaa.err)
df <- data.frame(microbe_treatment, diet_treatment, tissue, spiked_avg, spiked_sterr)

df$tissue <- factor(df$tissue,
                        levels = c("sal", "fat", "hem"))

df.ba <- df %>%
  filter(diet_treatment == "benzoic")

df.sa <- df %>%
  filter(diet_treatment == "salicylic")

df.iaa <- df %>%
  filter(diet_treatment == "indole")

#####Making Graphs of spiked tissues using only spiked diet####
ba.bar <- ggplot(df.ba, aes(x = tissue, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#9ecae1", "#08306b"), name = "")+
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  labs(x = "Tissue", y= "BA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
ba.bar

sa.bar <- ggplot(df.sa, aes(x = tissue, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#a1d99b",  "#00441b"), name = "") +
  labs(x = "Tissue", y= "SA Concentration (ng/g)") +
    scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
sa.bar

iaa.bar <- ggplot(df.iaa, aes(x = tissue, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#fdae6b","#7f2704"), name = "")+
  labs(x = "Tissue", y= "IAA Concentration (ng/g)") +
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
iaa.bar

