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
#data set for control diet only
c.data <- data.1 %>%
  select(1:4, 9, 12, 15)%>%
  filter (diet_type == "control")


##subset the data for our averages and standard errors
#by microbe treatment
axenic.data <- c.data %>%
  filter(microbe_level == "axenic")

xenic.data <- c.data %>%
  filter(microbe_level == "nonaxenic")


##subset the data by tissue
#all salivary glands
sal.ax <- axenic.data%>%
  filter(tissue == "sal")

sal.x <- axenic.data%>%
  filter(tissue == "sal")


#all fat body
fat.ax <- axenic.data%>%
  filter(tissue == "fat")

fat.x <- axenic.data%>%
  filter(tissue == "fat")

#all hemolymph
hem.ax <- axenic.data%>%
  filter(tissue == "hem")

hem.x <- axenic.data%>%
  filter(tissue == "hem")

##now calculate the averages
#for BA
sal.ax.ba <- mean(sal.ax$ba_avg)
sal.x.ba<- mean(sal.x$ba_avg)
fat.ax.ba <- mean(fat.ax$ba_avg)
fat.x.ba<- mean(fat.x$ba_avg)
hem.ax.ba <- mean(hem.ax$ba_avg)
hem.x.ba<- mean(hem.x$ba_avg)

#for SA
sal.ax.sa <- mean(sal.ax$sa_avg)
sal.x.sa <- mean(sal.x$sa_avg)
fat.ax.sa <- mean(fat.ax$sa_avg)
fat.x.sa <- mean(fat.x$sa_avg)
hem.ax.sa <- mean(hem.ax$sa_avg)
hem.x.sa <- mean(hem.x$sa_avg)

#for IAA
sal.ax.iaa <- mean(sal.ax$iaa_avg)
sal.x.iaa <- mean(sal.x$iaa_avg)
fat.ax.iaa <- mean(fat.ax$iaa_avg)
fat.x.iaa <- mean(fat.x$iaa_avg)
hem.ax.iaa <- mean(hem.ax$iaa_avg)
hem.x.iaa <- mean(hem.x$iaa_avg)

##calculate standard error
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

#standard error for BA
sal.ax.ba.err <- st.error(sal.ax$ba_avg)
sal.x.ba.err <- st.error(sal.x$ba_avg)
fat.ax.ba.err <- st.error(fat.ax$ba_avg)
fat.x.ba.err <- st.error(fat.x$ba_avg)
hem.ax.ba.err <- st.error(hem.ax$ba_avg)
hem.x.ba.err <- st.error(hem.x$ba_avg)
#for SA
sal.ax.sa.err <- st.error(sal.ax$sa_avg)
sal.x.sa.err <- st.error(sal.x$sa_avg)
fat.ax.sa.err <- st.error(fat.ax$sa_avg)
fat.x.sa.err <- st.error(fat.x$sa_avg)
hem.ax.sa.err <- st.error(hem.ax$sa_avg)
hem.x.sa.err <- st.error(hem.x$sa_avg)

#for IAA
sal.ax.iaa.err <- st.error(sal.ax$iaa_avg)
sal.x.iaa.err <- st.error(sal.x$iaa_avg)
fat.ax.iaa.err <- st.error(fat.ax$iaa_avg)
fat.x.iaa.err <- st.error(fat.x$iaa_avg)
hem.ax.iaa.err <- st.error(hem.ax$iaa_avg)
hem.x.iaa.err <- st.error(hem.x$iaa_avg)

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
  #this dataframe is what will go into the publication

df.ba <- df %>%
  filter(diet_treatment == "benzoic")

df.sa <- df %>%
  filter(diet_treatment == "salicylic")

df.iaa <- df %>%
  filter(diet_treatment == "indole")
