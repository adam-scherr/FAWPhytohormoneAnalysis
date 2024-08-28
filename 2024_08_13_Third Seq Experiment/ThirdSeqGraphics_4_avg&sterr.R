### Making Graphics for Third (Final) Sequestration Experiment, based on average and standard error
### Adam Scherr, August 28th, 2024
### Version 4

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

#make datasets for each phytohormone, just based on the spiked diet, NOT unspiked/control diet
ba.data <- data.1 %>%
  select(1:6) %>%
  filter(diet_type == "benzoic")

sa.data <- data.1 %>%
  select(1:3, 7:9) %>%
  filter(diet_type == "salicylic")

iaa.data <- data.1 %>%
  select(1:3, 10:12) %>%
  filter(diet_type=="indole")


#####Making graphs using averages and standard errors
###first, we must make the averages and standard errors data
##subset the data for our averages and standard errors
#for BA
axenic.data.ba <- ba.data %>%
  filter(microbe_level == "axenic")

xenic.data.ba <- ba.data %>%
  filter(microbe_level == "xenic")

emun.data.ba <- ba.data %>%
  filter(microbe_level == "E. mundtii")

#for SA
axenic.data.sa <- sa.data %>%
  filter(microbe_level == "axenic")

xenic.data.sa <- sa.data %>%
  filter(microbe_level == "xenic")

emun.data.sa <- sa.data %>%
  filter(microbe_level == "E. mundtii")

#for IAA
axenic.data.iaa <- iaa.data %>%
  filter(microbe_level == "axenic")

xenic.data.iaa <- iaa.data %>%
  filter(microbe_level == "xenic")

emun.data.iaa <- iaa.data %>%
  filter(microbe_level == "E. mundtii")

##now calculate the averages
#for BA
axenic_avg.ba <- mean(axenic.data.ba$ba_avg)
xenic_avg.ba <- mean(xenic.data.ba$ba_avg)
emun_avg.ba <- mean(emun.data.ba$ba_avg)

#for SA
axenic_avg.sa <- mean(axenic.data.sa$sa_avg)
xenic_avg.sa <- mean(xenic.data.sa$sa_avg)
emun_avg.sa <- mean(emun.data.sa$sa_avg)

#for IAA
axenic_avg.iaa <- mean(axenic.data.iaa$iaa_avg)
xenic_avg.iaa <- mean(xenic.data.iaa$iaa_avg)
emun_avg.iaa <- mean(emun.data.iaa$iaa_avg)


##calculate standard error
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

#standard error for BA
axenic_sterr.ba <- st.error(axenic.data.ba$ba_avg)
xenic_sterr.ba <- st.error(xenic.data.ba$ba_avg)
emun_sterr.ba <- st.error(emun.data.ba$ba_avg)

#for SA
axenic_sterr.sa <- st.error(axenic.data.sa$sa_avg)
xenic_sterr.sa <- st.error(xenic.data.sa$sa_avg)
emun_sterr.sa <- st.error(emun.data.sa$sa_avg)

#for IAA
axenic_sterr.iaa <- st.error(axenic.data.iaa$iaa_avg)
xenic_sterr.iaa <- st.error(xenic.data.iaa$iaa_avg)
emun_sterr.iaa <- st.error(emun.data.iaa$iaa_avg)

##lastly, making data frames combining all these data
microbe_treatment <- c("Axenic","Axenic","Axenic", "Xenic","Xenic","Xenic",
                       "E. mundtii","E. mundtii","E. mundtii")
diet_treatment <- c("benzoic", "salicylic", "indole","benzoic", "salicylic", "indole",
                    "benzoic", "salicylic", "indole")
spiked_avg <- c(axenic_avg.ba, axenic_avg.sa, axenic_avg.iaa,
                xenic_avg.ba,xenic_avg.sa, xenic_avg.iaa,
                emun_avg.ba,emun_avg.sa,emun_avg.iaa)
spiked_sterr <- c(axenic_sterr.ba,axenic_sterr.sa,axenic_sterr.iaa,
                  xenic_sterr.ba,xenic_sterr.sa,xenic_sterr.iaa,
                  emun_sterr.ba, emun_sterr.sa, emun_sterr.iaa)
df <- data.frame(microbe_treatment, diet_treatment, spiked_avg, spiked_sterr)

df.ba <- df %>%
  filter(diet_treatment == "benzoic")

df.sa <- df %>%
  filter(diet_treatment == "salicylic")

df.iaa <- df %>%
  filter(diet_treatment == "indole")

#####Making Graphs using averages and standard errors
ba.bar <- ggplot(df.ba, aes(x = microbe_treatment, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8)+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), width=0.6)+
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"), name = "") +
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
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), width=0.6)+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"), name = "") +
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
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), width=0.6)+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"), name = "")+
  labs(x = "Microbe Treatment", y= "IAA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
iaa.bar


