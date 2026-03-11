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

#####Making Graphs of spiked salivary glands using only spiked diet####
ba.bar <- ggplot(df.ba, aes(x = microbe_treatment, y = spiked_avg, fill = microbe_treatment))+
  geom_col(width=0.8)+
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.6, linewidth=1)+
  scale_fill_manual(values = c("#9ecae1", "#08306b"), name = "")+
  labs(x = "Microbe Treatment", y= "BA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =25),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=30),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
ba.bar


#BA and Unspiked diet dataset####
BA_simple <- data.sal %>%
  filter(diet_type != "salicylic", diet_type != "indole") %>%
  select(1:9)


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


#Now you need to get the averages and standard errors of each of these groups:
 #1 Control diet - axenic
 #2 Control diet - xenic
 #3 spiked diet - axenic
 #4 spiked diet - xenic

#make data frames for each grouping:
df.ax.ba <- BA_simple %>%
  filter(microbe_level == "axenic", diet_type == "benzoic")

df.x.ba <- BA_simple %>%
  filter(microbe_level == "nonaxenic", diet_type == "benzoic")

df.ax.con <- BA_simple %>%
  filter(microbe_level == "axenic", diet_type == "control")

df.x.con <- BA_simple %>%
  filter(microbe_level == "nonaxenic", diet_type == "control")

#make means from those data frames
ax.ba <- mean(df.ax.ba$ba_avg)
x.ba<- mean(df.x.ba$ba_avg)
ax.con <- mean(df.ax.con$ba_avg)
x.con<- mean(df.x.con$ba_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.ba.err <- st.error(df.ax.ba$ba_avg)
x.ba.err <- st.error(df.x.ba$ba_avg)
ax.con.err <- st.error(df.x.con$ba_avg)
x.con.err <- st.error(df.x.con$ba_avg)

#combine them all into one dataframe
microbe_treatment.1 <- c("Axenic", "Xenic", "Axenic", "Xenic")
diet_treatment.1 <- c("benzoic", "benzoic", "control", "control")
ba_avg <- c(ax.ba, x.ba, ax.con, x.con)
sterr <- c(ax.ba.err, x.ba.err, ax.con.err, x.con.err)

df.ba.control <- data.frame(microbe_treatment.1, diet_treatment.1, ba_avg, sterr)
View(df.ba.control)

#make levels so that it graphs as control diet first, then BA-spiked diet
data.1$tissue <- factor(data.1$tissue,
                        levels = c("sal", "fat", "hem"))
df.ba.control$diet_treatment.1 <- factor(df.ba.control$diet_treatment.1,
                                         levels = c("control", "benzoic"))
#final graph of the whole shebang
BA_simple_bar <- ggplot(df.ba.control, aes(x = diet_treatment.1, y = ba_avg, fill = microbe_treatment.1))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= ba_avg-sterr, ymax= ba_avg+sterr), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#9ecae1", "#08306b"), name = "")+
  scale_x_discrete(labels = c("Unspiked", "BA-Spiked"))+
 labs(x = "Diet Type", y= "BA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_simple_bar

##SA and Unspiked diet dataset####
SA_simple <- data.sal %>%
  filter(diet_type != "benzoic", diet_type != "indole")%>%
  select(1:8,12)

#Let's visualize these graphs
#graph for SA
SA_simple_sal <- ggplot(SA_simple, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "", 
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "",
                     labels = c("Axenic", "Xenic"))+
  labs(x = "Diet Treatment", y = "SA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "SA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_simple_sal

#Now you need to get the averages and standard errors of each of these groups:
#1 Control diet - axenic
#2 Control diet - xenic
#3 spiked diet - axenic
#4 spiked diet - xenic

#make data frames for each grouping:
df.ax.sa <- SA_simple %>%
  filter(microbe_level == "axenic", diet_type == "salicylic")

df.x.sa <- SA_simple %>%
  filter(microbe_level == "nonaxenic", diet_type == "salicylic")

df.ax.cons <- SA_simple %>%
  filter(microbe_level == "axenic", diet_type == "control")

df.x.cons <- SA_simple %>%
  filter(microbe_level == "nonaxenic", diet_type == "control")

#make means from those data frames
ax.sa <- mean(df.ax.sa$sa_avg)
x.sa<- mean(df.x.sa$sa_avg)
ax.cons <- mean(df.ax.cons$sa_avg)
x.cons<- mean(df.x.cons$sa_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.sa.err <- st.error(df.ax.sa$sa_avg)
x.sa.err <- st.error(df.x.sa$sa_avg)
ax.cons.err <- st.error(df.x.cons$sa_avg)
x.cons.err <- st.error(df.x.cons$sa_avg)

#combine them all into one dataframe
microbe_treatment.1s <- c("Axenic", "Xenic", "Axenic", "Xenic")
diet_treatment.1s <- c("salicylic", "salicylic", "control", "control")
sa_avg <- c(ax.sa, x.sa, ax.cons, x.cons)
sterrs <- c(ax.sa.err, x.sa.err, ax.cons.err, x.cons.err)

df.sa.control <- data.frame(microbe_treatment.1s, diet_treatment.1s, sa_avg, sterrs)
View(df.sa.control)

#final graph of the whole shebang
SA_simple_bar <- ggplot(df.sa.control, aes(x = diet_treatment.1s, y = sa_avg, fill = microbe_treatment.1s))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= sa_avg-sterrs, ymax= sa_avg+sterrs), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "")+
  scale_x_discrete(labels = c("Unspiked", "SA-Spiked"))+
  labs(x = "Diet Type", y= "SA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_simple_bar

#IAA and Unspiked diet dataset####
IAA_simple <- data.sal %>%
  filter(diet_type != "salicylic", diet_type != "benzoic") %>%
  select(1:8, 15)


#Let's visualize these graphs
#graph for BA
IAA_simple_sal <- ggplot(IAA_simple, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "", 
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "",
                     labels = c("Axenic", "Xenic"))+
  labs(x = "Diet Treatment", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "IAA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_simple_sal

#Now you need to get the averages and standard errors of each of these groups:
#1 Control diet - axenic
#2 Control diet - xenic
#3 spiked diet - axenic
#4 spiked diet - xenic

#make data frames for each grouping:
df.ax.iaa <- IAA_simple %>%
  filter(microbe_level == "axenic", diet_type == "indole")

df.x.iaa <- IAA_simple %>%
  filter(microbe_level == "nonaxenic", diet_type == "indole")

df.ax.coni <- IAA_simple %>%
  filter(microbe_level == "axenic", diet_type == "control")

df.x.coni <- IAA_simple %>%
  filter(microbe_level == "nonaxenic", diet_type == "control")

#make means from those data frames
ax.iaa <- mean(df.ax.iaa$iaa_avg)
x.iaa<- mean(df.x.iaa$iaa_avg)
ax.coni <- mean(df.ax.coni$iaa_avg)
x.coni<- mean(df.x.coni$iaa_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.iaa.err <- st.error(df.ax.iaa$iaa_avg)
x.iaa.err <- st.error(df.x.iaa$iaa_avg)
ax.coni.err <- st.error(df.x.coni$iaa_avg)
x.coni.err <- st.error(df.x.coni$iaa_avg)

#combine them all into one dataframe
microbe_treatment.1i <- c("Axenic", "Xenic", "Axenic", "Xenic")
diet_treatment.1i <- c("indole", "indole", "control", "control")
iaa_avg <- c(ax.iaa, x.iaa, ax.coni, x.coni)
sterri <- c(ax.iaa.err, x.iaa.err, ax.coni.err, x.coni.err)

df.iaa.control <- data.frame(microbe_treatment.1i, diet_treatment.1i, iaa_avg, sterri)
View(df.iaa.control)

#final graph of the whole shebang
IAA_simple_bar <- ggplot(df.iaa.control, aes(x = diet_treatment.1i, y = iaa_avg, fill = microbe_treatment.1i))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= iaa_avg-sterri, ymax= iaa_avg+sterri), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "")+
  scale_x_discrete(labels = c("Unspiked", "IAA-Spiked"))+
  labs(x = "Diet Type", y= "IAA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_simple_bar


              