### Making Graphics for Third Sequestration Experiment, combined Control and Spiked Avgs and errors
### Adam Scherr, March 31, 2026
### Version 1

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
  filter(diet_type == "benzoic" | diet_type == "control")

sa.data <- data.1 %>%
  select(1:3, 7:9) %>%
  filter(diet_type == "salicylic"| diet_type == "control")

iaa.data <- data.1 %>%
  select(1:3, 10:12) %>%
  filter(diet_type=="indole"| diet_type == "control")

####BA and Control####
#Let's visualize these graphs
#graph for BA
BA_simple_sal <- ggplot(ba.data, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"))+
  scale_color_manual(values = c("#9ecae1", "#2171b5", "#08306b"))+
  #labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
 # scale_y_log10(n.breaks = 6)+
  theme_clean() +
  #scale_x_discrete(labels = c("Control", "BA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_simple_sal


#Now you need to get the averages and standard errors of each of these groups:
#1 Control diet - axenic
#2 Control diet - E. mundtii
#3 Control diet - xenic
#4 spiked diet - axenic
#5 spiked diet - E. mundti
#6 spiked diet - xenic

#make data frames for each grouping:
df.ax.ba <- ba.data %>%
  filter(microbe_level == "axenic", diet_type == "benzoic")

df.e.ba <- ba.data %>%
  filter(microbe_level == "E. mundtii", diet_type == "benzoic")

df.x.ba <- ba.data %>%
  filter(microbe_level == "xenic", diet_type == "benzoic")

df.ax.con <- ba.data %>%
  filter(microbe_level == "axenic", diet_type == "control")

df.e.con <- ba.data %>%
  filter(microbe_level == "E. mundtii", diet_type == "control")

df.x.con <- ba.data %>%
  filter(microbe_level == "xenic", diet_type == "control")

#make means from those data frames
ax.ba <- mean(df.ax.ba$ba_avg)
x.ba<- mean(df.x.ba$ba_avg)
e.ba <- mean(df.e.ba$ba_avg)
ax.con <- mean(df.ax.con$ba_avg)
x.con<- mean(df.x.con$ba_avg)
e.con <- mean(df.e.con$ba_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.ba.err <- st.error(df.ax.ba$ba_avg)
x.ba.err <- st.error(df.x.ba$ba_avg)
e.ba.err <- st.error(df.e.ba$ba_avg)
ax.con.err <- st.error(df.ax.con$ba_avg)
x.con.err <- st.error(df.x.con$ba_avg)
e.con.err <- st.error(df.e.con$ba_avg)

#combine them all into one dataframe
microbe_treatment.b <- c("Axenic", "Xenic", "E. mundtii", "Axenic", "Xenic", "E. mundtii")
diet_treatment.b <- c("benzoic", "benzoic", "benzoic", "control", "control", "control")
ba_avg <- c(ax.ba, x.ba, e.ba, ax.con, x.con, e.con)
sterr <- c(ax.ba.err, x.ba.err, e.ba.err, ax.con.err, x.con.err, e.con.err)

df.ba.control <- data.frame(microbe_treatment.b, diet_treatment.b, ba_avg, sterr)
View(df.ba.control)

#make levels so that it graphs as control diet first, then BA-spiked diet
df.ba.control$diet_treatment.b <- factor(df.ba.control$diet_treatment.b,
                                         levels = c("control", "benzoic"))
#final graph of the whole shebang
BA_simple_bar <- ggplot(df.ba.control, aes(x = diet_treatment.b, y = ba_avg, fill = microbe_treatment.b))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= ba_avg-sterr, ymax= ba_avg+sterr), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#9ecae1", "#2171b5", "#08306b"), name = "")+
  scale_x_discrete(labels = c("Unspiked", "BA-Spiked"))+
  labs(x = "Diet Type", y= "BA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_simple_bar

####SA and Control####
#Let's visualize these graphs
#graph for SA
SA_simple_sal <- ggplot(sa.data, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  scale_color_manual(values = c("#a1d99b", "#41ab5d", "#00441b"))+
  #labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
  # scale_y_log10(n.breaks = 6)+
  theme_clean() +
  #scale_x_discrete(labels = c("Control", "BA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_simple_sal


#Now you need to get the averages and standard errors of each of these groups:
#1 Control diet - axenic
#2 Control diet - E. mundtii
#3 Control diet - xenic
#4 spiked diet - axenic
#5 spiked diet - E. mundti
#6 spiked diet - xenic

#make data frames for each grouping:
df.ax.sa <- sa.data %>%
  filter(microbe_level == "axenic", diet_type == "salicylic")

df.e.sa <- sa.data %>%
  filter(microbe_level == "E. mundtii", diet_type == "salicylic")

df.x.sa <- sa.data %>%
  filter(microbe_level == "xenic", diet_type == "salicylic")

df.ax.cons <- sa.data %>%
  filter(microbe_level == "axenic", diet_type == "control")

df.e.cons <- sa.data %>%
  filter(microbe_level == "E. mundtii", diet_type == "control")

df.x.cons <- sa.data %>%
  filter(microbe_level == "xenic", diet_type == "control")

#make means from those data frames
ax.sa <- mean(df.ax.sa$sa_avg)
x.sa<- mean(df.x.sa$sa_avg)
e.sa <- mean(df.e.sa$sa_avg)
ax.cons <- mean(df.ax.cons$sa_avg)
x.cons <- mean(df.x.cons$sa_avg)
e.cons <- mean(df.e.cons$sa_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.sa.err <- st.error(df.ax.sa$sa_avg)
x.sa.err <- st.error(df.x.sa$sa_avg)
e.sa.err <- st.error(df.e.sa$sa_avg)
ax.con.errs <- st.error(df.ax.cons$sa_avg)
x.con.errs <- st.error(df.x.cons$sa_avg)
e.con.errs <- st.error(df.e.cons$sa_avg)

#combine them all into one dataframe
microbe_treatment.s <- c("Axenic", "Xenic", "E. mundtii", "Axenic", "Xenic", "E. mundtii")
diet_treatment.s <- c("salicylic", "salicylic", "salicylic", "control", "control", "control")
sa_avg <- c(ax.sa, x.sa, e.sa, ax.cons, x.cons, e.cons)
sterr.s <- c(ax.sa.err, x.sa.err, e.sa.err, ax.con.errs, x.con.errs, e.con.errs)

df.sa.control <- data.frame(microbe_treatment.s, diet_treatment.s, sa_avg, sterr.s)
View(df.sa.control)

#make levels so that it graphs as control diet first, then BA-spiked diet
df.sa.control$diet_treatment.s <- factor(df.sa.control$diet_treatment.s,
                                         levels = c("control", "salicylic"))
#final graph of the whole shebang
SA_simple_bar <- ggplot(df.sa.control, aes(x = diet_treatment.s, y = sa_avg, fill = microbe_treatment.s))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= sa_avg-sterr.s, ymax= sa_avg+sterr.s), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#a1d99b", "#41ab5d", "#00441b"), name = "")+
  scale_x_discrete(labels = c("Unspiked", "SA-Spiked"))+
  labs(x = "Diet Type", y= "SA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_simple_bar

####IAA and Control####
#Let's visualize these graphs
#graph for SA
IAA_simple_sal <- ggplot(iaa.data, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  scale_color_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"))+
  #labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
  # scale_y_log10(n.breaks = 6)+
  theme_clean() +
  #scale_x_discrete(labels = c("Control", "BA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_simple_sal


#Now you need to get the averages and standard errors of each of these groups:
#1 Control diet - axenic
#2 Control diet - E. mundtii
#3 Control diet - xenic
#4 spiked diet - axenic
#5 spiked diet - E. mundti
#6 spiked diet - xenic

#make data frames for each grouping:
df.ax.iaa <- iaa.data %>%
  filter(microbe_level == "axenic", diet_type == "indole")

df.e.iaa <- iaa.data %>%
  filter(microbe_level == "E. mundtii", diet_type == "indole")

df.x.iaa <- iaa.data %>%
  filter(microbe_level == "xenic", diet_type == "indole")

df.ax.coni <- iaa.data %>%
  filter(microbe_level == "axenic", diet_type == "control")

df.e.coni <- iaa.data %>%
  filter(microbe_level == "E. mundtii", diet_type == "control")

df.x.coni <- iaa.data %>%
  filter(microbe_level == "xenic", diet_type == "control")

#make means from those data frames
ax.iaa <- mean(df.ax.iaa$iaa_avg)
x.iaa<- mean(df.x.iaa$iaa_avg)
e.iaa <- mean(df.e.iaa$iaa_avg)
ax.coni <- mean(df.ax.coni$iaa_avg)
x.coni <- mean(df.x.coni$iaa_avg)
e.coni <- mean(df.e.coni$iaa_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.iaa.err <- st.error(df.ax.iaa$iaa_avg)
x.iaa.err <- st.error(df.x.iaa$iaa_avg)
e.iaa.err <- st.error(df.e.iaa$iaa_avg)
ax.con.erri <- st.error(df.ax.coni$iaa_avg)
x.con.erri <- st.error(df.x.coni$iaa_avg)
e.con.erri <- st.error(df.e.coni$iaa_avg)

#combine them all into one dataframe
microbe_treatment.i <- c("Axenic", "Xenic", "E. mundtii", "Axenic", "Xenic", "E. mundtii")
diet_treatment.i <- c("indole", "indole", "indole", "control", "control", "control")
iaa_avg <- c(ax.iaa, x.iaa, e.iaa, ax.coni, x.coni, e.coni)
sterr.i <- c(ax.iaa.err, x.iaa.err, e.iaa.err, ax.con.erri, x.con.erri, e.con.erri)

df.iaa.control <- data.frame(microbe_treatment.i, diet_treatment.i, iaa_avg, sterr.i)
View(df.iaa.control)

#make levels so that it graphs as control diet first, then BA-spiked diet
df.iaa.control$diet_treatment.i <- factor(df.iaa.control$diet_treatment.i,
                                         levels = c("control", "indole"))
#final graph of the whole shebang
IAA_simple_bar <- ggplot(df.iaa.control, aes(x = diet_treatment.i, y = iaa_avg, fill = microbe_treatment.i))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= iaa_avg-sterr.i, ymax= iaa_avg+sterr.i), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"), name = "")+
  scale_x_discrete(labels = c("Unspiked", "IAA-Spiked"))+
  labs(x = "Diet Type", y= "IAA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_simple_bar


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
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), width=0.6,
                linewidth=1)+
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
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.6, linewidth=1)+
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
  geom_errorbar(aes(ymin= spiked_avg-spiked_sterr, ymax= spiked_avg+spiked_sterr), 
                width=0.6, linewidth=1)+
  scale_fill_manual(values = c("#fdae6b", "#f16913" ,"#7f2704"), name = "")+
  labs(x = "Microbe Treatment", y= "IAA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
iaa.bar


