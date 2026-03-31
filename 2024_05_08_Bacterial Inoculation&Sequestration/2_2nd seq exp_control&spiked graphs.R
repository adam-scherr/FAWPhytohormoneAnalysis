### Analysis of phytohormone contents in salivary glands of 2nd sequestration experiment
### Adam Scherr, March 15, 2026
### Version 2, making avg and st. error graphs for control and spiked diet


#upload your datasets
data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_05_08_Bacterial%20Inoculation%26Sequestration/Inoc_Seq_Data.csv")

#load in packages
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(ggthemes)

####setting up datasets####
#change the object types so the right things are characters/factors
str(data)
data.1 <- data %>%
  select(!3) %>%
  mutate_at(c("microbe_level", "diet_type"), as.factor) %>%
  filter(microbe_level != "Frass Com.")

data.1$diet_type <- factor(data.1$diet_type, 
                           levels = c("control", "benzoic", "salicylic", "indole"))
str(data.1)

#now make the negative values into zeros
data.1$ba_avg[data.1$ba_avg < 0] <- 0
data.1$sa_avg[data.1$sa_avg < 0] <- 0
data.1$iaa_avg[data.1$iaa_avg < 0] <- 0


###BA and Control####
BA_simple <- data.1 %>%
  filter(diet_type != "salicylic", diet_type != "indole") %>%
  select(1:4)


#Let's visualize these graphs
#graph for BA
BA_simple_sal <- ggplot(BA_simple, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"))+
  scale_color_manual(values = c("#74a9cf", "#023858"))+
  labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
  scale_y_log10(n.breaks = 6)+
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
#3 spiked diet - axenic
#4 spiked diet - E. mundti

#make data frames for each grouping:
df.ax.ba <- BA_simple %>%
  filter(microbe_level == "Axenic", diet_type == "benzoic")

df.x.ba <- BA_simple %>%
  filter(microbe_level == "E. mundtii", diet_type == "benzoic")

df.ax.con <- BA_simple %>%
  filter(microbe_level == "Axenic", diet_type == "control")

df.x.con <- BA_simple %>%
  filter(microbe_level == "E. mundtii", diet_type == "control")

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
ax.con.err <- st.error(df.ax.con$ba_avg)
x.con.err <- st.error(df.x.con$ba_avg)

#combine them all into one dataframe
microbe_treatment.b <- c("Axenic", "E. mundtii", "Axenic", "E. mundtii")
diet_treatment.b <- c("benzoic", "benzoic", "control", "control")
ba_avg <- c(ax.ba, x.ba, ax.con, x.con)
sterr <- c(ax.ba.err, x.ba.err, ax.con.err, x.con.err)

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


####SA and Control####
SA_simple <- data.1 %>%
  filter(diet_type != "benzoic", diet_type != "indole") %>%
  select(1:3, 7)


#Let's visualize these graphs
#graph for SA
SA_simple_sal <- ggplot(SA_simple, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"))+
 # labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
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
#3 spiked diet - axenic
#4 spiked diet - E. mundti

#make data frames for each grouping:
df.ax.sa <- SA_simple %>%
  filter(microbe_level == "Axenic", diet_type == "salicylic")

df.x.sa <- SA_simple %>%
  filter(microbe_level == "E. mundtii", diet_type == "salicylic")

df.ax.cons <- SA_simple %>%
  filter(microbe_level == "Axenic", diet_type == "control")

df.x.cons <- SA_simple %>%
  filter(microbe_level == "E. mundtii", diet_type == "control")

#make means from those data frames
ax.sa <- mean(df.ax.sa$sa_avg)
x.sa<- mean(df.x.sa$sa_avg)
ax.cons <- mean(df.ax.cons$sa_avg)
x.cons <- mean(df.x.cons$sa_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.sa.err <- st.error(df.ax.sa$sa_avg)
x.sa.err <- st.error(df.x.sa$sa_avg)
ax.con.errs <- st.error(df.ax.cons$sa_avg)
x.con.errs <- st.error(df.x.cons$sa_avg)

#combine them all into one dataframe
microbe_treatment.s <- c("Axenic", "E. mundtii", "Axenic", "E. mundtii")
diet_treatment.s <- c("salicylic", "salicylic", "control", "control")
sa_avg <- c(ax.sa, x.sa, ax.cons, x.cons)
sterr <- c(ax.sa.err, x.sa.err, ax.con.errs, x.con.errs)

df.sa.control <- data.frame(microbe_treatment.s, diet_treatment.s, sa_avg, sterr)
View(df.sa.control)

#make levels so that it graphs as control diet first, then BA-spiked diet
df.sa.control$diet_treatment.s <- factor(df.sa.control$diet_treatment.s,
                                         levels = c("control", "salicylic"))
#final graph of the whole shebang
SA_simple_bar <- ggplot(df.sa.control, aes(x = diet_treatment.s, y = sa_avg, fill = microbe_treatment.s))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= sa_avg-sterr, ymax= sa_avg+sterr), 
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

####IAA and Control####
IAA_simple <- data.1 %>%
  filter(diet_type != "benzoic", diet_type != "salicylic") %>%
  select(1:3, 10)


#Let's visualize these graphs
#graph for SA
IAA_simple_sal <- ggplot(IAA_simple, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"))+
  # labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
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
#3 spiked diet - axenic
#4 spiked diet - E. mundti

#make data frames for each grouping:
df.ax.iaa <- IAA_simple %>%
  filter(microbe_level == "Axenic", diet_type == "indole")

df.x.iaa <- IAA_simple %>%
  filter(microbe_level == "E. mundtii", diet_type == "indole")

df.ax.coni <- IAA_simple %>%
  filter(microbe_level == "Axenic", diet_type == "control")

df.x.coni <- IAA_simple %>%
  filter(microbe_level == "E. mundtii", diet_type == "control")

#make means from those data frames
ax.iaa <- mean(df.ax.iaa$iaa_avg)
x.iaa<- mean(df.x.iaa$iaa_avg)
ax.coni <- mean(df.ax.coni$iaa_avg)
x.coni <- mean(df.x.coni$iaa_avg)

#make standard erros from those data frames
#function for standard error
st.error <- function(x) sd(x) / sqrt(length(x))

ax.iaa.err <- st.error(df.ax.iaa$iaa_avg)
x.iaa.err <- st.error(df.x.iaa$iaa_avg)
ax.con.erri <- st.error(df.ax.coni$iaa_avg)
x.con.erri <- st.error(df.x.coni$iaa_avg)

#combine them all into one dataframe
microbe_treatment.i <- c("Axenic", "E. mundtii", "Axenic", "E. mundtii")
diet_treatment.i <- c("indole", "indole", "control", "control")
iaa_avg <- c(ax.iaa, x.iaa, ax.coni, x.coni)
sterr <- c(ax.iaa.err, x.iaa.err, ax.con.erri, x.con.erri)

df.iaa.control <- data.frame(microbe_treatment.i, diet_treatment.i, iaa_avg, sterr)
View(df.iaa.control)

#make levels so that it graphs as control diet first, then BA-spiked diet
df.iaa.control$diet_treatment.i <- factor(df.iaa.control$diet_treatment.i,
                                         levels = c("control", "indole"))
#final graph of the whole shebang
IAA_simple_bar <- ggplot(df.iaa.control, aes(x = diet_treatment.i, y = iaa_avg, fill = microbe_treatment.i))+
  geom_col(width=0.8, position = position_dodge())+
  geom_errorbar(aes(ymin= iaa_avg-sterr, ymax= iaa_avg+sterr), 
                width=0.5, linewidth=0.7, position=position_dodge(0.8))+
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "")+
  scale_x_discrete(labels = c("Unspiked", "SA-Spiked"))+
  labs(x = "Diet Type", y= "IAA Concentration (ng/g)") +
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_simple_bar
