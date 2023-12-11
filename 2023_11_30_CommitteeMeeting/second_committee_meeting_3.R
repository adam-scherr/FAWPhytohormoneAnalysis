###Analysis of BA, SA, and IAA Sequestration data, controled study
####For Committee Meeting and possibly general publication 
###Adam Scherr, 11/27/2023
###Version 3, 12/5/2023

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
#View(data.sal)

#2. data set for fat body
data.fat <- data.1 %>%
  filter(tissue == "fat")
#View(data.fat)

#3. data set for hemolymph
data.hem <- data.1 %>%
  filter(tissue == "hem")
#View(data.hem)

#4. data set for benzoic acid spiked diet only
data.ba <- data.1 %>%
  filter (diet_type == "benzoic")
#View(data.ba)

#5 data set sa spiked diet only
data.sa <- data.1 %>%
  filter(diet_type == "salicylic")
#View(data.sa)

#6 data set for iaa spiked diet only
data.iaa <- data.1 %>%
  filter(diet_type == "indole")
#View(data.iaa)

#7. data set for control diet only
data.control <- data.1 %>%
  filter(diet_type == "control")
#View(data.control)

##LET THE GRAPHING BEGIN!!
##firstly, the 9 graphs showing phytohormone levels in each tissue


#### salivary gland phytohormone levels #####
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

#### fat body phytohormone levels #####
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


#### hemolymph phytohormone levels #####
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

#### tissue graphs ####
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
#### legible graphs ####
##Graphs with legible concentrations for the salivary glands
##First, make datasets to work with by removing certain diet types
data.sansba <- data.sal %>%
  filter(diet_type != "benzoic")

data.sanssa <- data.sal %>%
  filter(diet_type != "salicylic")

data.sansiaa <- data.sal %>%
  filter(diet_type != "indole")

#now make the graphs
#legible BA levels in salivary glands
BA_legible <- ggplot(data.sansba, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "BA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "SA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_legible

#legible SA levels in salivary glands
SA_legible <- ggplot(data.sanssa, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "SA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "IAA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_legible

#legible IAA levels in salivary glands
IAA_legible <- ggplot(data.sansiaa, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_legible

#### STATS TIME! ####
###TIME to do some statistical analysis, woooo!!
#Goal: Use a two-way MANOVA (just to see the effect of microbes and diet type)
    #although, could it be worth it to also run a MANOVA with tissue type as a factor?

#run manova so you can test the residuals for all the assumptions
#lets run this test on just the salivary gland data. Don't complicate it yet
manova.test <- manova(cbind(ba_avg, iaa_avg, sa_avg) ~ microbe_level + diet_type, data=data.sal)
summary(manova.test)
  #p<0.05 only for microbe_level. Not for microbe_level

#now run anova on that manova output to get an idea of how each thing is affected
summary.aov(manova.test)
  #ba_avg: p<0.05 only for diet_type
  #iaa_avg: p<0.05 for both microbe_level and diet_type
  #sa_avg: p<0.05 only for diet_type

#test multivariate normality using Mardia's Skewness and Kurtosis Test
library(mvnormalTest)
grouped_data <- data.sal %>%
  group_by(microbe_level, diet_type)

mardia(grouped_data[, c("ba_avg", "iaa_avg", "sa_avg")])$mv.test
#getting a result of "NO" for all three outputs. This is not normally distributed data. Cool. (/s)
#time to perform a nonparametric test! I haven't really found a suitable test for this yet. Perhaps the committee has some insights

##We're running a PERMANOVA using Euclidean distance measures
##Preparing to run a PERMANOVA
#Used this resource to inform my decisions:
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/distance-measures/

#According to the article, Euclidean and Manhattan distance measures should work
#for my data because it is 1) symmetric (two variables that have a value of zero is
#just as meaningful as two variables that share a non-zero value) and 
#2)dependent variables are quantitative/continuous data (as opposed to categorical or binary)

#Euclidean distance measures need to be on the same scale. IAA goes 0-16,000 while BA
#goes from 0-10,000. I think this is still pretty much within the same scale 
#(not even off by a whole factor of 10)

#let's do this PERMANOVA, baby!
library(vegan)

#create species data matrix (our species are avg phytohormone concentrations), 
#so remove the "high" and "low" phytohormone rows
perm.data <- select(data.sal, !contains("high")) %>%
  select(!contains("low"))
#View(perm.data)

#only select the 3 columns containing your phytohormone concentration averages
species <- perm.data[, 9:11] 

#calculate the distances in your matrix
euc.distance <- vegdist(species, "euclidean", na.rm = TRUE)

#now run the permanova
permanova <- adonis2(euc.distance ~ microbe_level+diet_type,
                     permutations = 999,
                     method = "euclidean",
                     data=perm.data)

permanova #Significance! Pr(>F)<0.05 for BOTH microve_level and diet_ype. R squared values tell us that 
#microbe_level explains 2.42% of the variation, while diet_type explains 82.27% of the variation
#Now we need to use posthoc to see where specifically these effects are happening

##Run a posthoc test now. Used the following source:
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/#Implementation%20in%20R

#make the beta dispersion, maybe for just one group at a time?

#Microbe_level group:
betadisper_data1 <- betadisper(d=dist(perm.data[,9:11]), group = perm.data$microbe_level, 
                               type = "centroid")

betadisper_data1$distances
plot(betadisper_data1) #there is A LOT of overlap. Not a good sign for significance

#run anova on this to see significance
anova(betadisper_data1) #I'm seeing p=0.1027, NOT significant

#run this with adonis as a permutation based test (instead of anova), see if that makes more sense
adonis2(dist(betadisper_data1$distances) ~ perm.data$microbe_level+perm.data$diet_type)
#I got significance for both microbe_level (p=0.041) and diet_type (p=0.003),
#but for some reason, this test was purposely not run in the file "entsoc_2023_poster_v2" which 
  #I am basing this code off of. Why might that be?
  #OH, I see. It's because I defined the dataset "betadisper_data1" using microbe_level as my group,
    #NOT using a combination of microbe_level and diet_type. Let's try that out now

#microbe_level and diet_type as a single grouped dataset:
betadisper_data2 <- betadisper(d=dist(perm.data[,9:11]), group = perm.data$microbe_leve+perm.data$diet_type,
                               type= "centroid")
#well that spits out an error. cool. I guess that's why I didn't use the adonis2 output...lame


#now run pairwise comparison to tease out significance for the microbe_level grouped dataset. 
  #This is Tukey HSD
TukeyHSD(betadisper_data1)
#P-value is 0.1027, NO SIGNIFICANCE :(

#diet_type group:
betadisper_data3 <- betadisper(d=dist(perm.data[,9:11]), group = perm.data$diet_type, 
                               type = "centroid")
plot(betadisper_data3)
anova(betadisper_data3) #we got significance!
TukeyHSD(betadisper_data3) 
#p<0.05 for interaction between SA diet and control diet
#p<0.05 for interaction between IAA diet and control diet

##those are the only two significant interactions. I suppose that means that
  #we can't say any other graph has significance? THAT MAKES NO SENSE! The BA diet definitely
  #causes a wildly different level of phytohormones when compared to the control diet (as one
  #example of an interaction that is supposedly "not significant").

#I think we can chalk up the weird statistical results as the fault of the incorrect statistical
#analysis for the job. Flor Acevedo used one-way ANOVAs for some of her 2019 paper's statistics.
#But crucially, for analyzing the levels of phytohormones in the caterpillar saliva,
#she used a "two factor factorial design" in Minitab 18. WHAT THE FIRETRUCK DOES THAT EVEN MEAN?
#perhaps this resource will be of service:
#https://online.stat.psu.edu/stat503/lesson/5/5.1 

####Redoing the Stats (Dec. 5th, 2023)####
##start by making datasets that encapsulate just the data you want
#you want salivary gland data, then each data set will be restricted to the spiked-diet of interest
#and a control

#BA and BA-spiked diet dataset
BA_simple <- data.sal %>%
  filter(diet_type != "salicylic", diet_type != "indole") %>%
  select(1:9)

#SA and SA-spiked diet dataset
SA_simple <- data.sal %>%
  filter(diet_type != "benzoic", diet_type != "indole") %>%
  select(1:8, 12)

#IAA and IAA_spiked diet dataset
IAA_simple <- data.sal %>%
  filter(diet_type != "benzoic",  diet_type != "salicylic") %>%
  select(1:8, 15)

#Let's visualize these graphs
BA_simple_box <- ggplot(BA_simple, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "")+
  scale_color_manual(values = c("#74a9cf", "#023858"), name = "")+
  labs(x = "Diet Type", y = "BA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
BA_simple_box

SA_simple_box <- ggplot (SA_simple, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "")+
  scale_color_manual(values = c("#66c2a4", "#00441b"), name = "")+
  labs(x = "Diet Type", y = "SA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
SA_simple_box

IAA_simple_box <- ggplot (IAA_simple, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#ef6548", "#7f0000"), name = "")+
  scale_color_manual(values = c("#ef6548", "#7f0000"), name = "")+
  labs(x = "Diet Type", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
IAA_simple_box

###lets dig into the stats of it all now
##BA_simple stats:
#start by seeing if the data looks normally distrubted
hist(BA_simple$ba_avg) #data looks skewed, not bell curved
qqnorm(BA_simple$ba_avg) #doesn't look much like a straight, diagonal line. Not a good sign of normality
qqline(BA_simple$ba_av)

####my attempt at PERMANOVA that did not work####
#just to be sure, lets do a little Shapiro-Wilk test
shapiro.test(BA_simple$ba_avg) #p<0.05, so data is NOT normally distributed

#test the multivariate normality using Mardia's multivariate skewness and kurtosis test
library(mvnormalTest)
grouped_data_BA <- BA_simple %>%
  group_by(microbe_level, diet_type)

mardia(grouped_data_BA[, c("ba_avg")])$mv.test
  #we're getting a NO for multivariate normality. A shame. Let's do PERMANOVA
library(vegan)

#create species data matrix (our species are avg phytohormone concentrations)
#only select the column containing your phytohormone concentration average
species_ba <- grouped_data_BA[, 9] 

#calculate the distances in your matrix
euc.distance_ba <- vegdist(species_ba, "euclidean", na.rm = TRUE)

#now run the permanova
permanova_ba <- adonis2(euc.distance_ba ~ microbe_level+diet_type,
                     permutations = 999,
                     method = "euclidean",
                     data=grouped_data_BA)

permanova_ba #Significance only for diet_type. No significance for microbe_level

##Run a posthoc test now. Used the following source:
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/#Implementation%20in%20R

#make the beta dispersion, maybe for just one group at a time?
#Microbe_level group:
betadisper_ba1 <- betadisper(d=dist(grouped_data_BA[,8:9]), group = grouped_data_BA$microbe_level, 
                               type = "centroid")

betadisper_ba1$distances
plot(betadisper_ba1) #there is A LOT of overlap. Not a good sign for significance

#run anova on this to see significance
anova(betadisper_data1) #I'm seeing p=0.1027, NOT significant

#run this with adonis as a permutation based test (instead of anova), see if that makes more sense
adonis2(dist(betadisper_data1$distances) ~ perm.data$microbe_level+perm.data$diet_type)
#I got significance for both microbe_level (p=0.041) and diet_type (p=0.003),
#but for some reason, this test was purposely not run in the file "entsoc_2023_poster_v2" which 
#I am basing this code off of. Why might that be?
#OH, I see. It's because I defined the dataset "betadisper_data1" using microbe_level as my group,
#NOT using a combination of microbe_level and diet_type. Let's try that out now

#microbe_level and diet_type as a single grouped dataset:
betadisper_data2 <- betadisper(d=dist(perm.data[,9:11]), group = perm.data$microbe_leve+perm.data$diet_type,
                               type= "centroid")
#well that spits out an error. cool. I guess that's why I didn't use the adonis2 output...lame


#now run pairwise comparison to tease out significance for the microbe_level grouped dataset. 
#This is Tukey HSD
TukeyHSD(betadisper_data1)
#P-value is 0.1027, NO SIGNIFICANCE :(

#diet_type group:
betadisper_data3 <- betadisper(d=dist(perm.data[,9:11]), group = perm.data$diet_type, 
                               type = "centroid")
plot(betadisper_data3)
anova(betadisper_data3) #we got significance!
TukeyHSD(betadisper_data3) 
#p<0.05 for interaction between SA diet and control diet
#p<0.05 for interaction between IAA diet and control diet

##those are the only two significant interactions. I suppose that means that
#we can't say any other graph has significance? THAT MAKES NO SENSE! The BA diet definitely
#causes a wildly different level of phytohormones when compared to the control diet (as one
#example of an interaction that is supposedly "not significant").

####my attempt at Kruskal-Wallis####
kruskal.test(microbe_level ~ ba_avg, data = BA_simple)
 #Kruskal-Wallis test will only analze one explanatory variable's effect on one response variable
  #(but you wanted it to look at two explanatory variable's effects on one response variable)
#I'VE GOT IT! Just boil it down to a t-test or nonparametric equivalent. 
#Just look if the axenic and nonaxenic of the SAME DIET are different from each other!

####T-Test####
##let's start by making datasets, like we did last time
#dataset for BA
t.test.ba <- data.sal %>%
  filter(diet_type == "benzoic") %>%
  select(1:9)
t.test.ba

#now visualize the data
plot(t.test.ba$ba_avg~t.test.ba$microbe_level)

#observe if it looks normal
plot(t.test.ba$ba_avg)
qqnorm(t.test.ba$ba_avg)
shapiro.test(t.test.ba$ba_avg) #p=0.22, that's normally distributed!!! FINALLY

#group data to check the three assumptions: no outliers, normal distribution, equal variance
#resource for the following t-test info: "https://www.datanovia.com/en/lessons/t-test-in-r/" 
grouped.ba <- t.test.ba %>%
  group_by(microbe_level)
#are there outliers?
library(rstatix)
grouped.ba %>%
  identify_outliers(ba_avg) #I think no outliers?

#is it normal?
grouped.ba %>%
  shapiro_test(ba_avg) #we've got normality

#are variances equal?
library(car)
#levene_test(grouped.ba, ba_avg ~ microbe_level)
leveneTest(ba_avg ~ microbe_level, data = grouped.ba)
  #p=0.2 which is greater than 0.05, thus the variances ARE equal (the variances are 
      #not different to a statistically significant degree)

#all assumptions being met, we will use an independent, Student's t-test
#R's t_test() function defaults to a Welch's t-test, but since you're variances are equal,
#we will go with the student's t-test. We just have to tell R in the t_test() what we want
t.test(ba_avg ~ microbe_level, grouped.ba, var.equal = TRUE)
  #we get a p-value of 0.3258. Axenic and nonaxenic are NOT statistically different!

###now lets do that for SA
#make the SA dataset
t.test.sa <- data.sal %>%
  filter(diet_type == "salicylic") %>%
  select(1:8, 12)

#visualize the data
plot(t.test.sa$sa_avg ~ t.test.sa$microbe_level)

#observe to see if it looks normally distributed
qqnorm(t.test.sa$sa_av) #looks pretty normal, tbh

#make the data in grouped form for the next three assumptions tests
grouped.sa <- t.test.sa %>%
  group_by(microbe_level)

#are there outliers?
grouped.sa %>%
  identify_outliers(sa_avg) #I think no outliers?

#is it normally distributed?
grouped.sa %>%
  shapiro_test(sa_avg) #they are normally distrivuted, p>0.05

#equal variance?
leveneTest(sa_avg ~ microbe_level, data = grouped.sa)
  #p>0.05, they DO have equal variances

#run the student's t-test!
t.test(sa_avg ~ microbe_level, data = grouped.sa, var.equal = TRUE)
  #p = 0.5166, axenic and nonaxenic are NOT statistically different


###lastly, lets do it for IAA
#make the dataset
t.test.iaa <- data.sal %>%
  filter(diet_type == "indole") %>%
  select(1:8, 15)

#visualize the data
plot(t.test.iaa$iaa_avg ~ t.test.iaa$microbe_level)

#visualize the normality
qqnorm(t.test.iaa$iaa_avg)

#make the grouped dataset to run the next three assumptions tests
grouped.iaa <- t.test.iaa %>%
  group_by(microbe_level)

#are there outliers?
grouped.iaa %>%
  identify_outliers(iaa_avg) #maybe no outliers? I think this function doesn't actually work

#is it normal distribution?
grouped.iaa %>%
  shapiro_test(iaa_avg) #P>0.05, we have normal distribution

#is variance equal?
leveneTest(iaa_avg~microbe_level, data=grouped.iaa)
  #p>0.05, we have equal variance

##lets run our student's t-test
t.test(iaa_avg~microbe_level, data=grouped.iaa, var.equal=TRUE)
  #p = 0.01647, axenic and nonaxenic values ARE statistically different!