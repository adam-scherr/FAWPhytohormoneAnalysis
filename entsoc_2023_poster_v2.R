###BA, SA, and IAA sequestration data, from preliminary data
###For Ent Soc National 2023 Poster
###Adam Scherr, 10/30/2023

prelim_data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/SA_IAA_BA_prelim_data")
View(prelim_data)
library(dplyr)

#change object type so numeric data is numeric and treatments are factors
prelim_data <- prelim_data %>%
  mutate_at(c("ba_avg", "ba_high", "ba_low",
              "sa_avg", "sa_high", "sa_low",
              "iaa_avg", "iaa_high", "iaa_low"), as.numeric) %>%
  mutate_at(c("microbe_level", "diet_type"), as.factor) %>%
  filter(microbe_level != "diet", microbe_level != "mix")

str(prelim_data) #success, things are the object types that they should be

#now change the levels so graphs aren't alphabetically ordered
prelim_data$diet_type <- factor(prelim_data$diet_type, 
                                 levels = c("control", "benzoic", "SA", "IAA"))

levels(prelim_data$diet_type) #it worked, looks good

#make preliminary box and whisker plot to see the spread
library(ggplot2)
BA_box <- ggplot(prelim_data, aes(x = microbe_level, y= ba_avg, 
                                  color = diet_type))+
  geom_boxplot()+
  geom_jitter(size = 3, width = 0.1)
BA_box

SA_box <- ggplot(prelim_data, aes(x = microbe_level, y= sa_avg, 
                                  color = diet_type))+
  geom_boxplot()+
  geom_jitter(size = 3, width = 0.01)
SA_box

IAA_box <- ggplot(prelim_data, aes(x = microbe_level, y= iaa_avg, 
                                   color = diet_type))+
  geom_boxplot()+
  geom_jitter(size = 3, width = 0.1)
IAA_box

##################
##let's make some nicer graphs, the kind for your poster
library(ggdark)
library(ggthemes)
#benzoic = purple, #1
#control = orange, #2
#IAA = teal, #3
#brown = SA, #4
#color scheme options: autumn, purple, pink
#fall: grandbudapest1 in wesanderson package

#SA Scatter Plot
SA_point1 <- ggplot(prelim_data, aes(x= microbe_level, y = sa_avg, color = diet_type))+
  geom_pointrange(aes(ymin = sa_low, ymax = sa_high),
                  position=position_jitter(width=0.2), size =1.5)+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan", "sienna4"), name="Diet Type")+
  labs(x= "", y= "SA Concentration (ng/g)")+
  theme_clean(base_size=14)
SA_point1

#IAA Graph Scatter Plot
IAA_point1 <- ggplot(prelim_data, aes(x= microbe_level, y = iaa_avg, color = diet_type))+
  geom_pointrange(aes(ymin = iaa_low, ymax = iaa_high),
                  position=position_jitter(width=0.2), size =1.5)+
  scale_color_manual(values = c("darkmagenta", "chocolate1","cyan", "sienna4"), name="Diet Type")+
  labs(x= "", y= "IAA Concentration (ng/g)")+
  theme_clean(base_size=14)
IAA_point1

#Box and whisker plots, export with width 1000 and height 650
BA_box <- ggplot(prelim_data, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "grey34", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#fde0dd", "#7a0177"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#fde0dd", "#7a0177"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "BA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =25),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=30),
        legend.title = element_blank(),
        legend.text = element_text(size=25),
        legend.position = c(0.88,0.85))
BA_box

SA_box <- ggplot(prelim_data, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "grey34", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#fde0dd", "#7a0177"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#fde0dd", "#7a0177"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "SA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =25),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=30),
        legend.title = element_blank(),
        legend.text = element_text(size=25),
        legend.position = c(0.88, 0.85))
SA_box

IAA_box <- ggplot(prelim_data, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "grey34", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#fde0dd", "#7a0177"), name = "", 
                    labels = c("Axenic", "Nonaxenic"))+
  scale_color_manual(values = c("#fde0dd", "#7a0177"), name = "", 
                     labels = c("Axenic", "Nonaxenic"))+
  labs(x = "Diet Type", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA", "SA", "IAA"))+
  theme(axis.text = element_text(size =25),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=30),
        legend.title = element_blank(),
        legend.text = element_text(size=25),
        legend.position = c(0.88, 0.23))
IAA_box
  
#################
#lets take a crack at some statistics
  #Goal: Use a two-way MANOVA or its nonparametric equivalent (PERMANOVA)

#run manova so you can test the residuals for all the assumptions
manova.test <- manova(cbind(ba_avg, iaa_avg, sa_avg) ~ microbe_level + diet_type, data=prelim_data)
summary(manova.test)
#getting p-values below 0.05 for diet_type, but not for microbe_level
#now run anova on that manova output to get an idea of how each thing is affected
summary.aov(manova.test)
#ba_avg: only diet type has significant effect
#iaa_avg: BOTH microbiome and diet type have Pr(>F) under 0.05, so significant
#sa_avg: only diet type has signficiant effect

#test multivariate normality using Mardia's Skewness and Kurtosis Test
library(mvnormalTest)
grouped_data <- prelim_data %>%
  group_by(microbe_level, diet_type)

mardia(grouped_data[, c("ba_avg", "iaa_avg", "sa_avg")])$mv.test
#so much for that plan. We totally failed both the Skewness and Kurtosis test (result = NO)
#Looks like it's time to perform a non-parametric test!

##We're running a PERMANOVA using Euclidean distance measures
##Preparing to run a PERMANOVA
#Used this resource to inform my decisions:
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/distance-measures/

#According to the article, Euclidean and Manhattan distance measures should work
#for my data because it is 1) symmetric (two variables that have a value of zero is
#just as meaningful as two variables that share a non-zero value) and 
#2)dependent variables are quantitative/continuous data (as opposed to categorical or binary)

#Euclidean distance measures need to be on the same scale. IAA goes 0-20,000 while BA
  #goes from 0-2,500. I think this is still pretty much within the same scale 
  #(only off by a factor of 10 as opposed to a factor of like 100)

#let's do this PERMANOVA, baby!
library(vegan)

#create species data matrix (our species are avg phytohormone concentrations)
perm.data <- select(prelim_data, !contains("high")) %>%
  select(!contains("low"))
View(perm.data)

species <- perm.data[, 5:7] #only select the 3 columns containing your phytohormone concentration averages

#calculate the distances in your matrix
euc.distance <- vegdist(species, "euclidean", na.rm = TRUE)

#now run the permanova
permanova <- adonis2(euc.distance ~ microbe_level+diet_type,
                     permutations = 999,
                     method = "euclidean",
                     data=perm.data)
permanova #we're seeing significance in BOTH independent variables! Pr(>F) is less than 0.05!
  #R-squared values tell us that microbe_level explains 2.20% of the variation, and
  #diet_type explains 89.6% of the variation

##Run a posthoc test now. Used the following source:
  #https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/#Implementation%20in%20R

#make the beta dispersion, maybe for just one group at a time?
a
#Microbe_level group:
betadisper_data1 <- betadisper(d=dist(perm.data[,5:7]), group = perm.data$microbe_level, 
                              type = "centroid")

betadisper_data1$distances
plot(betadisper_data1)

#run anova on this to see significance
anova(betadisper_data1) #I'm seeing p=0.4762, no significance

#run this with adonis as a permutation based test (instead of anova), see if that makes more sense
#adonis2(dist(betadisper_data$distances) ~ perm.data$microbe_level+perm.data$diet_type)
  #I got no significance for microbe_level, but YES SIGNIFIANCE for diet_type!!

#now run pairwise comparison to tease out significance. This is Tukey HSD
TukeyHSD(betadisper_data1)
  #P-value is 0.4762, NO SIGNIFICANCE

#diet_type group:
betadisper_data2 <- betadisper(d=dist(perm.data[,5:7]), group = perm.data$diet_type, 
                               type = "centroid")
plot(betadisper_data2)
anova(betadisper_data2) #we got significance!
TukeyHSD(betadisper_data2) #looks to me like we have significance where it counts

#What you did for stats:
#ran PERMANOVA on non-parametric data, that's when you found 89.6% of variance explained by diet_type
#ran post-hoc based on betadispersion and TukeyHSD to make pairwise comparisons

#you can't just use a multivariate kruskal wallis because R doesn't have one of those
#and it is a complex statistical process
