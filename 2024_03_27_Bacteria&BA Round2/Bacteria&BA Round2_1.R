###Analysis of Phytohormones in a second BA run AND in Bacteria grown in 2XYT broth
###For Thesis Proposal (and final manuscript, probably)
### Adam Scherr, 3/29/2024
### Version 1

##upload your datasets
bact_data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_03_27_Bacteria%26BA%20Round2/Bacteria_Data_2024_03_27.csv")
ba_data <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2024_03_27_Bacteria%26BA%20Round2/BA_Round%202_2024_03_27.csv")

#load in packages
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(ggthemes)

#fix all that empty space that's there in these data files
bact_df <- bact_data %>%
  filter (ba_avg > 0) %>%
  select(1:10)

ba_df<- ba_data %>%
  filter(ba_avg >0) %>%
  select(1:7)

#make the objects the correct classes
str(bact_df)
str(ba_df)
ba_df <- ba_df %>%
  mutate_at(c("microbe_level", "tissue"), as.factor)
str(ba_df)
  
##bacterial graphs
bact_ba <- ggplot(bact_df, aes(x = id, y = ba_avg))+
  geom_bar(stat = "identity", fill = "blue")
bact_ba

bact_sa <- ggplot(bact_df, aes(x = id, y = sa_avg))+
  geom_bar(stat = "identity", fill = "green")
bact_sa

bact_iaa <- ggplot(bact_df, aes(x = id, y = iaa_avg))+
  geom_bar(stat = "identity", fill = "orange")
bact_iaa

grid.arrange(bact_ba, bact_sa, bact_iaa, ncol =1)

##ba data graph
BA_round2 <- ggplot(ba_df, aes(x = tissue, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
  scale_fill_manual(values = c("#74a9cf", "#023858"))+
  scale_color_manual(values = c("#74a9cf", "#023858")) 
  
BA_round2 #now THAT is not looking significant in the slightest

##running the stats on the BA round 2 data, just to be extra sure

#check for normality
hist(ba_df$ba_avg) #looks roughly normal, not ideal but good enough


#using Gaussian because that's the kind of data we have; which is a linear model
ba_lm <- lm(log(ba_avg) ~ microbe_level * tissue, ba_df)

#check if the linear model fits assumptions
par(mfrow = c(2, 2))
plot(ba_lm)
#the assumptions look good! The top left graph residuals look pretty well distributed
#and the top right graph showing normality is a pretty straight line, which means pretty normal distribution

#using emmeans for post-hoc analysis
library(emmeans)
summary(ba_lm)
emmeans(ba_lm, pairwise ~ (microbe_level*tissue), type = "response")
  #all p-values are above 0.05 EXCEPT for the comparison of fat body phytohormone levels
  #to salivary gland phytohormone levels. We even have p=1.0000 in the comparison of
  #axenic vs nonaxenic salivary glands. Wild how similarly they match up


# BA_tissue <- ggplot(data.ba, aes(x = tissue, y = ba_avg, fill = microbe_level))+
#   geom_boxplot(width=0.8, color = "black", alpha = 0.82)+
#   geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color= microbe_level)) +
#   scale_fill_manual(values = c("#74a9cf", "#023858"), name = "", 
#                     labels = c("Axenic", "Nonaxenic"))+
#   scale_color_manual(values = c("#74a9cf", "#023858"), name = "", 
#                      labels = c("Axenic", "Nonaxenic"))+
#   theme_clean() +
#   labs(x = "Tissue", y = "BA Concentration (ng/g)")+
#   scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
#   theme(axis.text = element_text(size =15),
#         panel.grid.major = element_line(size = 2),
#         axis.title = element_text(size=20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=15))


