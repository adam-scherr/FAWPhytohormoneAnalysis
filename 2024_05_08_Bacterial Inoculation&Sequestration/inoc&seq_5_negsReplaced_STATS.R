### Analysis of phytohormone contents in salivary glands of bacteria-inoculated FAW
### Adam Scherr, May 9th, 2024
### Version 5, FC-inoculated excluded AND negatives replaced with zeros


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

#make datasets for each phytohormone in spiked diet
ba.data <- data.1 %>%
  select(1:6) %>%
  filter(diet_type == "benzoic")

sa.data <- data.1 %>%
  select(1:3, 7:9) %>%
  filter(diet_type == "salicylic")

iaa.data <- data.1 %>%
  select(1:3, 10:12) %>%
  filter(diet_type=="indole")

#datasets for each phytohormone in unspiked (control) diet
unspiked.ba.data <- data.1 %>%
  select(1:6) %>%
  filter(diet_type == "control")

unspiked.sa.data <- data.1 %>%
  select(1:3, 7:9)%>%
  filter(diet_type == "control")

unspiked.iaa.data <- data.1 %>%
  select(1:3, 10:12)%>%
  filter(diet_type == "control")

####graphs of the data####
##lets see what these datasets look like graphically
  #still need to have axis titles changes,legends removed, and theme added
ba_box <- ggplot(ba.data, aes(x = microbe_level, y=ba_avg, fill=microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#08306b"), name = "")+
  scale_color_manual(values = c("#9ecae1", "#08306b"), name="") +
  labs(x= "Microbe Treatment", y= "BA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
  
ba_box

sa_box <- ggplot(sa.data, aes(x = microbe_level, y = sa_avg, fill = microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b",  "#00441b"), name = "")+
  scale_color_manual(values = c("#a1d99b",  "#00441b"), name = "")+
  labs(x = "Microbe Treatment", y = "SA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
sa_box

iaa_box <- ggplot(iaa.data, aes(x=microbe_level, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b",  "#7f2704"), name="")+
  scale_color_manual(values = c("#fdae6b", "#7f2704"), name="")+
  labs(x = "Microbe Treatment", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
iaa_box


##graphs of the unspiked diet data
unspiked_ba_box <- ggplot(unspiked.ba.data, aes(x = microbe_level, y=ba_avg, fill=microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1",  "#08306b"))+
  scale_color_manual(values = c("#9ecae1", "#08306b"))+
  labs(x = "Microbe Treatment", y = "BA Concentration (ng/g)")+
  theme_clean() +
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,335))

unspiked_sa_box <- ggplot(unspiked.sa.data, aes(x = microbe_level, y = sa_avg, fill = microbe_level)) +
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#a1d99b",  "#00441b"))+
  scale_color_manual(values = c("#64bc6e",  "#64bc6e"))

unspiked_iaa_box <- ggplot(unspiked.iaa.data, aes(x=microbe_level, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#fdae6b",  "#7f2704"))+
  scale_color_manual(values = c("#fd8d3c", "#fd8d3c"))

unspiked_ba_box
unspiked_sa_box
unspiked_iaa_box
####statistical analysis_One-way ANOVA (just comparing axenic and E. mundtii within same diet####
#Use one-way ANOVA because we have one categorical explanatory variable (microbe_level) and want to 
  #compare continuous response variable (phytohormone concentration) between 3 
#groups (axenic, Enterococcus mundtii -inoculated, and Frass Community-inoculated)

##testing ANOVA assumptions
#resource: https://www.statology.org/one-way-anova-r/

##checking assumptions for BA data
ba.model <- aov(ba_avg ~ microbe_level, data = ba.data)

#QQ plot for BA
par(mfrow = c(2, 2))
plot(ba.model)
 #QQ-plot looks like it fits the diagonal line very well, data is normal
 #the distribution looks roughly equal, but a Levene's test will show for certain

#Levene's test for equal variances for BA
library(car)
leveneTest(ba_avg ~ microbe_level, data = ba.data)
  #p>0.05, so we DO have equal variance, we passed both our assumptions for BA

##running ANOVA for BA data
summary(ba.model)
  #p = 0.0856, nearly significant for p<0.05 significance
  #I think I should run the Tukey-HSD posthoc test anyway to see where the significance is


##checking assumptions for SA data
sa.model <- aov(sa_avg ~ microbe_level, data = sa.data)

#QQ-plot for SA
par(mfrow = c(2, 2))
plot(sa.model)
  #QQ-plot is not perfectly fitting along diagonal line, let's double check with Shapiro-Wilk test

#Shapiro-Wilk test for SA
shapiro.test(sa.data$sa_avg)
  #p=0.83008, far above 0.05. We have normality

#Levene's test for equal variances for SA
leveneTest(sa_avg~microbe_level, data = sa.data)
  #p=0.564, far above 0.05. Variances are equal

##running ANOVA for SA data
summary(sa.model)
  #p=0.599. There is no significance of microbe_level on SA concentration


##checking assumptions for IAA data
iaa.model <- aov(iaa_avg~microbe_level, data = iaa.data)

##QQ-plot for IAA
par(mfrow = c(2, 2))
plot(iaa.model)

#Shapiro-Wilk test for IAA
shapiro.test(iaa.data$iaa_avg)
  #p = 0.2412, above 0.05 so we have normality

#Levene's test for IAA
leveneTest(iaa_avg ~ microbe_level, data = iaa.data)
  #p=0.2302, above 0.05 so variances are equal

##running ANOVA for IAA data
summary(iaa.model)
  #p=0.0848, above 0.05, microbe_level is NOT SIGNIFICANT for IAA concentration
  #this difference becomes significiant when you add in the frass 
    #community inoculated data. Just so you know

####statistical analysis_Post Hoc Tukey's Test####
##post-hoc for BA data (nearly significant p-value)
TukeyHSD(ba.model, conf.level=0.95)
  #p=0.085569 for the difference between E. mundtii and axenic

##post-hoc for SA data
TukeyHSD(sa.model, conf.level=0.95)
  #p is far greater than 0.05. NO SIGNIFICANCE

#post-hoc for IAA data
TukeyHSD(iaa.model, conf.level=0.95)
  #p=0.0848201; this interaction is even less significant (p=0.11) when you include
    #the frass community data

####unspiked diet ANOVA and posthoc####
#unspiked.ba.data, unspiked.sa.data, unspiked.iaa.data
##checking assumptions for BA data
unspiked.ba.model <- aov(ba_avg ~ microbe_level, data = unspiked.ba.data)

#QQ plot for BA
par(mfrow = c(2, 2))
plot(unspiked.ba.model)
#QQ-plot looks like it fits the diagonal line very well, data is normal
#the distribution looks roughly equal, too
summary(unspiked.ba.model)
  #we have significance! p=0.00325

##assumptions and anova for SA data
unspiked.sa.model <- aov(sa_avg~microbe_level, data = unspiked.sa.data)
par(mfrow=c(2,2))
plot(unspiked.sa.model)
  #looks roughly normal and varianced equally
summary(unspiked.sa.model)
  #NOT significant, p=0.29

##assumptions and anova for IAA data
unspiked.iaa.model <- aov(iaa_avg~microbe_level, data = unspiked.iaa.data)
plot(unspiked.iaa.model)
  #looks roughly normal and equally varianced
summary(unspiked.iaa.model)
  #NOT significant, p=0.261

##post hoc analysis using TukeyHSD
TukeyHSD(unspiked.ba.model, conf.level = 0.95)
  #SIGNIFICANCE for interaction between E. mund and axenic (p=0.003251)
TukeyHSD(unspiked.sa.model, conf.level=0.95) #no significance (p=0.29)
TukeyHSD(unspiked.iaa.model, conf.level=0.95) #no significance (p= 0.261)


####Two-way ANOVA for phytohormone-spiked and unspiked salivary glands####
data.1$diet_type <- factor(data.1$diet_type,
                           levels = c("control", "benzoic", "salicylic", "indole"))

unspikedandba.data <- data.1 %>%
  filter(diet_type != "indole") %>%
  filter(diet_type != "salicylic")

unspikedandsa.data <- data.1 %>%
  filter(diet_type != "indole") %>%
  filter(diet_type != "benzoic")

unspikedandiaa.data <- data.1 %>%
  filter(diet_type != "benzoic")%>%
  filter(diet_type != "salicylic")




##make pretty graphs of it all

unspikedandba_box <- ggplot(unspikedandba.data, aes(x = diet_type, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha=0.82)+
  geom_point(size = 4, alpha =1, position=position_dodge(width=0.8), aes(color=microbe_level))+
  scale_fill_manual(values = c("#9ecae1",  "#08306b"))+
  scale_color_manual(values =c("#9ecae1",  "#08306b"))+
  labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
  scale_y_log10(n.breaks = 6)+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "BA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
unspikedandba_box

unspikedandsa_box <- ggplot(unspikedandsa.data, aes(x = diet_type, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha=0.82)+
  geom_point(size = 4, alpha =1, position=position_dodge(width=0.8), aes(color=microbe_level))+
  scale_fill_manual(values = c("#a1d99b",  "#00441b"))+
  scale_color_manual(values =c("#a1d99b",  "#00441b"))+
 labs(x = "Diet Treatment", y = "SA Concentration (ng/g)")+
  scale_y_log10(n.breaks = 6)+
  theme_clean() +
 scale_x_discrete(labels = c("Control", "SA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
unspikedandsa_box

unspikedandiaa_box <- ggplot(unspikedandiaa.data, aes(x = diet_type, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha=0.82)+
  geom_point(size = 4, alpha =1, position=position_dodge(width=0.8), aes(color=microbe_level))+
  scale_fill_manual(values = c("#fdae6b",  "#7f2704"))+
  scale_color_manual(values =c("#fdae6b",  "#7f2704"))+
  labs(x = "Diet Treatment", y = "IAA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Control", "IAA-Spiked"))+
  scale_y_log10(n.breaks = 6)+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
unspikedandiaa_box

##now run an ANOVA on this. I trust you!
##assumptions testing
#for BA and unspiked diet
unspikedandba.model <- aov(ba_avg ~ microbe_level*diet_type, data = unspikedandba.data)
par(mfrow=c(2,2))
plot(unspikedandba.model) #not the best normality, but variance looks pretty good

#for SA and unspiked diet
unspikedandsa.model <- aov(sa_avg ~ microbe_level*diet_type, data = unspikedandsa.data)
plot(unspikedandsa.model) #not the best normality, but variance looks pretty good

#for IAA and unspiked diet
unspikedandiaa.model <- aov(iaa_avg ~ microbe_level*diet_type, data = unspikedandiaa.data)
plot(unspikedandiaa.model) #not the best normality, but variance looks pretty good

##get the two-way anova outputs
summary(unspikedandba.model)
  #significance for microbe_level AND diet_type (and nearly their interaction, p=0.0578)

summary(unspikedandsa.model)
  #significance for diet_type only (p=0.000355)

summary(unspikedandiaa.model)
  #significance for diet_type and NEARLY significant for microbe_level and 
    #interaction betwixt them (p=0.0510 and p=0.0533 respectively)


##post hoc with TukeyHSD, followed by compact letter display output
library(lme4)
library(emmeans)
library(multcomp)

#BA and unspiked diet data
TukeyHSD(unspikedandba.model, conf.level=0.95)
cld(emmeans(unspikedandba.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean   SE df lower.CL upper.CL .group
# Axenic        control      124 1284  8    -2837     3085  a    
# E. mundtii    control      275 1284  8    -2687     3236  a    
# Axenic        benzoic     7492 1284  8     4531    10453   b   
# E. mundtii    benzoic    13326 1284  8    10365    16287    c 


#SA and unspiked diet data
TukeyHSD(unspikedandsa.model, conf.level=0.95)
cld(emmeans(unspikedandsa.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)

#output:
# microbe_level diet_type  emmean   SE df lower.CL upper.CL .group
# Axenic        control     -40.8 4141  8    -9589     9507  a    
# E. mundtii    control      63.8 4141  8    -9484     9612  a    
# Axenic        salicylic 22148.2 4141  8    12600    31696   b   
# E. mundtii    salicylic 26865.1 4141  8    17317    36413   b 


#IAA and unspiked diet data
TukeyHSD(unspikedandiaa.model, conf.level=0.95)
cld(emmeans(unspikedandiaa.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"), 
    Letters =letters, alpha = 0.05)
      #this finding is also found when you include FC-inoculated data

#output:
# microbe_level diet_type  emmean   SE df lower.CL upper.CL .group
# E. mundtii    control      41.3 3609  8    -8281     8363  a    
# Axenic        control     144.9 3609  8    -8177     8467  a    
# E. mundtii    indole    26121.6 3609  8    17799    34444   b   
# Axenic        indole    42574.6 3609  8    34252    50897    c 


