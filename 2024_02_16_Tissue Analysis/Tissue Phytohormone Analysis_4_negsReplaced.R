###Analysis of Phytohormone in Different Tissues in CONTROL DIET
###For Thesis Proposal (and final manuscript, probably)
### Adam Scherr, 5/29/2024
### Version 4, 09/13/2024

#### Making Separate Datasets for Analysis ####
#start by importing the dataset
df <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2023_11_30_CommitteeMeeting/EntSoc_SeqData_2%20-%20tidy_data.csv")

#remove anything unnecessary: 
  #labeled, everything but avg concentrations, unneeded sample qualities, blanks
  #and make all the objects what you want them to be
library(dplyr)
library(ggplot2)
library(lme4)
library(ggthemes)
df.1 <- df %>%
  filter (microbe_level != "blank", nat_lab != "lab") %>%
  select(1:5, 9, 12, 15)

str(df.1) #now you know what needs to be made the correct object

df.1 <- df.1 %>%
  mutate_at(c("avg_ba_conc", "avg_sa_conc", "avg_iaa_conc"), as.numeric) %>%
  mutate_at(c("microbe_level", "diet_type", "tissue"), as.factor)
str(df.1)
View(df.1)

#now just make the names a little easier to type out
colnames(df.1)[c(6:8)] <- c("ba_avg", "sa_avg", "iaa_avg")

df.1$tissue <- factor(df.1$tissue,
                        levels = c("sal", "fat", "hem"))

fit.1 <- lm(ba_avg ~ microbe_level + tissue + diet_type, data = df.1)
summary(fit.1)

#now make the negative values into zeros
df.1$ba_avg[df.1$ba_avg < 0] <- 0
df.1$sa_avg[df.1$sa_avg < 0] <- 0
df.1$iaa_avg[df.1$iaa_avg < 0] <- 0


##now divide the dataset into sections for each phytohormone concentration
#tissues with high levels of benzoic acid
df.ba <- df.1 %>%
  filter(diet_type == "benzoic")

#tissues spiked with sa
df.sa <- df.1 %>%
  filter(diet_type == "salicylic")

#tissues spiked with iaa
df.iaa <- df.1 %>%
  filter(diet_type == "indole")

#tissues with control diet
df.control <- df.1 %>%
  filter(diet_type == "control") %>% 
  filter (id != "NC2H_1")
    #NC2H_1 and NC2H I think represent the same hemolymph sample
    #NC2H is prepared like all the other samples (100 microliters of hemolymph with 300 microliters PBS)
    #NC2H_1 was prepared as just 300 microliters of hemolymph, NO PBS added

##graphs with the control data for SA, BA, and IAA
control_babox <- ggplot(df.control, aes(x = tissue, y = ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha=0.82)+
  geom_point(size = 4, alpha =1, position=position_dodge(width=0.8), aes(color=microbe_level))+
  scale_fill_manual(values = c("#74a9cf", "#023858"), name = "",
                               labels = c("Axenic", "Xenic"))+
  scale_color_manual(values =c("#74a9cf", "#023858"), name = "",
                     labels = c("Axenic", "Xenic"))+
  labs(x = "Tissue", y = "BA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
control_babox

control_sabox <- ggplot(df.control, aes(x = tissue, y = sa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha=0.82)+
  geom_point(size = 4, alpha =1, position=position_dodge(width=0.8), aes(color=microbe_level))+
  scale_fill_manual(values = c("#66c2a4", "#00441b"), name = "",
                    labels = c("Axenic", "Xenic"))+
  scale_color_manual(values =c("#66c2a4", "#00441b"), name = "",
                     labels = c("Axenic", "Xenic"))+
  labs(x = "Tissue", y = "SA Concentration (ng/g)")+
  theme_clean() +
  scale_x_discrete(labels = c("Sal. Gland", "Fat Body", "Hemolymph"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
control_sabox

control_iaabox <- ggplot(df.control, aes(x = tissue, y = iaa_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha=0.82)+
  geom_point(size = 4, alpha =1, position=position_dodge(width=0.8), aes(color=microbe_level))+
  scale_fill_manual(values = c("#ef6548", "#7f0000"))+
  scale_color_manual(values =c("#ef6548", "#7f0000"))+
  # labs(x = "Diet Treatment", y = "BA Concentration (ng/g)")+
  theme_clean() +
  # scale_x_discrete(labels = c("Control", "BA-Spiked"))+
  theme(axis.text = element_text(size =15),
        panel.grid.major = element_line(size = 2),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
control_iaabox


#### Statistical Analysis Resources and Overview####
##Notes
#https://www.datacamp.com/tutorial/generalized-linear-models
#https://www.statology.org/interpret-glm-output-in-r/
#https://www.rensvandeschoot.com/tutorials/lme4/
#https://www.r-bloggers.com/2014/04/checking-glm-model-assumptions-in-r/
#I believe you want Gaussian (Normal) Distribution family? But you're data probs isn't normal so I don't know
      #would this mean you can just use a linear model, not necessarily a GLM?
      #https://www.youtube.com/watch?v=ddCO2714W-o
  #Binomial is for Yes/No data, Poisson is for frequency data (within a given interval) or count data


##Let's see what our data looks like, starting with the df.ba
hist(df.ba$ba_avg) #not normal, maybe a log-normal distribution? Poisson?
hist(df.sa$sa_avg) #POSSIBLY normal, with a lean to the left
hist(df.iaa$iaa_avg) #also possibly normal, not as skewed

hist(df.control$ba_avg) #possibly normal
hist(df.control$sa_avg) #skewed, left has bigger peak than the right
hist(log(df.control$sa_avg)) #logged sa_avg looks better for normality

#to tell which family (Gaussian, poison, etc. works best, use the residuals)
  #residuals are formed based on a linear regression (AKA line of best fit) for the data points
  #Residual = observed value - predicted value; essentially, how far off was your line of best fit
  #You will have one residual for each data point; better R-squared means smaller absolute values of residuals 
        #(because negative residual just means the prediction was too high for the true value)
  #Fun fact: the sum of all residuals will equal zero, and the mean of all residuals is 0

#assumptions for a GLM are based on this blog post:
  #https://www.r-bloggers.com/2014/04/checking-glm-model-assumptions-in-r/


###GLM for BA Concentration Data
library(lme4)


####GLMs for Interaction of microbe_level and tissue for each phytohormone####

#you need the emmeans package for the post-hoc tests (Tukey kramer, I think?)
library(emmeans)
library(multcomp)

##doing linear models followed by emmeans posthoc for control diet data
##starting with control diet ba data
cba_glm <- glm(ba_avg ~ microbe_level*tissue, family = gaussian(), df.control)

#check assumptions on control_ba data
par(mfrow=c(2,2))
plot(cba_glm)
  #looks relatively normal with QQ plot. Looks pretty much homoscedastic

#summary followed by emmeans posthoc on control_ba data
summary(cba_glm) #one significant interaction between fat body and other tissues
emmeans_cba <- emmeans(cba_glm, pairwise ~ (microbe_level*tissue), type = "response")
emmeans_cba

cld(emmeans(cba_glm, pairwise ~ microbe_level*tissue, adjust = "emmeans"),
    Letters =letters, alpha = 0.05) ##this one is right

cld(emmeans_cba, Letters = letters, alpha = 0.05) ##this one is also right


##now for control diet SA data
csa_glm <- glm(sa_avg ~ microbe_level*tissue, family = gaussian(), df.control)

#check assumptions
plot(csa_glm) #looks pretty normal and homoscedastic

#summary and emmeans posthoc
summary(csa_glm) #no significant interacitons
emmeans_csa <- emmeans(csa_glm, pairwise ~ (microbe_level*tissue), type = "response")

cld(emmeans(csa_glm, pairwise ~ microbe_level*tissue, adjust = "emmeans"),
    Letters =letters, alpha = 0.05) ##this one is right

cld(emmeans_csa, Letters = letters, alpha = 0.05) 


##for BA spiked diet tissue data
ba_glm <- glm(ba_avg ~ microbe_level + tissue + (microbe_level*tissue), 
              family = gaussian(), df.ba)
#check if assumptions are met
par(mfrow=c(2,2))
plot(ba_glm)
   #looks pretty homoscedastic, normality could be better

#look at the summary and emmeans posthoc
summary(ba_glm)
emmeans(ba_glm, pairwise ~ (microbe_level*tissue), type = "response")
#no significance anywhere

cld(emmeans(ba_glm, pairwise ~ microbe_level*tissue, adjust = "tukey"),
    Letters = letters, alpha=0.05)
  #output:
# microbe_level tissue emmean   SE df lower.CL upper.CL .group
# nonaxenic     hem      1286 1697 12    -2412     4985  a    
# axenic        fat      1509 1697 12    -2189     5208  a    
# axenic        hem      1621 1697 12    -2077     5320  a    
# axenic        sal      2261 1697 12    -1437     5960  a    
# nonaxenic     sal      5413 1697 12     1715     9112  a    
# nonaxenic     fat      5554 1697 12     1856     9253  a    



###now let's do this for SA data
sa_glm <- glm(sa_avg ~ microbe_level + tissue + (microbe_level*tissue), 
              family = gaussian(), df.sa)
#check if assumptions are met
par(mfrow=c(2,2))
plot(sa_glm)
  #looks pretty homoscedastic and a pretty straightline along the Q-Q plot

#look at the summary and emmeans posthoc
summary(sa_glm)
emmeans(sa_glm, pairwise ~ (microbe_level*tissue), type = "response")
  #axenic fat and nonaxenic hem under 0.05, and we are getting some at 0.09 and 0.10

cld(emmeans(sa_glm, pairwise ~ microbe_level*tissue, adjust = "tukey"),
    Letters = letters, alpha=0.05)
  #output:
# microbe_level tissue emmean   SE df lower.CL upper.CL .group
# nonaxenic     hem      3402 2340 12    -1696     8500  a    
# axenic        hem      5265 2340 12      167    10363  ab   
# nonaxenic     fat     10203 2340 12     5105    15301  ab   
# nonaxenic     sal     11238 2340 12     6140    16336  ab   
# axenic        sal     13570 2340 12     8472    18668  ab   
# axenic        fat     14889 2340 12     9791    19987   b   

###now lets do this for IAA data
iaa_glm <- glm(iaa_avg ~ microbe_level + tissue + (microbe_level*tissue),
               family = gaussian(), df.iaa)
#check if assumptions are met
par(mfrow=c(2,2))
plot(iaa_glm)
  #not the best fit for a striaght line on the Q-Q plot, but homoscedasticity looks kinda reasonable

#look at the summary and emmeans posthoc
summary(iaa_glm)
emmeans(iaa_glm, pairwise ~ (microbe_level*tissue), type = "response")
  #lots of signifiance! SCORE

cld(emmeans(iaa_glm, pairwise ~ microbe_level*tissue, adjust = "tukey"),
    Letters = letters, alpha=0.05)
  #output:
# microbe_level tissue emmean   SE df lower.CL upper.CL .group
# axenic        hem      4162 1649 12      569     7755  a    
# nonaxenic     hem      5996 1649 12     2403     9589  a    
# axenic        sal      7131 1649 12     3538    10724  a    
# axenic        fat     11800 1649 12     8207    15393  ab   
# nonaxenic     sal     16318 1649 12    12725    19911   bc  
# nonaxenic     fat     22953 1649 12    19361    26546    c  