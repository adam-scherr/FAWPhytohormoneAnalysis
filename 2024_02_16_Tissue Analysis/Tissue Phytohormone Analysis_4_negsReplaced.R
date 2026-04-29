###Analysis of Phytohormone in Different Tissues in ALL FOUR DIET TYPES
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
  dplyr::select(1:5, 9, 12, 15)

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
####Control Diet Tissue Analysis####
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
# output:
# microbe_level tissue emmean   SE df lower.CL upper.CL .group
# nonaxenic     sal      12.7 37.9 12   -69.88     95.3  a    
# axenic        sal      23.5 37.9 12   -59.08    106.1  a    
# nonaxenic     hem      69.4 37.9 12   -13.22    152.0  ab   
# axenic        hem      84.4 37.9 12     1.83    167.0  ab   
# nonaxenic     fat     236.3 37.9 12   153.74    318.9   bc  
# axenic        fat     281.4 37.9 12   198.80    364.0    c  

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
#output:
# microbe_level tissue emmean   SE df lower.CL upper.CL .group
# axenic        hem       0.0 32.9 12    -71.8     71.8  a    
# axenic        fat      11.1 32.9 12    -60.7     82.8  ab   
# axenic        sal      11.5 32.9 12    -60.3     83.3  ab   
# nonaxenic     hem      31.4 32.9 12    -40.4    103.1  ab   
# nonaxenic     sal      85.5 32.9 12     13.8    157.3  ab   
# nonaxenic     fat     157.1 32.9 12     85.3    228.8   b 
##now for control diet IAA data
hist(df.control$iaa_avg) #looks like Poisson
hist(log(df.control$iaa_avg)) #eh, not quit more normal
  #but this doesn't matter as long as the glm fits assumptions

ciaa_glm <- glm(iaa_avg ~ microbe_level*tissue, family = gaussian(), df.control)

#check assumptions
plot(ciaa_glm) #looks pretty normal and homoscedastic

#summary and emmeans posthoc
summary(ciaa_glm) #no significant interacitons
emmeans_ciaa <- emmeans(ciaa_glm, pairwise ~ (microbe_level*tissue), type = "response")

cld(emmeans(ciaa_glm, pairwise ~ microbe_level*tissue, adjust = "emmeans"),
    Letters =letters, alpha = 0.05) ##this one is right

cld(emmeans_ciaa, Letters = letters, alpha = 0.05) 
#output:
# microbe_level tissue emmean   SE df lower.CL upper.CL .group
# axenic        hem      0.00 7.62 12  -16.602     16.6  a    
# axenic        sal      0.00 7.62 12  -16.602     16.6  a    
# nonaxenic     hem      0.00 7.62 12  -16.602     16.6  a    
# nonaxenic     fat      2.32 7.62 12  -14.284     18.9  a    
# nonaxenic     sal      6.84 7.62 12   -9.763     23.4  a    
# axenic        fat     17.21 7.62 12    0.609     33.8  a   
###Control+BA-Spiked Diet Combined####
#We are seeing if BA-spiked diet yields significantly higher SA levels
  #compared to control diet

#start by making a combined control-diet and BA-Spiked-diet data set
df.conba <- df.1 %>%
  filter(diet_type == "benzoic" | diet_type == "control") %>%
  filter (id != "NC2H_1")

#visualize the data; we are only interested in sa_avg for this data set
hist(df.conba$sa_avg)


#make the model, add diet_type as an explanatory variable
conba_glm <- glm(sa_avg ~ microbe_level*tissue*diet_type, family = gaussian(), df.conba)

#check assumptions
plot(conba_glm) #looks approximately normal and homoscedastic; it could be worse

#summary and emmeans posthoc
summary(conba_glm) #no significant interacitons
emmeans_conba <- emmeans(conba_glm, ~ (microbe_level*tissue*diet_type), type = "response")
#output:
# cld(emmeans_conba, Letters = letters, alpha = 0.05) 
# microbe_level tissue diet_type emmean   SE df lower.CL upper.CL .group
# axenic        hem    control      0.0 28.5 24   -58.82     58.8  a    
# axenic        hem    benzoic      0.0 28.5 24   -58.82     58.8  a    
# axenic        fat    control     11.1 28.5 24   -47.75     69.9  a    
# axenic        sal    control     11.5 28.5 24   -47.34     70.3  a    
# axenic        sal    benzoic     21.3 28.5 24   -37.52     80.1  ab   
# nonaxenic     hem    benzoic     21.9 28.5 24   -36.90     80.7  ab   
# nonaxenic     hem    control     31.4 28.5 24   -27.45     90.2  ab   
# nonaxenic     sal    benzoic     37.1 28.5 24   -21.69     96.0  ab   
# axenic        fat    benzoic     60.5 28.5 24     1.66    119.3  ab   
# nonaxenic     sal    control     85.5 28.5 24    26.71    144.4  ab   
# nonaxenic     fat    benzoic     93.7 28.5 24    34.93    152.6  ab   
# nonaxenic     fat    control    157.1 28.5 24    98.24    215.9   b   

##BA spiked diet tissue data####
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



###SA Spiked Diet####
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

###IAA Spiked Diet####
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