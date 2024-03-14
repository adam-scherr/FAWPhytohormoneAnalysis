###Analysis of Phytohormone in Different Tissues
###For Thesis Proposal (and final manuscript, probably)
### Adam Scherr, 2/16/2024
### Version 1

#### Making Separate Datasets for Analysis ####
#start by importing the dataset
df <- read.csv("https://raw.githubusercontent.com/adam-scherr/FAWPhytohormoneAnalysis/main/2023_11_30_CommitteeMeeting/EntSoc_SeqData_2%20-%20tidy_data.csv")

#remove anything unnecessary: 
  #labeled, everything but avg concentrations, unneeded sample qualities, blanks
  #and make all the objects what you want them to be
library(dplyr)
library(ggplot2)
library(lme4)
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

fit.1 <- lm(avg_ba_conc ~ microbe_level + tissue + diet_type, data = df.1)
summary(fit.1)
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

#to tell which family (Gaussian, poison, etc. works best, use the residuals)
  #residuals are formed based on a linear regression (AKA line of best fit) for the data points
  #Residual = observed value - predicted value; essentially, how far off was your line of best fit
  #You will have one residual for each data point; better R-squared means smaller absolute values of residuals 
        #(because negative residual just means the prediction was too high for the true value)
  #Fun fact: the sum of all residuals will equal zero, and the mean of all residuals is 0

#assumptions for a GLM are based on this blog post:
  #https://www.r-bloggers.com/2014/04/checking-glm-model-assumptions-in-r/

hist(log(df.ba$ba_avg))
###GLM for BA Concentration Data
library(lme4)
gaussian_glm <- glm(log(ba_avg) ~ microbe_level + tissue, family = gaussian(),df.ba)
gaussian_glm

poisson_glm <- glm(ba_avg ~ microbe_level + tissue, family = poisson(), df.ba)
warnings()
poisson_glm

summary(gaussian_glm)
summary(poisson_glm) #poisson and logged gaussian are NOT the same

df.ba.n <- df.ba %>%
  mutate(log_ba = log(ba_avg))
logged_gaussian_glm <- glm(log_ba ~ microbe_level + tissue, family = gaussian(), df.ba.n)
summary(logged_gaussian_glm)
summary(gaussian_glm) #yeah, these two are the same

#so we're going to use Gaussian, which is the same as just running a linear model
ba_glm <- glm(log(ba_avg) ~ microbe_level + tissue, family = gaussian(), df.ba)
ba_lm <- lm(log(ba_avg) ~ microbe_level + tissue, df.ba)

#now lets see if these two are the same
summary(ba_glm)
summary(ba_lm)
  #these two are the same. I suppose an lm is all you need because it's Gaussian Distribution

par(mfrow = c(2, 2))
plot(ba_lm)
  #the assumptions look good! The top left graph residuals look pretty well distributed
  #and the top right graph showing normality is a pretty straight line, which means pretty normal distribution

####GLMs for Interaction of microbe_level and tissue for each phytohormone####

#you need the emmeans package for the post-hoc tests (Tukey kramer, I think?)
library(emmeans)

#interaction between them
par(mfrow=c(1,1))
fit.ba <- lm(log(ba_avg) ~ microbe_level * tissue, data = df.ba)
summary(fit.ba)
car::qqPlot(resid(fit.ba))
hist(resid(fit.ba))
car::Anova(fit.ba, type = "III")
summary(aov(fit.ba))

emmeans_ba <- emmeans(fit.ba, ~ microbe_level * tissue, type = "response")
pairs_ba <- pairs(emmeans_ba)

fit.iaa <- lm(iaa_avg ~ microbe_level * tissue, data = df.iaa)
summary(fit.iaa)
car::qqPlot(resid(fit.iaa))
hist(resid(fit.iaa))
car::Anova(fit.iaa, type = "III")

#tukeyHSD(fit.iaa)
emmeans_iaa.1 <- emmeans(fit.iaa, pairwise ~ microbe_level * tissue)
emmeans_iaa.2 <- emmeans(fit.iaa,pairwise ~ c(microbe_level,tissue))
pairs_iaa <- pairs(emmeans_iaa)
summary(emmeans_iaa.1)
summary(fit.iaa)

#trying to do compact letter display
library(multcomp)
#cld(emmeans_iaa$emmeans, adjust="tukey",Letters = letters,alpha=0.05) #compact letter display is what
                                                          #you would use to make the letters over the graph bars

cld(emmeans_iaa.1, Letters = letters)
pwpm(emmeans_iaa.1$emmeans) ## <- THIS IS HELPFUL, the output shows what is significant to what
#pairs(emmeans_iaa.1$emmeans)

ba_interaction_glm <- glm(log(ba_avg) ~ microbe_level + tissue + (microbe_level:tissue),
                       family = gaussian(), df.ba)
summary(ba_interaction_glm)
par(mfrow=c(2,2))
plot(ba_interaction_glm)

#posthoc on ba interaction glm
emmeans(ba_interaction_glm, pairwise ~ (microbe_level), type = "response")
emmeans(ba_interaction_glm, pairwise ~ (tissue), type = "response")
emmeans(ba_interaction_glm, pairwise ~ (microbe_level*tissue), type = "response")
  #we're getting no signifiance for any interactions. Now we know!


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