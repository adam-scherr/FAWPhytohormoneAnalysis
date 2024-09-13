###Analysis of BA, SA, and IAA Sequestration data, controlled study
####adapted from Second Committee Meeting Code sessions
###Adam Scherr, 11/27/2023
###Version 3, 09/13/2024

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

data.1$diet_type <- factor(data.1$diet_type, 
                           levels = c("control", "benzoic", "salicylic", "indole"))
data.1$tissue <- factor(data.1$tissue,
                        levels = c("sal", "fat", "hem"))

#now make the negative values into zeros
data.1$ba_avg[data.1$ba_avg < 0] <- 0
data.1$sa_avg[data.1$sa_avg < 0] <- 0
data.1$iaa_avg[data.1$iaa_avg < 0] <- 0

data.2 <- data.1 %>%
  select(1:4, 9, 12, 15)

data.sal <- data.2 %>%
  filter(tissue == "sal")



##set up datasets for each phytohormone-control diet comparison
#BA and control
ba.sal <- data.sal %>%
  filter(diet_type != "indole") %>%
  filter(diet_type != "salicylic")

#SA and control
sa.sal <- data.sal %>%
  filter(diet_type != "indole")
 # filter(diet_type != "benzoic")

#IAA and control
iaa.sal <- data.sal %>%
  filter(diet_type != "benzoic") %>%
  filter(diet_type != "salicylic")

#example graph
library(ggplot2)
ba.box <- ggplot(ba.sal, aes(x = diet_type, y=ba_avg, fill = microbe_level))+
  geom_boxplot(width=0.8, color = "black", alpha = 0.82) +
  geom_point(size = 4, alpha = 1, position=position_dodge(width=0.8), aes(color = microbe_level))+
  scale_fill_manual(values = c("#9ecae1", "#3182bd", "#08306b"))+
  scale_color_manual(values = c("#2171b5", "#2171b5", "#2171b5")) + 
  theme_bw()
ba.box

####lets try running ANOVA on this####
#resource: http://www.sthda.com/english/wiki/two-way-anova-test-in-r#google_vignette
#goal: use a two-way anova to show how microbe_level and diet_type affect phytohormone concentration
library(lme4)
library(emmeans)
library(multcomp)
##using BA as first two-way ANOVA example
ba.model <- aov(ba_avg ~ microbe_level*diet_type, data = ba.sal)
summary(ba.model)
plot(ba.model)

ba.model.1 <- aov(ba_avg^0.1 ~ microbe_level*diet_type, data = ba.sal)
plot(ba.model.1)

#check our assumptions
par(mfrow = c(2, 2))
plot(ba.model.scale)

#QQ-plot and variance don't look great
#We have so few data points that it's easy to skew the normality and variance
#So let's go through with it since our small number of points is skewing normality and variance so much
summary(ba.model)
  #p=0.0261 for diet_type. Only diet_type has SIGNIFICANT effect on phytohormone concentration
summary(ba.model.1)
#two-way ANOVA on SA data
sa.model <- aov(sa_avg ~ microbe_level*diet_type, data = sa.sal)

#check assumptions
par(mfrow = c(2,2))
plot(sa.model)
  #also not great equal variance, and not great normality, but once again the small
  #number of points messes with these assumptions easily
sa.model.1 <- aov(sa_avg^0.3 ~ microbe_level*diet_type, data = sa.sal)
plot(sa.model.1)

summary(sa.model)
  #SIGNIFICANT effect of diet_type on phytohormone concentration. p<0.001
summary(sa.model.1)

#two-way ANOVA on IAA data
iaa.model <- aov(iaa_avg ~microbe_level*diet_type, data = iaa.sal)

#check assumptions
par(mfrow(2,2))
plot(iaa.model)
  #variance is better, but normality is not optimal
summary(iaa.model)
  #we have significnace for ALL values!! The things that significantly affect
  #IAA concentrations are microbe_level, diet_type, and the interaction of the two
  #p=0.0040, p<0.001, p=0.0042 respectively
iaa.model.1 <- aov(iaa_avg^0.3 ~microbe_level*diet_type, data = iaa.sal)
plot(iaa.model.1)
summary(iaa.model.1)
###now do post hoc analysis for these two-way ANOVAs
#BA post hoc
TukeyHSD(ba.model, conf.level=0.95) #somehow no significance anywhere. Like. What.
#plot(TukeyHSD(ba.model, conf.level=.95), las = 2)
t.test(ba_avg ~ diet_type, data = ba.sal, var.equal = TRUE) #significance, p=0.024

cld(emmeans(ba.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean   SE df lower.CL upper.CL .group
# nonaxenic     control     12.7 1408  8    -3235     3260  a    
# axenic        control     23.5 1408  8    -3224     3271  a    
# axenic        benzoic   2261.3 1408  8     -986     5509  a    
# nonaxenic     benzoic   5413.1 1408  8     2166     8660  a 

#see if the transformed BA model is any different
TukeyHSD(ba.model.1, conf.level=0.95)
cld(emmeans(ba.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean   SE df lower.CL upper.CL .group
# nonaxenic     control    0.480 0.35  8   -0.327     1.29  a    
# axenic        control    0.944 0.35  8    0.137     1.75  ab   
# axenic        benzoic    2.151 0.35  8    1.345     2.96   b   
# nonaxenic     benzoic    2.266 0.35  8    1.460     3.07   b 
        #it IS  different. We now have significant effect of diet type/kinda microbe_level

#SA post hoc
TukeyHSD(sa.model, conf.level=0.95)
  #all interactions between SA and control/BA diets are significant (p<0.05)

cld(emmeans(sa.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
 # microbe_level diet_type  emmean   SE df lower.CL upper.CL .group
 # axenic        control      11.5 1340 12    -2908     2931  a    
 # axenic        benzoic      21.3 1340 12    -2898     2941  a    
 # nonaxenic     benzoic      37.1 1340 12    -2882     2957  a    
 # nonaxenic     control      85.5 1340 12    -2834     3005  a    
 # nonaxenic     salicylic 11238.1 1340 12     8319    14158   b   
 # axenic        salicylic 13569.8 1340 12    10650    16489   b 

#see if the transformed SA model is any different
TukeyHSD(sa.model.1, conf.level=0.95)
cld(emmeans(sa.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean    SE df lower.CL upper.CL .group
# axenic        control     1.56 0.931 12   -0.472     3.58  a    
# axenic        benzoic     2.42 0.931 12    0.391     4.45  a    
# nonaxenic     benzoic     2.60 0.931 12    0.574     4.63  a    
# nonaxenic     control     2.85 0.931 12    0.818     4.87  a    
# nonaxenic     salicylic  16.17 0.931 12   14.142    18.20   b   
# axenic        salicylic  17.34 0.931 12   15.309    19.36   b  
        #the transformed model is not any different

#IAA post hoc
TukeyHSD(iaa.model, conf.level = 0.95)
  #all interactions significant EXCEPT nonaxenic control and axenic control are the same (p=0.999)

cld(emmeans(iaa.model, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
# microbe_level diet_type   emmean   SE df lower.CL upper.CL .group
# axenic        control       0.00 1156  8    -2665     2665  a    
# nonaxenic     control       6.84 1156  8    -2658     2672  a    
# axenic        indole     7131.06 1156  8     4466     9796   b   
# nonaxenic     indole    16317.88 1156  8    13653    18983    c  

#see if the transformed SA model is any different
TukeyHSD(iaa.model.1, conf.level = 0.95)
cld(emmeans(iaa.model.1, pairwise ~ microbe_level*diet_type, adjust = "tukey"),
    Letters =letters, alpha = 0.05)
  #output:
# microbe_level diet_type emmean    SE df lower.CL upper.CL .group
# axenic        control    0.000 0.658  8   -1.518     1.52  a    
# nonaxenic     control    0.825 0.658  8   -0.693     2.34  a    
# axenic        indole    14.228 0.658  8   12.710    15.75   b   
# nonaxenic     indole    18.297 0.658  8   16.779    19.81    c  
          #transformed model is not any different

#####Trying the suggested linear models instead of ANOVAs####
ba.model <- aov(ba_avg ~ microbe_level*diet_type, data = ba.sal)
summary(ba.model)
plot(ba.model)

hist(ba.sal$ba_avg) #not at all normal looking

#trying a linear model
ba.lm <- lm(ba_avg~microbe_level*diet_type, data=ba.sal)
hist(resid(ba.lm))
car::qqPlot(resid(ba.lm)) #bad qq plot, lets try transforming the data

ba.lm.1 <- lm(ba_avg^0.3~microbe_level*diet_type, data=ba.sal)
hist(resid(ba.lm.1))
car::qqPlot(resid(ba.lm.1)) #that looks so much better

#running the posthoc analysis on the model
car::Anova(ba.lm.1, type = "II")
ba.emmeans <- emmeans::emmeans(ba.lm.1, pairwise ~ c(microbe_level, diet_type),
                               adjust = "tukey")
multcomp::cld(ba.emmeans, Letters = letters)
## fix negative values
ba.sal$ba_avg_norm <- pmax(ba.sal$ba_avg,0)

## check assumptions
hist(ba.sal$ba_avg_norm) # not super normal, straight anova would be tricky

## Use an lm instead
ba.model <- lm(ba_avg_norm ~ microbe_level*diet_type, data = ba.sal)
hist(resid(ba.model))
car::qqPlot(resid(ba.model)) # looks terrible

ba.model <- lm((ba_avg_norm)^(1/3) ~ microbe_level*diet_type, data = ba.sal)
hist(resid(ba.model))
car::qqPlot(resid(ba.model)) # much better

car::Anova(ba.model, type="II")
ba.emmeans <- emmeans::emmeans(ba.model, pairwise ~ c(microbe_level,diet_type), adjust="tukey")
multcomp::cld(ba.emmeans, Letters = letters)
# keep in mind that emmeans does not easily recognize cube root, so for exact values you may need to back transform

#####unspiked data phytohormone comparisons####
unspiked.sal <- data.sal %>%
  filter(diet_type == "control")

##graphs of this data
unspiked.ba.point <- ggplot(unspiked.sal, aes(x = microbe_level, y =ba_avg))+
  geom_point(size=4, aes(color = microbe_level))+
  scale_color_manual(values = c("blue", "darkblue"))
unspiked.ba.point

unspiked.sa.point <- ggplot(unspiked.sal, aes(x = microbe_level, y =sa_avg))+
  geom_point(size=4, aes(color = microbe_level))+
  scale_color_manual(values = c("green", "darkgreen"))
unspiked.sa.point

unspiked.iaa.point <- ggplot(unspiked.sal, aes(x = microbe_level, y=iaa_avg))+
  geom_point(size = 4, aes(color = microbe_level))+
  scale_color_manual(values = c("pink", "maroon"))
unspiked.iaa.point

##running one-way ANOVAs, just to see
#BA
unspiked.ba.anova <- aov(ba_avg~microbe_level, data = unspiked.sal)
par(mfrow=c(2,2))
plot(unspiked.ba.anova)
summary(unspiked.ba.anova) #no significance, p=0.394
TukeyHSD(unspiked.ba.anova, conf.level=0.95) #no sig, p=0.394

unspiked.sa.anova <- aov(sa_avg~microbe_level, data = unspiked.sal)
plot(unspiked.sa.anova)
summary(unspiked.sa.anova) #no sig, p=0.201
TukeyHSD(unspiked.sa.anova, conf.level=0.95) #no sig, p=0.201

unspiked.iaa.anova <- aov(iaa_avg~microbe_level, data = unspiked.sal)
plot(unspiked.iaa.anova)
summary(unspiked.iaa.anova) #no sig, p=0.288
TukeyHSD(unspiked.iaa.anova, conf.level=0.95) #no sig, p=0.288
