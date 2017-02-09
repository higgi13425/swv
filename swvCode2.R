#load libraries
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(prediction)
library(margins)
library(ggplot2)
library(broom)
library(intubate)
library(magrittr)
library(ggbeeswarm)
library(corrr)

#clear environment and set working directory
rm(list=ls())

#read in excel file
df<- read_excel("RTNBS30_SWV_tidy.xlsx", sheet=1)

#convert swv measurements to numeric
df$swv<-as.numeric(df$swv)

#rename weight as wall thickness
colnames(df)[2]<- "wallthick"

#reset perf to a 0/1 dichotomous variable
df$perf = df$perf-1

#create clock12 variable with 12 o'clock=1, others (9,3) = 0
df$clock12<-0
df$clock12[df$clock>9]<-1

#filter out those that are perf=1 or swv=NA
df <- df %>% 
  filter(perf==0) %>% 
  filter(!is.na(swv))

# test whether there are significant differences by clock position - first a table
df %>%  group_by(clock) %>% summarise(mean_swv = mean(swv, na.rm=T), sd_swv=sd(swv, na.rm=T))

#then anolysis of variance
p<- aov(swv~clock, data=df)
summary(p)
#clock position definitely an issue

#filter for those that are clock12==1
df12 <- df %>% 
  filter(clock12==1)

#rerun mean, SD and anova only on 12 o'clock measurements
# test whether there are significant differences by clock position - first a table
df %>%  group_by(weeks) %>% summarise(mean_swv = mean(swv, na.rm=T), sd_swv=sd(swv, na.rm=T))
#then anolysis of variance
p<- aov(swv~weeks, data=df)
summary(p)


#linear modeling of df12
#linear modeling
lfit <- lm(swv ~ strain + category  + trial + Animal + wallthick, data=df12) # .13
lfit <- lm(swv ~ strain + category + Animal + wallthick, data=df12) # .1324
lfit <- lm(swv ~ strain + category + Animal, data=df12) # .1335
lfit <- lm(swv ~ strain * category  + Animal , data=df12) # .1346
lfit <- lm(swv ~ strain * weeks  + Animal , data=df12) # .1401
lfit <- lm(swv ~ strain + weeks + Animal , data=df12) # .135
lfit <- lm(swv ~ strain * weeks * wallthick, data=df12) # .1282
lfit <- lm(swv ~ strain * weeks *Animal , data=df12) # .1637

summary(lfit)
m <- margins(lfit)
summary(m)
plot(m[[1]])
#cplot(lfit, x="Weeks", se.type="shade")
#persp(lfit, xvar= "strain", yvar = "Weeks")
#image(lfit, xvar= "strain", yvar = "Weeks")

#mixed modeling with df12
mfit <- lmer(swv ~ strain * weeks + (trial | Animal/trial), data=df12)
summary(mfit)
coefs <- data.frame(coef(summary(mfit)))
coefs$p.z <- 2* (1-pnorm(abs(coefs$t.value)))
coefs


################################
# ttests
#3 vs 6, using only 12 oclock, all strains
df12 %>% filter(weeks>2) %>% group_by(weeks) %>% mutate(as.factor(weeks)) %>% 
ntbt(t.test, swv ~ weeks)

#0 vs 1, using only 12 oclock, , all strains
df12 %>% filter(weeks<3) %>% group_by(weeks) %>% mutate(as.factor(weeks)) %>% 
  ntbt(t.test, swv ~ weeks)

#0 vs 6, using only 12 oclock, all strains
df12 %>% filter(weeks<1 | weeks>4) %>% group_by(weeks) %>% mutate(as.factor(weeks)) %>% 
  ntbt(t.test, swv ~ weeks)

#1 vs 3, using only 12 oclock, all strains
df12 %>% filter(weeks==1 | weeks==3)  %>% group_by(weeks) %>% mutate(as.factor(weeks)) %>% 
  ntbt(t.test, swv ~ weeks)

#0 vs 3, using only 12 oclock, all strains
df12 %>% filter(weeks==0 | weeks==3)  %>% group_by(weeks) %>% mutate(as.factor(weeks)) %>% 
  ntbt(t.test, swv ~ weeks)

#1 vs 6, using only 12 oclock, all strains
df12 %>% filter(weeks==1 | weeks==6)  %>% group_by(weeks) %>% mutate(as.factor(weeks)) %>% 
  ntbt(t.test, swv ~ weeks)

#filter out mutiple trials in order to analyze biomarkers. Use mean value for SWV
df12sum <- df12 %>%  group_by(Animal) %>% 
  mutate(meanswv= mean(swv)) %>% 
  filter(trial==1)

#means by strain
df12 %>% group_by(weeks, strain) %>% summarize(meanswv=mean(swv))

#Beeswarm-Box plots swv by Weeks, select strain
df12graph <-df12sum %>% filter(strain==10)
g <- ggplot(df12graph, aes(factor(weeks), swv))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("Shear Wave Speed by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="Shear Wave Speed (m/sec") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots Collagen by Weeks
g <- ggplot(df12sum, aes(factor(weeks),  col1A1))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("Collagen by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="Collagen") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots IGF by Weeks
g <- ggplot(df12sum, aes(factor(weeks),  IGF))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("IGF by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="IGF") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots IL1 by Weeks
g <- ggplot(df12sum, aes(factor(weeks),  IL1))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("IL1 by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="IL1") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots CTGF by Weeks
g <- ggplot(df12sum, aes(factor(weeks),  CTGF))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("CTGF by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="CTGF") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots TGFb by Week
g <- ggplot(df12sum, aes(factor(weeks),  TGF))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("TGFb by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="TGFb") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots Akt3 by Weeks
g <- ggplot(df12sum, aes(factor(weeks),  AKT3))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("Akt3 by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="Akt3") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots Axl by Weeks
g <- ggplot(df12sum, aes(factor(weeks),  AXL))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("Axl by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="AXL") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)

#Beeswarm-Box plots aSMA by Weeks
g <- ggplot(df12sum, aes(factor(weeks),  aSMA_prot))
g +geom_boxplot() + geom_quasirandom(varwidth = T, cex=1.5) +   ggtitle("aSMA by Weeks of TNBS Treatment") +
  labs(x="Weeks of TNBS Treatments") + labs(y="aSMA") +
  theme(plot.title = element_text(family = "Arial", color = "black", face="bold", size=20, hjust=0.5)) +
  theme(axis.title = element_text(family = "Arial", color = "black", face="bold", size=16)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=15, color="red", alpha= 0.5)



###################################
#test split and model slopes

td <-df12 %>% 
  group_by(weeks) %>% 
  do(tidy(lm(swv ~ strain, .))) %>% 
  filter(term=="strain") %>% 
  select(c(weeks, estimate, std.error)) 
td

td2 <-df12 %>% 
  filter(!is.na(swv)) %>% 
  group_by(weeks, Animal) %>% 
  do(tidy(lm(swv ~ strain, .))) %>% 
  filter(term=="strain") %>% 
  select(c(weeks, Animal, term, estimate)) %>% 
  filter(weeks>2)
td2


ggplot(td, aes(estimate, factor(weeks), color=weeks)) +
  geom_point(aes(size=6)) +
  scale_x_continuous(name= "Slope of Strain as predictor of SWV", limits= c(-0.0375,0.060)) +
  geom_errorbarh(xmin=td$estimate-1.96*td$std.error, xmax=td$estimate+1.96*td$std.error)

