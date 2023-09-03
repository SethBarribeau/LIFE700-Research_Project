setwd("~/4th Year Modules/LIFE700/Wimalasiri-Yapa")
library(vegan)

x<-read.csv('Wimalasiri_totalCounts_notEuk_taxLevel1.csv', sep=",")
# still some NAs
x[, 7:75][is.na(x[, 7:75])] <- 0

namesD <- as.data.frame(names(x[7:75]))
names(namesD)<-'sample'

rownames(x) <- x[,1] 
scounts <- x[,c(7:75)]

# Count per millions
scountsPerMillion <- cpm(scounts)
summary(scountsPerMillion)

meta <- read.csv('Wimalasiri-Yapa.metadata.csv', sep= ',', header=T)
meta <- meta[,c(1,11,12,13)]
names(meta)<-c('sample', 'time', 'temperature', 'infection')
table(meta$time)
table(meta$temperature)
table(meta$infection)
Infection<-factor(meta$infection,levels=c("Control","Chikungunya"))
table(Infection)
tempTreat <- factor(meta$temperature)

class(x)
str(x)

shannondiv<-diversity(x[,7:75],index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

#treatment subsets
Ch.18.3<-x[,c(8,35,46,56,66)]
Ch.18.7<-x[,c(62,63,64,65,67,68)]
Ch.28.3<-x[,c(69,70,71,72,73,74)]
Ch.28.7<-x[,c(51,52,53,54,55)]
Ch.32.3<-x[,c(7,9,10,13,24,75)]
Ch.32.7<-x[,c(57,58,59,60,61)]
Co.18.3<-x[,c(44,45,47,48,49,50)]
Co.18.7<-x[,c(26,27,28,29,30)]
Co.28.3<-x[,c(31,32,33,34,36,37)]
Co.28.7<-x[,c(18,19,20,21,22,23,25)]
Co.32.3<-x[,c(38,39,40,41,42,43)]
Co.32.7<-x[,c(11,12,14,15,16,17)]

#shannon each subset
shannondiv<-diversity(Ch.18.3,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Ch.18.7,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Ch.28.3,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Ch.28.7,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Ch.32.3,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Ch.32.7,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Co.18.3,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Co.18.7,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Co.28.3,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Co.28.7,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Co.32.3,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Co.32.7,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

# get right way around for glm
#shannon temperatures
shannondiv.T<-diversity(t(x[,7:75]), index="shannon")
head(shannondiv.T)
length(shannondiv.T)
dfDiv = as.data.frame(shannondiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
metaShort = meta[,c(1,11,12,13)]
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(shannondiv.T~temperature*infection*time, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=shannondiv.T, fill=Infection)) +
  scale_fill_manual(values = c("cornflowerblue","firebrick"))+
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',
position=position_dodge(1), dotsize=0.5, color='black'
pt<- ggplot(dfDiv, aes(x=as.factor(time), y=shannondiv.T, fill=Infection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","firebrick"))+
  ylab("Shannon's Index") +
  xlab("Time") +
  theme_bw() 
pt

#Simpsons temperatures
simpsonsdiv.T<-diversity(t(x[,7:75]), index="simpson")
head(simpsonsdiv.T)
length(simpsonsdiv.T)
dfDiv = as.data.frame(simpsonsdiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
metaShort = meta[,c(1,11,12,13)]
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(simpsonsdiv.T~temperature*infection*time, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=simpsonsdiv.T, fill=Infection)) +
  scale_fill_manual(values = c("cornflowerblue","firebrick"))+
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Simpson's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',
position=position_dodge(1), dotsize=0.5, color='black')
pt<- ggplot(dfDiv, aes(x=as.factor(time), y=simpsonsdiv.T, fill=Infection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","firebrick"))+
  ylab("Simpson's Index") +
  xlab("Time") +
  theme_bw() 
pt

#InvSimpsons temperatures
invsimpsonsdiv.T<-diversity(t(x[,7:75]), index="invsimpson")
head(invsimpsonsdiv.T)
length(invsimpsonsdiv.T)
dfDiv = as.data.frame(invsimpsonsdiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
metaShort = meta[,c(1,11,12,13)]
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(invsimpsonsdiv.T~temperature*infection*time, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=invsimpsonsdiv.T, fill=Infection)) +
  scale_fill_manual(values = c("cornflowerblue","firebrick"))+
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("InvSimpson's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',
position=position_dodge(1), dotsize=0.5, color='black'
pt<- ggplot(dfDiv, aes(x=as.factor(time), y=invsimpsonsdiv.T, fill=Infection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","firebrick"))+
  ylab("InvSimpson's Index") +
  xlab("Time") +
  theme_bw() 
pt



