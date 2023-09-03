setwd("~/4th Year Modules/LIFE700/Ferreira")

library(vegan)

x<-read.csv('totalCounts_notEuk_taxLevel1_integer_zeros.csv', sep=",")
# still some NAs
x[,6:149][is.na(x[,6:149])] <- 0
namesD <- as.data.frame(names(x[6:149]))
names(namesD)<-'sample'
rownames(x) <- x[,1] 
scounts <- x[,c(6:149)]

meta<-read.csv('Ferreira.metadata.csv',sep=',',header=T)
meta <- meta[,c(1,11,12,14)]
names(meta)<-c('sample', 'temperature', 'infection', 'time')
table(meta$time)
# error here one 28h sample is this correct?
meta$time[meta$time == 28] <- 48
tempTreat <- factor(meta$temperature)

class(x)
str(x)

shannondiv<-diversity(x[,6:149],index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

#treatment subsets
Z.24.20<-x[,c(20,21,23,24,25,26,27,28,29,30,31,32)]
Z.24.28<-x[,c(47,48,49,50,51,52,53,118,120,121,122,123)]
Z.24.36<-x[,c(73,74,75,76,78,79,80,81,82,87,98,109)]
Z.48.20<-x[,c(34,35,36,37,38,39,40,41,42,43,45,46)]
Z.48.28<-x[,c(70,71,72,124,125,126,127,128,129,131,132,133)]
Z.48.36<-x[,c(11,22,33,44,64,77,119,130,136,147,148,149)]
C.24.20<-x[,c(83,84,85,86,88,89,90,91,92,93,94,95)]
C.24.28<-x[,c(110,111,112,113,114,115,116,117,134,135,137,138)]
C.24.36<-x[,c(6,58,59,60,61,62,63,65,66,67,68,69)]
C.48.20<-x[,c(96,97,99,100,101,102,103,104,105,106,107,108)]
C.48.28<-x[,c(54,55,56,57,139,140,141,142,143,144,145,146)]
C.48.36<-x[,c(7,8,9,10,12,13,14,15,16,17,18,19)]

#shannon each subset
shannondiv<-diversity(Z.24.20,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Z.24.28,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Z.24.36,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Z.48.20,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Z.48.28,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(Z.48.36,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(C.24.20,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(C.24.28,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(C.24.36,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(C.48.20,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(C.48.28,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(C.48.36,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

# get right way around for glm
#shannon temperatures
shannondiv.T<-diversity(t(x[,6:149]), index="shannon")
head(shannondiv.T)
length(shannondiv.T)
dfDiv = as.data.frame(shannondiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
metaShort = meta[,c(1,11,12,14:17)]
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(shannondiv.T~temperature*infection*time, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=shannondiv.T, fill=infection)) +
  scale_fill_manual(values = c("cornflowerblue", "firebrick"))+
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',
position=position_dodge(1), dotsize=0.5, color='black'
pt<- ggplot(dfDiv, aes(x=as.factor(time), y=shannondiv.T, fill=infection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue", "firebrick"))+
  ylab("Shannon's Index") +
  xlab("Time") +
  theme_bw() 
pt

#Simpsons temperatures
simpsonsdiv.T<-diversity(t(x[,6:149]), index="simpson")
head(simpsonsdiv.T)
length(simpsonsdiv.T)
dfDiv = as.data.frame(simpsonsdiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
metaShort = meta[,c(1,11,12,14:17)]
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(simpsonsdiv.T~temperature*infection*time, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=simpsonsdiv.T, fill=infection)) +
  scale_fill_manual(values = c("cornflowerblue", "firebrick"))+
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Simpson's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',
position=position_dodge(1), dotsize=0.5, color='black')
pt<- ggplot(dfDiv, aes(x=as.factor(time), y=simpsonsdiv.T, fill=infection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue", "firebrick"))+
  ylab("Simpson's Index") +
  xlab("Time") +
  theme_bw() 
pt

#InvSimpsons temperatures
invsimpsonsdiv.T<-diversity(t(x[,6:149]), index="invsimpson")
head(invsimpsonsdiv.T)
length(invsimpsonsdiv.T)
dfDiv = as.data.frame(invsimpsonsdiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
metaShort = meta[,c(1,11,12,14:17)]
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(invsimpsonsdiv.T~temperature*infection*time, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=invsimpsonsdiv.T, fill=infection)) +
  scale_fill_manual(values = c("cornflowerblue", "firebrick"))+
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("InvSimpson's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',
position=position_dodge(1), dotsize=0.5, color='black'
pt<- ggplot(dfDiv, aes(x=as.factor(time), y=invsimpsonsdiv.T, fill=infection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue", "firebrick"))+
  ylab("InvSimpson's Index") +
  xlab("Time") +
  theme_bw() 
pt



