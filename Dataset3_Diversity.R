setwd("~/4th Year Modules/LIFE700/Murrieta")

library(vegan)

x<-read.csv('Murrieta_nrnt_notEuk.csv', sep=",")
# still some NAs
x[, 6:186,][is.na(x[, 6:186])] <- 0
namesD <- as.data.frame(names(x[6:186]))
names(namesD)<-'sample'
rownames(x) <- x[,1] 
scounts <- x[,c(6:186)]

meta <- read.csv('Murrieta.metadata.csv', sep= ',', header=T)
meta <- meta[,c(1,12)]
names(meta)<-c('sample', 'temperature')
table(meta$temperature)
tempTreat <- factor(meta$temperature)

#remove 25-35 treatment
meta <- meta[,c(1,2)]
names(meta)<-c('sample', 'temperature')
table(meta$temperature)
meta2=meta[meta$temperature != '25-35',]
names(meta2)<-c('sample', 'temperature')
table(meta2$temperature)
meta2$temperature=factor(meta2$temperature)
keep = scounts[, names(scounts) %in% meta2$sample] # here we drop all the columns that aren't in our new metadata

class(x)
str(x)

shannondiv<-diversity(x[,6:186],index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

#treatment subsets
T25<-scounts[,c(4,15,26,27,28,39,50,61,72,83,94,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,120,131,142,153,164,175)]
T25.35<-scounts[,c(118,119,121,122,123,124,125,126,127,128,129,130,132,133,134,135,136,137,138,139,140,141,143,144,145,146,147,148,149,150,151,152,154,155)]
T28<-scounts[,c(63,64,65,66,67,68,69,70,71,73,74,75,76,77,78,79,80,81,82,84,85,86,87,88,89,90,91,92,93,95,96,97,98,99,100,101,102)]
T32<-scounts[,c(21,22,23,24,25,29,30,31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,62)]
T35<-x[,c(1,2,3,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,156,157,158,159,160,161,162,163,165,166,167,168,169,170,172,173,174,176,177,178,179,180,181,182)]
#problem with T35

#shannon each subset
shannondiv<-diversity(T25,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(T25.35,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(T28,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(T32,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

shannondiv<-diversity(T35,index="shannon")
head(shannondiv)
par(mfrow = c(1, 1))
hist(shannondiv)

# get right way around for glm
#shannon temperatures
shannondiv.T<-diversity(t(keep), index="shannon")
head(shannondiv.T)
length(shannondiv.T)
dfDiv = as.data.frame(shannondiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(shannondiv.T~temperature, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=shannondiv.T)) +
  geom_boxplot(fill="firebrick")+
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',

#Simpsons temperatures
simpsonsdiv.T<-diversity(t(x[,6:186,]), index="simpson")
head(simpsonsdiv.T)
length(simpsonsdiv.T)
dfDiv = as.data.frame(simpsonsdiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(simpsonsdiv.T~temperature, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(meta$temperature), y=simpsonsdiv.T)) +
  geom_boxplot(fill="firebrick")+
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Simpson's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',

#InvSimpsons temperatures
invsimpsonsdiv.T<-diversity(t(x[,6:186]), index="invsimpson")
head(invsimpsonsdiv.T)
length(invsimpsonsdiv.T)
dfDiv = as.data.frame(invsimpsonsdiv.T)
dfDiv$sample=rownames(dfDiv)
dfDiv=dfDiv[,c(2,1)]
head(dfDiv)
dfDiv = merge(dfDiv, meta, by='sample')
testGlm <- glm(invsimpsonsdiv.T~temperature, data=dfDiv)
anova(testGlm, test='F')
# Strong temp fx, weak infection fx (but sig), strong time fx, infection x time interaction.
table(dfDiv$temperature)
library(ggplot2)
p<- ggplot(dfDiv, aes(x=as.factor(temperature), y=invsimpsonsdiv.T)) +
  geom_boxplot(fill="firebrick")+
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("InvSimpson's Index") +
  xlab("Temperature") +
  theme_bw() 
p  # + geom_dotplot(binaxis='y', stackdir='center',