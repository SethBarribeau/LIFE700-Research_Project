setwd("~/4th Year Modules/LIFE700")

logFCs<-read.csv('species_logFC.csv', header=T)
species<-factor(logFCs$Species,levels=c("Enterobacter ludwigii","Gardnerella vaginalis","Rothia mucilaginosa","Enterobacter roggenkampii"))
temperature<-factor(logFCs$Temperature,levels=c("18","20","25","28","32","35","36"))

sumlogFCs<-summary(logFCs)
sumlogFCs

plot(logFCs$ï..logFC~species+temperature,data=logFCs)

logFCsaov<-aov(logFCs$ï..logFC~temperature,data=logFCs)
logFCslinear<-lm(logFCs$ï..logFC~temperature,data=logFCs)
summary(logFCslinear)
plot(logFCslinear)

pairwise.t.test(logFCs$ï..logFC,temperature,p.adj="bonf")
TukeyHSD(logFCsaov)
plot(TukeyHSD(logFCsaov))

#example
install.packages("multcomp")
expt1 <- read.csv("species_logFC.csv",header=T)
amod <- aov(expt1$ï..logFC~expt1$Temperature,data=expt1)
summary(amod)
library("multcomp")
tmod <- glht(amod,linfct=mcp(expt1$Temperature="Tukey"))
summary(tmod)
TukeyHSD(amod, conf.level=0.95)
TukeyHSD(amod, conf.level=0.99)
T <- summary(tmod)$test$tstat
as.data.frame(T)
k <- amod$rank
v <- amod$df.residual
pValScheffe <- 1-pf(T**2/(k-1),k-1,v)
as.data.frame(pValScheffe)
contrasts <- rbind(
  "B - A" = c(-1,1,0,0), 
  "C - A" = c(-1,0,1,0), 
  "D - A" = c(-1,0,0,1),
  "C - B" = c(0,-1,1,0),
  "D - B" = c(0,-1,0,1),
  "D - C" = c(0,0,-1,1)
)
contrasts
bmod <- glht(amod,linfct=mcp(X=contrasts))
summary(bmod,test=adjusted("bonferroni"))
summary(bmod,test=adjusted("holm"))
contrasts <- rbind(
  "B - A" = c(-1,1,0,0), 
  "C - A" = c(-1,0,1,0), 
  "D - A" = c(-1,0,0,1)
)
bmod <- glht(amod,linfct=mcp(X=contrasts))
summary(bmod,test=adjusted("bonferroni"))
summary(bmod,test=adjusted("holm"))