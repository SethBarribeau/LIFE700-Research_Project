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
meta <- meta[,c(1,11,12,13,15)]
names(meta)<-c('sample', 'time', 'temperature', 'infection','treatment')
table(meta$time)
table(meta$temperature)
infection<-factor(meta$infection,levels=c("Control","Chikungunya"))
table(infection)
tempTreat <- factor(meta$temperature)

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

#transpose so samples are rows
t_scounts <- as.data.frame(t(scounts))
#rarefy data
min_depth = min(colSums(scounts))
t_scounts_rarefied <- as.data.frame(round(rrarefy(t_scounts, min_depth)))
#sqr root transform rarefied table and determine best method calculating distance matrix
sqrt_t_scounts_rarefied = sqrt(t_scounts_rarefied)
rank.tscounts <- rankindex(as.matrix(sqrt_t_scounts_rarefied), t_scounts_rarefied, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")
print(paste("The highest rank was given by the", names(sort(rank.tscounts, decreasing = TRUE)[1]), "method."))
scounts_dist = as.matrix((vegdist(t_scounts_rarefied, "bray")))

#NMDS
NMDS = metaMDS(scounts_dist)
#build a data frame with NMDS coordinates and metadata
#but first factor time and temp
meta$time=factor(meta$time,labels=c("3","7"))
meta$temperature=factor(meta$temperature,labels=c("18","28","32"))
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, infection=infection, time=meta$time, temperature=meta$temperature)
head(NMDS)
#how do samples fall in ordination space?
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=infection)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "NMDS Plot")
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=time)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=temperature)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")

#ANOSIM
#infection
anosim_infection = anosim(scounts_dist, infection)
anosim_infection # take a look at results
summary(anosim_infection)
plot(anosim_infection)
#time
anosim_time = anosim(scounts_dist, meta$time)
anosim_time # take a look at results
summary(anosim_time)
plot(anosim_time)
#temperature
anosim_temperature = anosim(scounts_dist, meta$temperature)
anosim_temperature # take a look at results
summary(anosim_temperature)
plot(anosim_temperature)
#all of these showed significance which I think means taht our communities were different frome achother

#permANOVA
#infection
adonis_infection = adonis(scounts_dist ~ infection, meta)
adonis_infection
summary(adonis_infection)
plot(adonis_infection)
#no note of the significance for infection
#time
adonis_time = adonis(scounts_dist ~ time, meta)
adonis_time 
summary(adonis_time)
plot(adonis_time)
#significant for time
#temperature
adonis_temperature = adonis(scounts_dist ~ temperature, meta)
adonis_temperature 
summary(adonis_temperature)
plot(adonis_temperature)
#significant for temperature