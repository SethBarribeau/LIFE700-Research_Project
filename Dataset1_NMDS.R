setwd("~/4th Year Modules/LIFE700/Ferreira")

library(vegan)

x<-read.csv('totalCounts_notEuk_taxLevel1_integer_zeros.csv', sep=",")
# still some NAs
x[,6:149][is.na(x[,6:149])] <- 0
namesD <- as.data.frame(names(x[6:149]))
names(namesD)<-'sample'
rownames(x) <- x[,1] 
scounts <- x[,c(6:149)]

#subset without 28C to replace scounts in future code
x2030<-x[,c(20,21,23,24,25,26,27,28,29,30,31,32,73,74,75,76,78,79,80,81,82,87,98,109,34,35,36,37,38,39,40,41,42,43,45,46,11,22,33,44,64,77,119,130,136,147,148,149,83,84,85,86,88,89,90,91,92,93,94,95,6,58,59,60,61,62,63,65,66,67,68,69,96,97,99,100,101,102,103,104,105,106,107,108,7,8,9,10,12,13,14,15,16,17,18,19)]
#subsets of each treatment
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
#subsets without time
Z.20<-x[,c(20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,45,46)]
Z.28<-x[,c(47,48,49,50,51,52,53,118,120,121,122,123,70,71,72,124,125,126,127,128,129,131,132,133)]
Z.36<-x[,c(73,74,75,76,78,79,80,81,82,87,98,109,11,22,33,44,64,77,119,130,136,147,148,149)]
C.20<-x[,c(83,84,85,86,88,89,90,91,92,93,94,95,96,97,99,100,101,102,103,104,105,106,107,108)]
C.28<-x[,c(110,111,112,113,114,115,116,117,134,135,137,138,54,55,56,57,139,140,141,142,143,144,145,146)]
C.36<-x[,c(6,58,59,60,61,62,63,65,66,67,68,69,7,8,9,10,12,13,14,15,16,17,18,19)]
#Subsets of just infection type
zika<-x[,c(20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,45,46,47,48,49,50,51,52,53,118,120,121,122,123,70,71,72,124,125,126,127,128,129,131,132,133,73,74,75,76,78,79,80,81,82,87,98,109,11,22,33,44,64,77,119,130,136,147,148,149)]
control<-x[,c(83,84,85,86,88,89,90,91,92,93,94,95,96,97,99,100,101,102,103,104,105,106,107,108,110,111,112,113,114,115,116,117,134,135,137,138,54,55,56,57,139,140,141,142,143,144,145,146,6,58,59,60,61,62,63,65,66,67,68,69,7,8,9,10,12,13,14,15,16,17,18,19)]

meta<-read.csv('Ferreira.metadata3 (1).csv',sep=',',header=T)
meta <- meta[,c(1,11,12,14,16,17)]
names(meta)<-c('sample','temperature','infection','time','treatment','inf-temp')
table(meta$time)
# error here one 28h sample is this correct?
meta$time[meta$time == 28] <- 48
tempTreat <- factor(meta$temperature)

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
meta$time=factor(meta$time,labels=c("24","48"))
meta$temperature=factor(meta$temperature,labels=c("20","28","36"))
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, infection=meta$infection, time=meta$time, temperature=meta$temperature)
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
anosim_infection = anosim(scounts_dist, meta$infection)
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