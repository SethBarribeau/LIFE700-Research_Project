setwd("~/4th Year Modules/LIFE700/Murrieta")

x<-read.csv('Murrieta_nrnt_notEuk.csv', sep=",")
# still some NAs
x[, 6:186][is.na(x[, 6:186])] <- 0
namesD <- as.data.frame(names(x[6:186]))
names(namesD)<-'sample'
rownames(x) <- x[,1] 
scounts <- x[,c(6:186)]

meta<- read.csv('Murrieta.metadata1.csv', sep= ',', header=T)
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

#treatment subsets
T25<-scounts[,c(4,15,26,27,39,50,61,72,83,94,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,120,131,142,153,164,175)]
T25.35<-scounts[,c(118,119,121,122,123,124,125,126,127,128,129,130,132,133,134,135,136,137,138,139,140,141,143,144,145,146,147,148,149,150,151,152,154,155)]
T28<-scounts[,c(63,64,65,66,67,68,69,70,71,73,74,75,76,77,78,79,80,81,82,84,85,86,87,88,89,90,91,92,93,95,96,97,98,99,100,101,102)]
T32<-scounts[,c(21,22,23,24,25,29,30,31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,62)]
T35<-x[,c(1,2,3,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,156,157,158,159,160,161,162,163,165,166,167,168,169,170,172,173,174,176,177,178,179,180,181,182)]
#problem with T35

library(vegan)

#transpose so samples are rows
t_scounts <- as.data.frame(t(keep))
#rarefy data
min_depth = min(colSums(keep))
t_scounts_rarefied <- as.data.frame(round(rrarefy(t_scounts, min_depth)))
#sqr root transform rarefied table and determine best method calculating distance matrix
sqrt_t_scounts_rarefied = sqrt(t_scounts_rarefied)
rank.tscounts <- rankindex(as.matrix(sqrt_t_scounts_rarefied), t_scounts_rarefied, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")
print(paste("The highest rank was given by the", names(sort(rank.tscounts, decreasing = TRUE)[1]), "method."))
scounts_dist = as.matrix((vegdist(t_scounts_rarefied, "horn")))

#NMDS
NMDS = metaMDS(scounts_dist)
#build a data frame with NMDS coordinates and metadata
#but first factor time and temp
meta2$temperature=factor(meta2$temperature,labels=c("25","28","32","35"))
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, temperature=meta2$temperature)
head(NMDS)
#how do samples fall in ordination space?
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=temperature)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")

#ANOSIM
#temperature
anosim_temperature = anosim(scounts_dist, meta2$temperature)
anosim_temperature # take a look at results
summary(anosim_temperature)
plot(anosim_temperature)
#all of these showed significance which I think means taht our communities were different frome achother

#permANOVA
#temperature
adonis_temperature = adonis(scounts_dist ~ temperature, meta)
adonis_temperature 
summary(adonis_temperature)
plot(adonis_temperature)
#significant for temperature