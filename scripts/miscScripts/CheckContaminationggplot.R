setwd("Documents/Erika-chytridProject/GIT/2024_cyanochytridMET/figTab/")

kingdomBlast <- read.csv("FINALTran_tableDiamondKingdom.csv")
phylumBlast <- read.csv("FINALTran_tableDiamondPhylum.csv")
scoresBlast <- read.table("diamond_1e-5pvalFINALTran.out")

head(phylumBlast)
head(kingdomBlast)
head(scoresBlast)

scoresBlast

bitscores = scoresBlast[c(1, 6, 14, 15)]
names(bitscores) = c("id", "percentSimilarity", "evalue", "bitscore")

fungi <- kingdomBlast[kingdomBlast$bestsumorder_kingdom %in% "Fungi",]
Fungi <- phylumBlast[phylumBlast$X_id %in% fungi$X_id,]
rm(fungi)


Fungi <- merge(Fungi, bitscores)
library(scales)
library(ggplot2)

ggplot(Fungi, aes(x=percentSimilarity, y=evalue, fill=assembly.reads_cov))+
  geom_point(aes(size=length), shape = 21, alpha = .7)+
  facet_wrap(.~bestsumorder_phylum)+
  scale_fill_gradient(low = "blue", high = "red")+
  scale_y_continuous(labels = scientific) +
  theme_bw()

ggplot(Fungi, aes(x=percentSimilarity, y=evalue))+
  geom_point(aes(size=length), shape = 21, alpha = .7)+
  facet_wrap(.~bestsumorder_phylum)+
  scale_y_continuous(labels = scientific) +
  theme_bw()ggplot(Fungi, aes(x=percentSimilarity, y=evalue))+
  geom_point(aes(size=assembly.reads_cov), shape = 21, alpha = .7)+
  facet_wrap(.~bestsumorder_phylum)+
  scale_y_continuous(labels = scientific) +
  theme_bw()

  

