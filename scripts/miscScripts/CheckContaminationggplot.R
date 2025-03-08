setwd("Documents/Erika-chytridProject/GIT/2024_cyanochytridMET/figTab/")

kingdomBlast <- read.csv("FINALTran_tableDiamondKingdom.csv")
phylumBlast <- read.csv("FINALTran_tableDiamondPhylum.csv")
scoresBlast <- read.table("diamond_1e-5pvalFINALTran.out")

head(phylumBlast)
head(kingdomBlast)
head(scoresBlast)

bitscores = scoresBlast[c(1, 6, 14, 15)]
names(bitscores) = c("id", "percentSimilarity", "evalue", "bitscore")

fungi <- kingdomBlast[kingdomBlast$bestsumorder_kingdom %in% "Fungi",]
Fungi <- phylumBlast[phylumBlast$X_id %in% fungi$X_id,]
rm(fungi)

Fungi <- merge(Fungi, bitscores)
library(scales)
library(ggplot2)

P1 <- ggplot(Fungi, aes(x=percentSimilarity, y=evalue, fill=assembly.reads_cov))+
  geom_point(aes(size=length), shape = 21, alpha = .7)+
  facet_wrap(.~bestsumorder_phylum)+
  scale_fill_gradient(low = "blue", high = "red")+
  scale_y_continuous(labels = scientific) +
  theme_bw()
P1

P2<-ggplot(Fungi, aes(x=percentSimilarity, y=evalue))+
  geom_point(aes(size=length), shape = 21, alpha = .7)+
  facet_wrap(.~bestsumorder_phylum)+
  scale_y_continuous(labels = scientific) +
  theme_bw()

P2

P3 <- ggplot(Fungi, aes(x=percentSimilarity, y=length, 
                        fill=assembly.reads_cov))+
  geom_point(aes(size=length), shape = 21, alpha = .7)+
  facet_wrap(.~bestsumorder_phylum)+
  scale_fill_gradient(low = "blue", high = "red")+
  theme_bw()
P3

p4 <- ggplot(Fungi, aes(x=gc, y=length, fill=assembly.reads_cov))+
  geom_point(aes(size=length), shape = 21, alpha = .7)+
  facet_wrap(.~bestsumorder_phylum)+
  scale_fill_gradient(low = "blue", high = "red")+
  geom_hline(yintercept = 5000)+
  theme_bw()
p4

Fungi$gc30to60 = (Fungi$gc <0.6 & Fungi$gc >0.3)

P5 <- ggplot(Fungi, aes(x=gc, y=evalue, fill=gc30to60))+
  geom_point(aes(size=length), shape = 21, alpha = .9)+
  facet_wrap(.~bestsumorder_phylum)+
  theme_bw()
P5
