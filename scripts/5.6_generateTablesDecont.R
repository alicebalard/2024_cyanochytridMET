setwd("~/Documents/Erika-chytridProject/GIT/2024_cyanochytridMET/scripts/")

Z1Z12 <- read.table("../../../Z1Z12.table.tsv", sep="\t", header = T)
In1In12 <- read.table("../../../In1In12.table.tsv", sep="\t", header = T)

tabFungiFamZ1Z12 <- Z1Z12[Z1Z12$bestsumorder_kingdom %in% "Fungi",] %>% 
  count(bestsumorder_family) %>%
  arrange(desc(n))

tabFungiFamIn1In12 <- In1In12[In1In12$bestsumorder_kingdom %in% "Fungi",] %>% 
  count(bestsumorder_family) %>%
  arrange(desc(n))

write.table(tabFungiFamZ1Z12, file = "../figTab/tabFungiFam_Z1Z12_e-25pval.txt", row.names = F)
write.table(tabFungiFamIn1In12, file = "../figTab/tabFungiFam_In1In12_e-25pval.txt", row.names = F)
