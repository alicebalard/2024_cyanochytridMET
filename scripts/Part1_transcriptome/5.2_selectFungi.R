## Select transcripts from both transcriptomes that are blastx-ing with fungi

library(dplyr)

assZ <- read.csv("/scratch/alicebalard/outData/diamondBlastX/assemblyZ_diamondNR_1e-5pval.out", sep = "\t", header = F)

assInCo <- read.csv("/scratch/alicebalard/outData/diamondBlastX/assemblyIn_cocult_diamondNR_1e-5pval.out", sep = "\t", header = F)

names(assZ) <- c("qseqid", "staxids", "bitscore", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sscinames", "sskingdoms", "skingdoms", "sphylums")

names(assInCo) <- c("qseqid", "staxids", "bitscore", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sscinames", "sskingdoms", "skingdoms", "sphylums")

assZ$skingdoms %>% table

assInCo$skingdoms %>% table 

## % of fungal reads in each assembly
sum(assZ$skingdoms %in% "Fungi") # 43644
sum(assZ$skingdoms %in% "Fungi")/nrow(assZ)*100 # 7.69%

sum(assInCo$skingdoms %in% "Fungi") # 36112
sum(assInCo$skingdoms %in% "Fungi")/nrow(assInCo)*100 # 4.35%

assZ_Fungi <- assZ[assZ$skingdoms %in% "Fungi",]
nrow(assZ_Fungi) # 43644

write.csv(assZ_Fungi["qseqid"], "/scratch/alicebalard/outData/diamondBlastX/assZ_nr_Fungi_transcripts", quote=F, row.names=F)

assInCo_Fungi <- assInCo[assInCo$skingdoms %in% "Fungi",]
nrow(assInCo_Fungi) # 36112

write.csv(assInCo_Fungi["qseqid"], "/scratch/alicebalard/outData/diamondBlastX/assInCo_nr_Fungi_transcripts", quote=F, row.names=F)
