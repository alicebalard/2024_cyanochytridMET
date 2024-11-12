## Select transcripts from both transcriptomes that are blastx-ing with fungi
assZ <- read.csv("/scratch/alicebalard/outData/diamondBlastX/assemblyZ_diamond_1e-5pval.out", sep = "\t", header = F)

assIn <- read.csv("/scratch/alicebalard/outData/diamondBlastX/assemblyIn_diamond_1e-5pval.out", sep = "\t", header = F)

names(assZ) <- c("qseqid", "staxids", "bitscore", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sscinames", "sskingdoms", "skingdoms", "sphylums")

names(assIn) <- c("qseqid", "staxids", "bitscore", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sscinames", "sskingdoms", "skingdoms", "sphylums")

assZ$skingdoms %>% table

assIn$skingdoms %>% table 

## % of fungal reads in each assembly
sum(assZ$skingdoms %in% "Fungi") # 37793
sum(assZ$skingdoms %in% "Fungi")/nrow(assZ)*100 # 6.97%

sum(assIn$skingdoms %in% "Fungi") # 37842
sum(assIn$skingdoms %in% "Fungi")/nrow(assIn)*100 # 4.01%

assZ_Fungi <- assZ[assZ$skingdoms %in% "Fungi",]
nrow(assZ_Fungi) # check

write.csv(assZ_Fungi["qseqid"], "/scratch/alicebalard/outData/diamondBlastX/assZ_Fungi_transcripts", quote=F, row.names=F)

assIn_Fungi <- assIn[assIn$skingdoms %in% "Fungi",]
nrow(assIn_Fungi) # check

write.csv(assIn_Fungi["qseqid"], "/scratch/alicebalard/outData/diamondBlastX/assIn_Fungi_transcripts", quote=F, row.names=F)
