
## Load annotation
newAnnot <- read.csv("/scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_simplified_GOKegg.tsv", sep = "\t")

## Remove everything before "Full"; Extract the first word after the first "^"; rm ;
extractBlastxTaxa <- function(x){
    ## Rm everything before "Full
    data = sub(".*?(Full=.*)", "\\1",
               x$sprot_Top_BLASTX_hit)
    ## Extract everything after the first "^"
    data = sub("^[^\\^]*\\^(.*)", "\\1", data)
    ## Extract everything before the first ";"
    x$blastxKingdom = sub("(.*?);.*", "\\1", data)
    ## Extract everything after the first ";" and before the second ";"
    x$blastxPhylum = sub("^[^;]*;\\s*([^; ]+).*", "\\1", data)
    #sub(" ", "", sub("^[^;]*;([^;]*);.*", "\\1", data))
    return(x)
}

newAnnot <- extractBlastxTaxa(newAnnot)

table(newAnnot$blastxKingdom)
##      .   Archaea  Bacteria Eukaryota   Viruses 
##  48688       192      5801     61449       165 

## NB duplicated isoforms
table(newAnnot$blastxKingdom, newAnnot$blastxPhylum)

nrow(newAnnot[newAnnot$blastxKingdom %in% "Eukaryota",])/
    nrow(newAnnot) *100
## 52.8% is annotated with Eukaryota

eukTransc <- newAnnot[newAnnot$blastxKingdom %in% "Eukaryota","transcript_id"]

nrow(newAnnot[newAnnot$blastxPhylum %in% "Fungi",])/
    nrow(newAnnot) *100
## 25% is annotated with Fungi

eukTransc <- newAnnot[newAnnot$blastxKingdom %in% "Eukaryota","transcript_id"]

length(eukTransc);length(unique(eukTransc))

fungTransc <- newAnnot[newAnnot$blastxPhylum %in% "Fungi","transcript_id"]

length(fungTransc);length(unique(fungTransc))

#################################################
## Subset our transcriptome for these transcripts
library("seqinr")

transc <- read.fasta("/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta")

transc1 <- transc[names(transc) %in% eukTransc]

write.fasta(transc1,
            names(transc1),
            file.out= "/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta", nbchar=50000)
## NB to upper case in bash:
## awk '/^>/ {print; next} {print toupper($0)}' Trinity_eukaryoteHits.fasta > temp
## mv temp Trinity_eukaryoteHits.fasta

write.fasta(transc[names(transc) %in% fungTransc],
            names(transc[names(transc) %in% fungTransc]),
            file.out= "/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_fungiHits.fasta", nbchar=50000)

###############################################
## Subset gene-trans-map for the Eukaryotic one
gtm <- read.table("/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta.gene_trans_map") 

write.table(gtm[gtm$V2 %in% eukTransc,], "/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.gene_trans_map", row.names=F, quote=F, col.names = F)

#################################################
## Subset annotation file for the Eukaryotic hits
write.table(newAnnot[newAnnot$blastxKingdom %in% "Eukaryota",], "/scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_filterEuk_simplified_GOKegg.tsv", sep = "\t", quote = F, row.names=F, col.names=T)
