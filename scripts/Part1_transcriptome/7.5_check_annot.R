oldAnnot <- read.csv("data/allFungiTrinot_simplified.tsv", sep = "\t")

newAnnot <- read.csv("data/assemblyMergedFungi_simplified.tsv", sep = "\t")

head(oldAnnot,1)




test <- newAnnot$sprot_Top_BLASTX_hit[1:10]

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

oldAnnot <- extractBlastxTaxa(oldAnnot)
newAnnot <- extractBlastxTaxa(newAnnot)

table(oldAnnot$blastxKingdom)
table(newAnnot$blastxKingdom)

table(oldAnnot$blastxKingdom, oldAnnot$blastxPhylum)
table(newAnnot$blastxKingdom, newAnnot$blastxPhylum)

table(newAnnot$blastxPhylum[
               newAnnot$blastxKingdom %in% "Eukaryota"])

table(newAnnot$blastxPhylum[
               newAnnot$blastxKingdom %in% "Bacteria"])
