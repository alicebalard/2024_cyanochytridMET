##############################
### Gene Ontology analysis ###
##############################

annotationChytrid$GO <- str_extract_all(annotationChytrid$gene_ontology_BLASTX, "GO:\\d+")

library(tidyr)
library(AnnotationDbi)

install.packages("AnnotationDbi")

# create gene universe from all the covered CpGs
gene_universe <- annotationChytrid %>%
  dplyr::filter(lengths(GO)!=0) %>% # rm non existing GO terms
  dplyr::select(c("X.gene_id", "GO")) %>%
  mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, it's from Interproscan after Maker. Sticklebacks (biomart) have 82701 IEA and 63 ISS.
  relocate("GO","go_linkage_type","X.gene_id") %>%
  unnest(GO) %>% # one GO per line (was a list before in this column)
  data.frame()

gene_universe$X.gene_id %>% unique %>% length #5320 genes

# Create gene set collection
goFrame <- GOFrame(gene_universe, organism="Rhizophydium megarrhizum")
goAllFrame <- GOAllFrame(goFrame)
gsc_universe <- GeneSetCollection(goAllFrame, setType = GOCollection())
