cd scripts/Part1_transcriptome

FASTA=../../gitignore/Trinity_eukaryoteHits.rmConta.fasta
ANNOT=../../gitignore/assemblyMergedFungi_filterEuk_simplified_GOKegg.tsv
OUT=../../gitignore/Rmegarrhizum_ChyKol2008_denovo_transcriptome_annotation.tsv

python3 filter_annotation_to_fasta.py --fasta $FASTA --annot $ANNOT --out $OUT

##== FASTA (decontaminated transcriptome) ==
##  transcripts (headers) :   16909
##  genes (unique)        :    3894
##== Annotation BEFORE filtering ==
##  rows                  :   56441
##  transcripts (unique)  :   20356
##  genes (unique)        :    4860
##== Annotation AFTER filtering (written to --out) ==
##  rows                  :   51340
##  transcripts (unique)  :   16909
##  genes (unique)        :    3894
##== Sanity checks ==
##  annotation transcripts removed as not-in-FASTA :    3447
##  FASTA transcripts WITHOUT an annotation row    :       0
##  kept transcripts NOT in FASTA (must be 0)      :       0
##  kept genes are a subset of FASTA genes         : True
##
##RESULT: OK

# transcripts & genes in the FASTA
grep -c '^>' $FASTA
grep '^>' $FASTA | sed 's/^>//;s/[[:space:]].*//;s/_i[0-9]\+$//' | sort -u | wc -l

# transcripts & genes kept in the filtered annotation (skip header)
tail -n +2 $OUT | cut -f2 | sort -u | wc -l
tail -n +2 $OUT | cut -f1 | sort -u | wc -l

## Rename transcriptome:
cp $FASTA ../../gitignore/Rmegarrhizum_ChyKol2008_denovo_transcriptome.fasta
