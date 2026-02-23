# 2024_cyanochytridMET

Script for the article "Herbicide metolachlor alters gene expression and disrupts the interaction between a bloom-forming cyanobacterium and its chytrid parasite"

Authors: Alice Balard, Jürgen F. H. Strassert, Justyna Wolinska, and Erika B. Martínez-Ruiz

This is the pipeline and the files it should create (note: some of the files are not present any longer in my folder, must have been erased in the clean)

# Pipeline description:

# Part 1. Transcriptome

## Step 1.1: remove rRNA
**script:** scripts/Part1_transcriptome/1.1_runSortmeRNA.sh

**input:**
- downloaded database of rRNA: rRNA_databases_v4 
- trimmed data: /scratch/alicebalard/RawData/Sample_trimm_paired.fastq.gz

**output:** /scratch/alicebalard/outData/sortmerna/Sample_non_rRNA_fwd(or rev).fq.gz 

## Step 1.2: check how well the rRNA removal worked
**script:** scripts/Part1_transcriptome/1.2_checkEfficiency.sh

**input:** input & output of previous step

**output:** a number of reads printed on screen before and after rRNA removal

## Step 2: quality check after rRNA removal
**script:** scripts/Part1_transcriptome/2_fastqc_afterSortmeRNA.sh

**input:** /scratch/alicebalard/outData/sortmerna/*non_rRNA*fwd.fq.gz

**output:** fastqc files in /scratch/alicebalard/outData/fastqc/afterSortmerna

## Step 3: merge quality check after rRNA removal for all samples
**script:** scripts/Part1_transcriptome/3_multiQC_afterSortmeRNA.sh

**input:** /scratch/alicebalard/outData/fastqc/afterSortmerna

**output:** /scratch/alicebalard/outData/multiQC/sortmeRNA

## Step 4.1: assemble a de novo transcriptome for the sole chytrid R. megarrhizum
**script:** scripts/Part1_transcriptome/4.1_transcriptomeAssemblyRmega.sh

**input:** /scratch/alicebalard/outData/sortmerna/Z'$i'_non_rRNA_fwd.fq.gz & Z'$i'_non_rRNA_rev.fq.gz (combined)

**output:** /scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta

## Step 4.2: assemble a de novo transcriptome for the infecting chytrid R. megarrhizum *** to change Jurgen adding the extra RNAseq from co-culture experiment***
**script:** scripts/Part1_transcriptome/4.2_transcriptomeAssemblyInfectedFilaments.sh

**input:** /scratch/alicebalard/outData/sortmerna/In'$i'_non_rRNA_fwd.fq.gz & In'$i'_non_rRNA_rev.fq.gz

**output:** /scratch/alicebalard/outData/assembly_In/trinity_out_dir/Trinity.fasta

## Step 4.3: check basic contig statistics with Trinity tools
**script:** scripts/Part1_transcriptome/4.3_checkAssemblyQuality_basicStats.sh

**input:**
- assembly 1: /scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta 
- assembly 2: /scratch/alicebalard/outData/assembly_In/trinity_out_dir/Trinity.fasta

**output:** print basic transcriptome stats

## Step 5.1: diamond blastx for both assemblies against uniprot *** to change Jurgen for nr ?***
**script:** scripts/Part1_transcriptome/5.1_diamondBlastX.sh

**input:**
- assembly 1: /scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta 
- assembly 2: /scratch/alicebalard/outData/assembly_In/trinity_out_dir/Trinity.fasta

**output:** 
- /scratch/alicebalard/outData/diamondBlastX/assemblyZ_diamond_1e-5pval.out
- /scratch/alicebalard/outData/diamondBlastX/assemblyIn_diamond_1e-5pval.out


## Step 5.2: Select transcripts from both transcriptomes that are blastx-ing with fungi
**script:** scripts/Part1_transcriptome/5.2_selectFungi.R

**input:**
- /scratch/alicebalard/outData/diamondBlastX/assemblyZ_diamond_1e-5pval.out
- /scratch/alicebalard/outData/diamondBlastX/assemblyIn_diamond_1e-5pval.out
  
**output:** 
- /scratch/alicebalard/outData/diamondBlastX/assZ_Fungi_transcripts
- /scratch/alicebalard/outData/diamondBlastX/assIn_Fungi_transcripts

## Step 5.3: align the original reads to their transcriptomes
**script:** scripts/Part1_transcriptome/5.3_alignments.sh

**input:**
- both assemblies (/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta or /scratch/alicebalard/outData/assembly_In/trinity_out_dir/Trinity.fasta)
- both sets of combined reads from step 4.1 (/scratch/alicebalard/outData/assembly/combined_left.fq and right, and /scratch/alicebalard/outData/assembly_In/combined_left.fq and right)
  
**output:** 
- /scratch/alicebalard/outData/alignments/assemblyZ.reads.bam
- /scratch/alicebalard/outData/alignments/assemblyIn.reads.bam

## Step 5.4: filter out fungi reads with bedtools bamtobed and bash script
**script:** scripts/Part1_transcriptome/5.4_filterReads.sh

**input:**
- fungi transcripts from 5.2 (/scratch/alicebalard/outData/diamondBlastX/assZ_Fungi_transcripts and In)
- bam files from 5.3 (/scratch/alicebalard/outData/alignments/assemblyZ.reads.bam and In)

**output:**
Z_OUT=/scratch/alicebalard/outData/alignments/assemblyZ.reads_filteredFungi
In_OUT=/scratch/alicebalard/outData/alignments/assemblyIn.reads_filteredFungi

## Step 5.5: subset the fasta reads to keep only those mapping to fungal trasncripts
**script:** scripts/Part1_transcriptome/5.5_makeSubsetFastaFungalReads.sh

**input:** 
FUNGIREADS_NAMES from 5.4 (/scratch/alicebalard/outData/alignments/assemblyZ.reads_filteredFungi and In)
combined reads from 4.1 (/scratch/alicebalard/outData/assembly/combined_left.fq, right, and In left & right)

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/assemblyZ.reads_filteredFungi_left.fq & right & same for In

## Step 6.1: assemble the final transcriptome with Trinity
**script:** scripts/Part1_transcriptome/6.1_transcriptomeAssemblyFinal.sh

**input:** filtered reads from 5.5 lft & right for Z and In

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta

## Step 6.2: check basic contig statistics with Trinity tools
**script:** scripts/Part1_transcriptome/6.2_checkAssemblyQuality_basicStats.sh

**input:** assembly from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta)

**output:** print basic transcriptome stats

## Step 6.3: BUSCO of the transcriptome
**script:** scripts/Part1_transcriptome/6.3_BUSCOFungiAssembly.sh

**input:** assembly from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta)

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/BUSCO_MergedFungiTrannscriptome

## Step 7.1. prepare trinotate for annotation
**script:** scripts/Part1_transcriptome/7.1_finalTransAnnotation_prepareTrinotate.sh

all installation done in /scratch/alicebalard/outData/annotation/Trinotate to prepare the following

## Step 7.2. TransDecoder done in its own script
**script:** scripts/Part1_transcriptome/7.2_finalTransAnnotation_Transdecoder.sh

**input:** 
- assembly from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta)
- gene trans map from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta.gene_trans_map)

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/annotation/transdecoder

## Step 7.3. SignalP done in its own script
**script:** scripts/Part1_transcriptome/7.3_finalTransAnnotation_SignalP.sh

**input:** output of 7.2 (/scratch/alicebalard/outData/assemblyMergedFungi/annotation/transdecoder/Trinity.fasta.transdecoder.pep)

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/annotation/sigP6outdir/

## Step 7.4. all the rest of the annotation done by Trinotate in one script
**script:** scripts/Part1_transcriptome/7.4_finalTransAnnotation_Trinotate.sh

**input:** 
- assembly from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta)
- coding_seqs.pep from 7.2 (/scratch/alicebalard/outData/assemblyMergedFungi/annotation/transdecoder/Trinity.fasta.transdecoder.pep)
- gene_to_trans_map from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta.gene_trans_map)

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_simplified_GOKegg.tsv

## Step 7.5. Select only the transcripts with eukaryote hits
**script:** scripts/Part1_transcriptome/7.5_check_annot.R

**input:**
- output for 7.4 (/scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_simplified_GOKegg.tsv)
- assembly from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta)
- gene_to_trans_map from 6.1 (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta.gene_trans_map)

**output:**
- transcriptome filtered for eukaryote hits (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta)
- subset of gene trans map for eukaryote (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.gene_trans_map)
- subset annotation file for the Eukaryotic hits (/scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_filterEuk_simplified_GOKegg.tsv)

## Step 7.6: check basic contig statistics with Trinity tools
**script:** scripts/Part1_transcriptome/7.6_checkAssembliesFilter2stats.sh

**input:** output of 7.5 transcriptome filtered for eukaryote hits (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta)

**output:** print basic transcriptome stats

## Step 7.7: BUSCO of the transcriptome filtered for eukaryote hits
**script:** scripts/Part1_transcriptome/7.7_BUSCOifFilter.sh

**input:** output of 7.5 transcriptome filtered for eukaryote hits (/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta)

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/BUSCO_MergedFungiTranscriptome_eukaryoteHits

## Step 8. Decontaminate post DEG (to do after part 2)
**script:** 8_cleanTranscriptomeAfterDEG.sh

**input:** listOfTranscriptContaminant_toRmFromChytridTranscriptome from part 2 S04

**output:** /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.rmDEGconta.fasta

# Part 2. DEG

## Step 1. Prepare combined transcriptome
**script:** S01_prepareCombinedTranscriptome.sh

**input:**
T_CHY=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta
T_CHY_GTM=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.gene_trans_map
T_CYA=/scratch/alicebalard/outData/mergedTransc/GCF_904830935.1_P._agardhii_No.976_rna_from_genomic.fna

**output:**
/scratch/erikamr/cyano_chytrid_met/data/combined_gene_trans_map_cds_final_hope.txt
/scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds_final_hope.fna

## Step 2. Mapping with bowtie 
**script:** S02_mappingBowtie.sh

**input:**
prepared in S01: 
ASSEMBLY_BOTH=/scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds_final_hope.fna
GTM=/scratch/erikamr/cyano_chytrid_met/data/combined_gene_trans_map_cds_final_hope.txt

manually prepared (removing bad quality samples):                                             
SAMPLE_FILE=/scratch/erikamr/cyano_chytrid_met/data/samples_file_remove.txt

**output:** found in OUTDIR=/scratch/erikamr/cyano_chytrid_met/data/out_trinity_align_rsem_final_hope

## Step 3. Generate count matrix
**script:** S03_make_count_matrix_RSEM.sh

**input:** 
Created in S01: GTM=/scratch/erikamr/cyano_chytrid_met/data/combined_gene_trans_map_cds_final_hope.txt
Manually made listing S02 results: QUANT=/scratch/erikamr/cyano_chytrid_met/data/quant_file_isoform_final_hope.txt

**output:** in: /scratch/erikamr/cyano_chytrid_met/scripts/RSEM_final_hope_matrices

## Step 4. DEG analysis in R
**script:** S04_fullAnalysis (& dataload.R): generate all results for DEG, as well as data/listOfTranscriptContaminant_toRmFromChytridTranscriptome useful in part 1 script 8

# Part 3. WGCNA
**script:** networkAnalysis.R Builds up on part 2
