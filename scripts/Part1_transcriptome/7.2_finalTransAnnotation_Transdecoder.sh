#!/bin/bash

#SBATCH --job-name=7.2_transdecoder
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB 
#SBATCH --time=1-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

cd /scratch/alicebalard/outData/annotation/Trinotate/TransDecoder-TransDecoder-v5.7.1

## Final filtered transcriptome:
TRAN=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta
GTM=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta.gene_trans_map
OUT=/scratch/alicebalard/outData/assemblyMergedFungi/annotation/transdecoder

## Step 1: extract the long open reading frames
perl TransDecoder.LongOrfs -t $TRAN --gene_trans_map $GTM --output_dir $OUT

## Step 3: predict the likely coding regions
perl TransDecoder.Predict -t $TRAN --output_dir $OUT
