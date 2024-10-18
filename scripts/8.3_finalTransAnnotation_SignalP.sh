#!/bin/bash

#SBATCH --job-name=8.3_signalP
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2GB 
#SBATCH --time=10:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module load Anaconda3

source ~/.bashrc
conda init --all

cd /scratch/alicebalard/outData/annotation/Trinotate

coding_seqs=/scratch/alicebalard/outData/annotation/Trinotate/TransDecoder-TransDecoder-v5.7.1/Trinity.filtered.Fungi.fasta.transdecoder.pep

conda activate myannot
signalp6 --fastafile $coding_seqs --output_dir sigP6outdir --format none --organism euk --mode fast --torch_num_threads 20
conda deactivate
