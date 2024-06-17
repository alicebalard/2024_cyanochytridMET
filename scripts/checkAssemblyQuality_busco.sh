#!/bin/sh --login
#SBATCH --job-name=BUSCO
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8GB 
#SBATCH --time=24:00:00
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu
#SBATCH --qos=standard              
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end

cd /scratch/alicebalard/outData/qualityAssembly/BUSCO

module add BUSCO/5.1.2-foss-2020b

TRANSC=/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta

## With more samples
busco -i $TRANSC -l fungi_odb10 -o output_directory -m transcriptome -c 20 -o "BUSCO_transcZ1to12"

## l = current fungi database (will be downloaded automatically)
## m = genome or transcriptome (in your case transcriptome)

## NOTE: BUSCO doesn't like slash signs in the fasta header, you may need to replace them before running BUSCO.



############################################################
## Previous results (only Z2 sample, before decontamination)
# ***** Results: *****
#
#        C:8.7%[S:3.2%,D:5.5%],F:13.7%,M:77.6%,n:758        
#        66      Complete BUSCOs (C)                        
#        24      Complete and single-copy BUSCOs (S)        
#        42      Complete and duplicated BUSCOs (D)         
#        104     Fragmented BUSCOs (F)                      
#        588     Missing BUSCOs (M)                         
#        758     Total BUSCO groups searched


