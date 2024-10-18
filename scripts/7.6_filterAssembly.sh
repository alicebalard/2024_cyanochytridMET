#!/bin/sh --login
#SBATCH --job-name=7.6_btk_filterAssembly
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1GB 
#SBATCH --time=4:00:00
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu
#SBATCH --qos=standard              
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end

module purge
module load Anaconda3
source ~/.bashrc
conda init --all
conda activate btk
module load GCCcore/12.3.0

ASSEMBLYDIR=/scratch/alicebalard/outData/assemblyFinal
ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.fasta 
BTKOUT=/scratch/alicebalard/outData/blobtools/FINALTran
cd $BTKOUT

## echo "Filter only transcripts belonging to the Fungi kingdom:"
##  
## ## keep only fungi
## blobtools filter --param bestsumorder_kingdom--Keys=Fungi --invert --fasta $ASSEMBLY $BTKOUT
##   
## ## output in ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.fasta. Rename output:
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.Fungi.fasta

ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.Fungi.fasta

## ########
## ## BUSCO
## echo "BUSCO:"
## module load BUSCO/5.1.2-foss-2020b
## busco -i $ASSEMBLY -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_Trinity.filtered.Fungi.fasta"
## 
## ###################
## ## Basic statistics
## echo "basic contig statistics"
## module load Trinity/2.10.0-foss-2019b-Python-3.7.4
## $TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY
## 
## ## Only chytrids and close fungal hits (e-value <1.10e-25)
## 
## ## 1. only Chytridiomycota
## blobtools filter --param bestsumorder_phylum--Keys=Chytridiomycota --invert --fasta $ASSEMBLY $BTKOUT
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.Fungi.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Chytridiomycota.filtered.fasta
## 
## ## 2. other fungi with low e-value
## 
## ### 2.1. other fungi
## blobtools filter --param bestsumorder_phylum--Keys=Chytridiomycota --fasta $ASSEMBLY $BTKOUT
## 
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.Fungi.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.FungiNotChytridiomycota.filtered.fasta
## 
## ### 2.2. with low p-value
## blobtools filter --param length--Min=5000 --fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.FungiNotChytridiomycota.filtered.fasta $BTKOUT
## 
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.FungiNotChytridiomycota.filtered.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.FungiNotChytridiomycota_5000+.filtered.fasta
## 
## ### merge
## cat $ASSEMBLYDIR/trinity_out_dir/Trinity.Chytridiomycota.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.FungiNotChytridiomycota_5000+.filtered.fasta > $ASSEMBLYDIR/trinity_out_dir/Trinity.ChytridAndFungi5k+.filtered.fasta
 

ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.ChytridAndFungi5k+.filtered.fasta

echo "Chytrid and fungi 5k+:"
########
## BUSCO
echo "BUSCO:"
module load BUSCO/5.1.2-foss-2020b
busco -i $ASSEMBLY -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_Trinity.filtered.ChytridAndFungi5k+.fasta"

###################
## Basic statistics
echo "basic contig statistics"
module load Trinity/2.10.0-foss-2019b-Python-3.7.4
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY











##  ##############################
##  ## Ascomycota - Basidiomycota = higher fungi (subkingdom Dikarya)
##  ## Mucoromycota - close to Dikarya
##  ## Chytridiomycota - Blastocladiomycota - Olpidiomycota - Cryptomycota - Microsporidia - Zoopagomycota = early diverged fungi
##  ###############################
##  ## Assembly keeping only the early diverged fungi (rm Ascomycota - Basidiomycota- Mucoromycota)
##  
## echo "Filter only transcripts belonging to the early diversed Fungi (remove Ascomycota - Basidiomycota- Mucoromycota):"
## 
## ## filter 1
## blobtools filter --param bestsumorder_phylum--Keys=Ascomycota --fasta $ASSEMBLY $BTKOUT
## 
## ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.Fungi.filtered.fasta
## 
## ## filter 2
## blobtools filter --param bestsumorder_phylum--Keys=Basidiomycota --fasta $ASSEMBLY $BTKOUT
## 
## ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.Fungi.filtered.filtered.fasta
## 
## ## echo "Without Asco and Basidio:"
## ## ########
## ## ## BUSCO
## ## echo "BUSCO:"
## ## module load BUSCO/5.1.2-foss-2020b
## ## busco -i $ASSEMBLY -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_Trinity.filtered.woAscoBasidio.fasta"
## ## 
## ## ###################
## ## ## Basic statistics
## ## echo "basic contig statistics"
## ## module load Trinity/2.10.0-foss-2019b-Python-3.7.4
## ## $TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY
## ## 
## ############# Rm Mucoro
## echo "Without Mucoro as well:"
## 
## ## filter 3
## blobtools filter --param bestsumorder_phylum--Keys=Mucoromycota --fasta $ASSEMBLY $BTKOUT
## 
## ## rename output:
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.Fungi.filtered.filtered.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.EarlydivFungi.fasta 
## 
## ########
## ## BUSCO
## echo "BUSCO:"
## module load BUSCO/5.1.2-foss-2020b
## ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.EarlydivFungi.fasta
## 
## busco -i $ASSEMBLY -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_Trinity.filtered.EarlydivFungi"
## 
## ###################
## ## Basic statistics
## echo "basic contig statistics"
## module load Trinity/2.10.0-foss-2019b-Python-3.7.4
## 
## $TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY
## 
########################
## Only chytrids & close (Blastocladiomycota & Olpidiomycota & Fungi-undef)

## blobtools filter --param bestsumorder_phylum--Keys=Chytridiomycota --invert --fasta $ASSEMBLY $BTKOUT
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Chytridiomycota.filtered.fasta
## 
## blobtools filter --param bestsumorder_phylum--Keys=Blastocladiomycota --invert --fasta $ASSEMBLY $BTKOUT
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Blastocladiomycota.filtered.fasta
## 
## blobtools filter --param bestsumorder_phylum--Keys=Olpidiomycota --invert --fasta $ASSEMBLY --output Olpidiomycota $BTKOUT
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Olpidiomycota.filtered.fasta
## 
## blobtools filter --param bestsumorder_phylum--Keys=Fungi-undef --invert --fasta $ASSEMBLY --output Fungi-undef $BTKOUT
## mv $ASSEMBLYDIR/trinity_out_dir/Trinity.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Fungi-undef.filtered.fasta
## 
## ## Merge
## cat $ASSEMBLYDIR/trinity_out_dir/Trinity.Chytridiomycota.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Blastocladiomycota.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Olpidiomycota.filtered.fasta $ASSEMBLYDIR/trinity_out_dir/Trinity.Fungi-undef.filtered.fasta > $ASSEMBLYDIR/trinity_out_dir/Trinity.allchytrids.filtered.fasta
## 
## ## BUSCO
## module add BUSCO/5.1.2-foss-2020b
## 
## busco -i $ASSEMBLYDIR/trinity_out_dir/Trinity.allchytrids.filtered.fasta -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_Trinity.filtered.allchytrids"
## 
## ###################
## ## Basic statistics
## echo "basic contig statistics"
## module load Trinity/2.10.0-foss-2019b-Python-3.7.4
## 
## $TRINITY_HOME/util/TrinityStats.pl $ASSEMBLYDIR/trinity_out_dir/Trinity.allchytrids.filtered.fasta
## 
## 
