#!/bin/bash

#SBATCH --job-name=GTaxdecont
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                             
#SBATCH --mem=60GB 
#SBATCH --time=24:00:00
#SBATCH --qos=standard              
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module purge 

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module load Python/3.9.6-GCCcore-11.2.0

cd /scratch/alicebalard/GTax/

## gtax installed from https://gtax.readthedocs.io/en/latest/installation.html
## Activate python environment
source gtax_venv/bin/activate

## datasets downloaded from https://gtax.readthedocs.io/en/latest/datasets.html
./datasets download genome taxon 2 --assembly-source refseq --dehydrated --filename bacteria_meta.zip

## read the zipped metadata file for each superkingdom and create the folders for hydration with the datasets command.
## This command will keep the reference genome for each taxa if it is available. If no reference genome is available,
## the latest assembly will be kept.
filter_metadata_zip

## Hydrate directories with datasets
### Archaea
./datasets rehydrate --directory archaea/

### Bacteria
./datasets rehydrate --directory bacteria/

### Viruses
./datasets rehydrate --directory viruses/

### Eukaryotes
./datasets rehydrate --directory eukaryotes/

## Create Gtax FASTA files
## After all data is downloaded, it will take few hours to finish, we can create the FASTA, indexes and TaxID maps for the databases.
gtax_database
