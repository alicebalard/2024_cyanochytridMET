#!/bin/bash

#SBATCH --job-name=dlUniprot
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --mem=100
#SBATCH --time=1-24:00:00                       
#SBATCH --qos=standard              

cd /scratch/alicebalard/outData/blobtools/

wget -c -q -O uniprot/reference_proteomes.tar.gz \
  ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
      -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
      awk '/tar.gz/ {print $9}')

echo "reference_proteomes.tar.gz has been DL"

cd uniprot
tar xf reference_proteomes.tar.gz

echo "reference_proteomes.tar.gz has been untared"

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

echo "reference_proteomes.fasta.gz has been created"

echo -e "accession\taccession.version\ttaxid\tgi" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

echo "reference_proteomes.taxid_map has been created. Now you need to run diamond makedb with diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd"

