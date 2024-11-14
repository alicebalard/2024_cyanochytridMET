module purge

## Prepare Trinotate
TRINOTATE_HOME=/scratch/alicebalard/outData/annotation/Trinotate

## Install Trinotate and associated programs
cd $TRINOTATE_HOME

## Add to my path
PATH=$PATH:$TRINOTATE_HOME

## https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required
## Trinity as module
module add Trinity/2.10.0-foss-2019b-Python-3.7.4
## Install Trinotate v4.0.2 in $TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2
## Add to my path
PATH=$PATH:$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2

## Install Transdecoder v5.7.1 $TRINOTATE_HOME/Transdecoder-Transdecoder-v5.7.1
## Add to my path
PATH=$PATH:$TRINOTATE_HOME/Transdecoder-Transdecoder-v5.7.1

## SQLite as module v 3.43.1
module load SQLite/3.43.1-GCCcore-13.2.0
## diamond as module v 2.0.13
module load BLAST+/2.13.0-gompi-2022a
module load DIAMOND/2.0.13-GCC-11.2.0

## Install Hmmer v 3.4 in $TRINOTATE_HOME/hmmer
## Add to my path
PATH=$PATH:$TRINOTATE_HOME/hmmer

## Install Infernal v1.1.5 in $TRINOTATE_HOME/infernal
## Add to my path
PATH=$PATH:$TRINOTATE_HOME/infernal/bin

## Install signalP v6.0 fast in a conda environment (myannot)
## Install tmhmm v2.0 DONE in same environment. NB some edits done following trinotate tuto page

export TRINOTATE_DATA_DIR=/scratch/alicebalard/outData/annotation/Trinotate/DATADIR
