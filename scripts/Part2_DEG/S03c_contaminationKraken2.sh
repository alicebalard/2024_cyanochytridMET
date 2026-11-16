#!/bin/bash
# S03c_contaminationKraken2.sh
# Taxonomic composition of the reads remaining after rRNA removal
# (input of Bowtie, same files as S03b), to estimate how many reads come
# from organisms other than P. agardhii and R. megarrhizum.
#
# Per library: random subsample of read pairs -> Kraken2 (paired) -> report.
# Then one summary table: % Planktothrix / other cyanobacteria / other bacteria /
# fungi / other eukaryotes / unclassified.
# NB: R. megarrhizum is not in standard Kraken2 databases -> most chytrid reads
# will be "unclassified" (or assigned to Fungi / Chytridiomycota at best).
#
# Run:  sbatch S03c_contaminationKraken2.sh          (array: one task per library)
#       then, after all tasks finished:  bash S03c_contaminationKraken2.sh summary
#SBATCH --job-name=S03c
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --array=1-36
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/scripts/logs_dir/%x.%A_%a.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/scripts/logs_dir/%x.%A_%a.err
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -uo pipefail

SAMPLE_FILE=/scratch/alicebalard/outData/RSEM/samples_file.txt
OUTDIR=/scratch/alicebalard/outData/kraken2_contamination
OUT=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/data/kraken2_composition.tsv

## Kraken2 (v.2.17.2) and seqtk (v.1.5-r133) installed in conda env kraken2_env

## path to a Kraken2 database: PlusPF = Refseq archaea, bacteria, viral, plasmid, human1, UniVec_Core, plus Refseq protozoa & fungi -8GB to fit in less memory
KDB=/scratch/alicebalard/db/k2_pluspf_8gb
NSUB=200000          # read pairs per library (enough for % estimates to +-0.2%)
mkdir -p $OUTDIR "$(dirname $OUT)"

## =============================================================
## Summary mode: combine all Kraken2 reports into one table
## =============================================================
if [ "${1:-}" = "summary" ]; then
    echo -e "sample\tpct_Planktothrix\tpct_other_Cyanobacteria\tpct_other_Bacteria\tpct_Fungi\tpct_other_Eukaryota\tpct_Viruses_Archaea\tpct_unclassified" > $OUT
    for rep in $(cut -f2 $SAMPLE_FILE); do
        r=$OUTDIR/${rep}.report
        [ -s $r ] || { echo "[S03c] missing report for $rep" >&2; continue; }
        ## Kraken2 report columns: 1 % clade, 4 rank code, 6 name (indented)
        awk -F'\t' -v rep=$rep '
            {name=$6; gsub(/^ +/, "", name)}
            $4=="U"                          {u=$1}
            $4=="D" && name=="Bacteria"      {bac=$1}
            $4=="D" && name=="Eukaryota"     {euk=$1}
            $4=="D" && (name=="Viruses" || name=="Archaea") {va+=$1}
            $4=="K" && name=="Fungi"         {fun=$1}
            $4=="P" && name=="Cyanobacteriota" {cya=$1}
            $4=="P" && name=="Cyanobacteria"   {cya=$1}
            $4=="G" && name=="Planktothrix"  {pla=$1}
            END {printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
                 rep, pla, cya-pla, bac-cya, fun, euk-fun, va, u}' $r >> $OUT
    done
    echo "[S03c] written $OUT"
    ## top 10 other genera per sample type (what are the contaminants?)
    for t in _chy_ _both_ _cyano_; do
        echo "== top genera (excluding Planktothrix), mean % over $t libraries"
        for r in $OUTDIR/*${t}*.report; do
            awk -F'\t' '$4=="G" {n=$6; gsub(/^ +/,"",n); if (n!="Planktothrix") print n"\t"$1}' $r
        done | awk -F'\t' -v k=$(ls $OUTDIR/*${t}*.report | wc -l) '{s[$1]+=$2} END {for (g in s) printf "%.2f\t%s\n", s[g]/k, g}' |
            sort -gr | head -10
    done
    exit 0
fi

## =============================================================
## Array mode: one library per task
## =============================================================
module purge
set +u                                   # conda/bashrc scripts use unset variables
source /trinity/shared/easybuild/software/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate kraken2_env
set -u                                   # back to strict mode for the rest of the script

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_FILE)
rep=$(echo "$line" | cut -f2); fq1=$(echo "$line" | cut -f3); fq2=$(echo "$line" | cut -f4)
echo "[S03c] $rep: $fq1 $fq2"

## same seed for R1 and R2 -> pairs stay together
seqtk sample -s 11 $fq1 $NSUB > $OUTDIR/${rep}_sub_1.fq
seqtk sample -s 11 $fq2 $NSUB > $OUTDIR/${rep}_sub_2.fq

kraken2 --db $KDB --threads ${SLURM_CPUS_PER_TASK:-8} --paired \
        --confidence 0.1 \
        --report $OUTDIR/${rep}.report --output /dev/null \
        $OUTDIR/${rep}_sub_1.fq $OUTDIR/${rep}_sub_2.fq

rm $OUTDIR/${rep}_sub_?.fq
echo "[S03c] $rep done"

##[S03c] written /home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/data/kraken2_composition.tsv                                                                                      
##== top genera (excluding Planktothrix), mean % over _chy_ libraries
##0.94    Hydrogenophaga
##0.17    Blastomonas
##0.11    Sphingorhabdus
##0.04    Lacibacter
##0.04    Homo
##0.02    Fluviicola
##0.02    Ferrovibrio
##0.02    Acinetobacter
##0.01    Pseudomonas
##0.01    Niveispirillum
##== top genera (excluding Planktothrix), mean % over _both_ libraries
##0.57    Hydrogenophaga
##0.18    Homo
##0.17    Acinetobacter
##0.11    Fluviicola
##0.11    Blastomonas
##0.05    Sphingorhabdus
##0.04    Cutibacterium
##0.03    Lacibacter
##0.03    Burkholderia
##0.03    Acidovorax
##== top genera (excluding Planktothrix), mean % over _cyano_ libraries
##0.46    Hydrogenophaga
##0.03    Pseudomonas
##0.02    Blastomonas
##0.01    Methylibium
##0.01    Lacibacter
##0.01    Homo
##0.01    Ferrovibrio
##0.01    Acinetobacter
##0.00    Xanthomonas
##0.00    Variovorax

