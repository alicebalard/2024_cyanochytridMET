#!/bin/bash
# S03b_readAssignmentStats.sh
# Reviewer comment C4: per-library read numbers needed to
# report the fraction of reads assigned to the chytrid and cyanobacterium.
# Bowtie 1 + RSEM pipeline (S01-S04). Run after S02 (RSEM results present).
#
# For every library in samples_file.txt it reports:
#   input_pairs   read pairs given to Bowtie (records in the R1 fastq, i.e.
#                 after trimming and SortMeRNA rRNA removal) -> denominator
#   aligned_pairs fragments aligned by Bowtie and used by RSEM
#                 (N1, line 1 of RSEM.stat/RSEM.cnt)
#   unique/multi  uniquely / multi-mapping fragments (line 2 of RSEM.cnt)
#
# Bowtie was run with --no-unal (S02): unaligned reads never reach RSEM, so
# N0 in RSEM.cnt is 0 and the denominator must be counted from the fastq.
#
# Output: data/read_assignment_stats.tsv, used in S05_fullAnalysis.R
# (section "Read assignment per library"), where the aligned read pairs are
# split into chytrid / cyanobacterium with the decontaminated count matrix.
#SBATCH --job-name=S03b
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=06:00:00
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -uo pipefail

## Same folders as S02 / S03 (Trinity writes one folder per library in DATADIR)
DATADIR=/scratch/alicebalard/outData/RSEM
SAMPLE_FILE=$DATADIR/samples_file.txt
OUT=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/read_assignment_stats.tsv
TMP=$DATADIR/read_assignment_tmp
mkdir -p $TMP

## ---------------------------------------------------------------
## One library: input read pairs (fastq) + RSEM alignment counts
## ---------------------------------------------------------------
stats_one() {   # $1 = library name, $2 = R1 fastq(.gz)
    local rep=$1 fq1=$2 n cnt aligned uniq multi
    case "$fq1" in
        *.gz) n=$(zcat "$fq1" | wc -l) ;;
        *)    n=$(wc -l < "$fq1") ;;
    esac
    cnt=$DATADIR/$rep/RSEM.stat/RSEM.cnt
    if [ -s "$cnt" ]; then
        aligned=$(awk 'NR==1 {print $2}' "$cnt")   # line 1: N0 N1 N2 N_tot
        uniq=$(awk 'NR==2 {print $1}' "$cnt")      # line 2: nUnique nMulti nUncertain
        multi=$(awk 'NR==2 {print $2}' "$cnt")
    else
        echo "[S03b] WARNING: no RSEM.cnt for $rep ($cnt)" >&2
        aligned=NA; uniq=NA; multi=NA
    fi
    echo -e "${rep}\t$(( n / 4 ))\t${aligned}\t${uniq}\t${multi}" > $TMP/${rep}.tsv
    echo "[S03b] $rep: input=$(( n / 4 )) aligned=$aligned"
}
export -f stats_one
export DATADIR TMP

## ---------------------------------------------------------------
## All libraries, 8 at a time (columns of samples_file: condition, replicate, fq1, fq2)
## ---------------------------------------------------------------
cut -f2,3 $SAMPLE_FILE | grep -v '^\s*$' |
    xargs -P ${SLURM_CPUS_PER_TASK:-8} -n 2 bash -c 'stats_one "$0" "$1"'

## ---------------------------------------------------------------
## Combine, in the order of samples_file
## ---------------------------------------------------------------
echo -e "sample\tinput_pairs\taligned_pairs\tunique\tmulti" > $OUT
for rep in $(cut -f2 $SAMPLE_FILE); do
    if [ -s $TMP/${rep}.tsv ]; then cat $TMP/${rep}.tsv >> $OUT
    else echo "[S03b] ERROR: no result for $rep" >&2; fi
done
rm -rf $TMP

N_LIBS=$(( $(wc -l < $OUT) - 1 ))
echo "[S03b] Done: $OUT ($N_LIBS libraries, expected $(grep -c . $SAMPLE_FILE))"
## Totals per sample type in read pairs (x2 = individual reads, as in the old text)
awk -F'\t' 'NR > 1 {t = ($1 ~ /_chy_/) ? "zoospore" : ($1 ~ /_both_/) ? "infected filament" : "uninfected cyanobacterium"
                    s[t] += $2}
            END {for (t in s) printf "[S03b] %s samples: %.1f M read pairs (%.1f M reads)\n", t, s[t]/1e6, 2*s[t]/1e6}' $OUT
