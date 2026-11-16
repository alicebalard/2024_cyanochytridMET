#!/bin/bash
# S03b_readAssignmentStats.sh
# Reviewer comment C2 (and new Fig. 2): per-library read numbers needed to
# report the fraction of reads assigned to the chytrid and cyanobacterium.
# Bowtie 1 + RSEM pipeline (S01-S04). Run after S02 (RSEM results present).
#
# For every library in samples_file.txt it reports:
#   trimmed_pairs read pairs after trimming, BEFORE SortMeRNA (records in the  ## NEW
#                 R1 *_trimm_paired.fastq.gz in RawData, input of 1.1_runSortmeRNA.sh)
#   input_pairs   read pairs given to Bowtie (records in the R1 fastq, i.e.
#                 after trimming and SortMeRNA rRNA removal) -> denominator
#   aligned_pairs fragments aligned by Bowtie and used by RSEM
#                 (N1, line 1 of RSEM.stat/RSEM.cnt)
#   unique/multi  uniquely / multi-mapping fragments (line 2 of RSEM.cnt)
#   pct_after_rRNA  100 * input_pairs / trimmed_pairs (% of pairs kept by SortMeRNA)  ## NEW
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
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -uo pipefail

## Same folders as S02 / S03 (Trinity writes one folder per library in DATADIR)
DATADIR=/scratch/alicebalard/outData/RSEM
SAMPLE_FILE=$DATADIR/samples_file.txt
OUT=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/data/read_assignment_stats.tsv
TMP=$DATADIR/read_assignment_tmp
## NEW: trimmed reads (input of SortMeRNA, see Part1_transcriptome/1.1_runSortmeRNA.sh)
RAWDIR=/scratch/alicebalard/RawData
TRIM_EXT=_trimm_paired.fastq.gz
mkdir -p $TMP "$(dirname $OUT)"     # create the output folder if it does not exist

## The RSEM folders of S02 are optional: if they were deleted (e.g. /scratch
## clean-up), aligned/unique/multi are NA. Only input_pairs is used in S05.
if ! find $DATADIR -path "*/RSEM.stat/RSEM.cnt" -print -quit 2> /dev/null | grep -q .; then
    echo "[S03b] NOTE: no RSEM output found below $DATADIR -> only input_pairs is computed" >&2
    NO_RSEM=1
else
    NO_RSEM=0
fi
export NO_RSEM

## ---------------------------------------------------------------
## One library: input read pairs (fastq) + RSEM alignment counts
## ---------------------------------------------------------------
stats_one() {   # $1 = library name, $2 = R1 fastq(.gz)
    local rep=$1 fq1=$2 n cnt aligned uniq multi
    case "$fq1" in
        *.gz) n=$(zcat "$fq1" | wc -l) ;;
        *)    n=$(wc -l < "$fq1") ;;
    esac
    ## NEW: read pairs after trimming, before rRNA removal
    ## library name -> raw sample name (control_chy_Z1 -> Z1, met_both_In7 -> In7)
    local short=${rep##*_} trim_fq ntrim trimmed pct
    trim_fq=$RAWDIR/${short}_1${TRIM_EXT}
    if [ -s "$trim_fq" ]; then
        ntrim=$(zcat "$trim_fq" | wc -l)
        trimmed=$(( ntrim / 4 ))
        pct=$(awk -v a=$(( n / 4 )) -v b=$trimmed 'BEGIN {printf "%.1f", 100 * a / b}')
    else
        echo "[S03b] WARNING: trimmed fastq not found for $rep: $trim_fq" >&2
        trimmed=NA; pct=NA
    fi
    ## end NEW
    ## RSEM folders are not always directly in DATADIR (S02 output layout):
    ## search for <library>/RSEM.stat/RSEM.cnt anywhere below DATADIR
    local hits nhits
    hits=$(find $DATADIR -path "*/${rep}/RSEM.stat/RSEM.cnt" 2> /dev/null)
    nhits=$(echo "$hits" | grep -c .)
    if [ "$nhits" -gt 1 ]; then
        echo "[S03b] ERROR: several RSEM.cnt for $rep:" >&2; echo "$hits" >&2
        cnt=""
    else
        cnt=$hits
    fi
    if [ -n "$cnt" ] && [ -s "$cnt" ]; then
        echo "[S03b] $rep: RSEM statistics from $cnt"
        aligned=$(awk 'NR==1 {print $2}' "$cnt")   # line 1: N0 N1 N2 N_tot
        uniq=$(awk 'NR==2 {print $1}' "$cnt")      # line 2: nUnique nMulti nUncertain
        multi=$(awk 'NR==2 {print $2}' "$cnt")
    else
        [ "$NO_RSEM" -eq 1 ] || echo "[S03b] WARNING: no RSEM.cnt found for $rep below $DATADIR" >&2
        aligned=NA; uniq=NA; multi=NA
    fi
    echo -e "${rep}\t${trimmed}\t$(( n / 4 ))\t${pct}\t${aligned}\t${uniq}\t${multi}" > $TMP/${rep}.tsv   ## NEW: 2 columns added
    echo "[S03b] $rep: trimmed=$trimmed after_rRNA=$(( n / 4 )) ($pct%) aligned=$aligned"
}
export -f stats_one
export DATADIR TMP RAWDIR TRIM_EXT   ## NEW: RAWDIR TRIM_EXT

## ---------------------------------------------------------------
## All libraries, 8 at a time (columns of samples_file: condition, replicate, fq1, fq2)
## ---------------------------------------------------------------
cut -f2,3 $SAMPLE_FILE | grep -v '^\s*$' |
    xargs -P ${SLURM_CPUS_PER_TASK:-8} -n 2 bash -c 'stats_one "$0" "$1"'

## ---------------------------------------------------------------
## Combine, in the order of samples_file
## ---------------------------------------------------------------
echo -e "sample\ttrimmed_pairs\tinput_pairs\tpct_after_rRNA\taligned_pairs\tunique\tmulti" > $OUT   ## NEW: header
for rep in $(cut -f2 $SAMPLE_FILE); do
    if [ -s $TMP/${rep}.tsv ]; then cat $TMP/${rep}.tsv >> $OUT
    else echo "[S03b] ERROR: no result for $rep" >&2; fi
done
rm -rf $TMP

N_LIBS=$(( $(wc -l < $OUT) - 1 ))
echo "[S03b] Done: $OUT ($N_LIBS libraries, expected $(grep -c . $SAMPLE_FILE))"
## NEW: totals per sample type, before (trimmed) and after rRNA removal, in read pairs
## (x2 = individual reads). % = pooled pairs after / pooled pairs before, as in the text.
## Columns: $2 trimmed_pairs, $3 input_pairs
awk -F'\t' 'NR > 1 {t = ($1 ~ /_chy_/) ? "zoospore" : ($1 ~ /_both_/) ? "infected filament" : "uninfected cyanobacterium"
                    a[t] += $3; if ($2 != "NA") {b[t] += $2; ab[t] += $3} else miss[t]++}
            END {for (t in a) {
                   printf "[S03b] %s samples: trimmed %.1f M pairs -> after rRNA removal %.1f M pairs (%.1f M reads)", t, b[t]/1e6, a[t]/1e6, 2*a[t]/1e6
                   if (b[t] > 0) printf ", %.1f%% kept", 100*ab[t]/b[t]
                   if (miss[t]) printf " [%d libraries without trimmed fastq]", miss[t]
                   printf "\n" }}' $OUT
## NEW: per-library range of % kept (median is clearer than the pooled % if libraries differ in size)
for t in _chy_ _both_ _cyano_; do
    awk -F'\t' -v t=$t 'NR > 1 && $1 ~ t && $4 != "NA" {print $4}' $OUT | sort -g |
        awk -v t=$t '{v[NR]=$1} END {if (NR) printf "[S03b] %s: %% kept per library median %.1f (range %.1f-%.1f, n=%d)\n", t, (NR%2 ? v[(NR+1)/2] : (v[NR/2]+v[NR/2+1])/2), v[1], v[NR], NR}'
done

## Results:

# Done: /home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/data/read_assignment_stats.tsv (36 libraries, expected 36)
#[S03b] infected filament samples: trimmed 500.6 M pairs -> after rRNA removal 440.6 M pairs (881.1 M reads), 88.0% kept
#[S03b] uninfected cyanobacterium samples: trimmed 482.0 M pairs -> after rRNA removal 455.0 M pairs (910.1 M reads), 94.4% kept
#[S03b] zoospore samples: trimmed 498.4 M pairs -> after rRNA removal 188.1 M pairs (376.1 M reads), 37.7% kept
#[S03b] _chy_: % kept per library median 34.4 (range 5.9-65.2, n=12)
#[S03b] _both_: % kept per library median 93.0 (range 46.1-99.1, n=12)
#[S03b] _cyano_: % kept per library median 95.0 (range 91.6-96.5, n=12)

