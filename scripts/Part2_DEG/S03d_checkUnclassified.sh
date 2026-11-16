#!/bin/bash
# S03d_checkUnclassified.sh
# What are the reads that Kraken2 cannot classify? (one library, C1 by default)
#  (a) Kraken2 WITHOUT --confidence on 50,000 read pairs -> % unclassified
#  (b) unclassified reads: length, poly-G / poly-A, duplication, top sequences
# Usage: sbatch S03d_checkUnclassified.sh        (or: sbatch S03d_checkUnclassified.sh In3)
# Results: in the .out log (logs_dir) + files in $T
#SBATCH --job-name=S03d
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/2024_cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -uo pipefail

LIB=${1:-C1}                                     # raw sample name (C1, Z3, In5...)
D=/scratch/alicebalard/outData/sortmerna
KDB=/scratch/alicebalard/db/k2_pluspf_8gb
T=/scratch/alicebalard/outData/kraken2_contamination/test_${LIB}
mkdir -p $T

set +u                                           # conda scripts use unset variables
source /trinity/shared/easybuild/software/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate kraken2_env
set -u

## 50,000 read pairs (same seed for R1 and R2 -> pairs stay together)
seqtk sample -s 11 $D/${LIB}_non_rRNA_fwd.fq.gz 50000 > $T/r1.fq
seqtk sample -s 11 $D/${LIB}_non_rRNA_rev.fq.gz 50000 > $T/r2.fq

echo "=========== (a) Kraken2 without confidence threshold ($LIB) ==========="
kraken2 --db $KDB --threads ${SLURM_CPUS_PER_TASK:-4} --paired \
        --report $T/noconf.report --output /dev/null \
        --unclassified-out $T/uncl#.fq \
        $T/r1.fq $T/r2.fq
echo "% unclassified (no threshold):"
grep -P "\tU\t" $T/noconf.report
echo "Top classified taxa (genus level):"
awk -F'\t' '$4=="G"' $T/noconf.report | sort -t$'\t' -k1,1gr | head -10

U=$T/uncl_1.fq
NU=$(( $(wc -l < $U) / 4 ))
echo
echo "=========== (b) unclassified reads (R1): $NU of 50000 ==========="
echo "-- most frequent read lengths (count length):"
awk 'NR%4==2 {print length($0)}' $U | sort -n | uniq -c | sort -nr | head -5
echo "-- length summary (min / median / max):"
awk 'NR%4==2 {print length($0)}' $U | sort -n | awk '{a[NR]=$1} END {print a[1], a[int((NR+1)/2)], a[NR]}'
echo "-- reads with poly-G (20 G):  $(awk 'NR%4==2' $U | grep -c 'GGGGGGGGGGGGGGGGGGGG')"
echo "-- reads with poly-A (20 A):  $(awk 'NR%4==2' $U | grep -c 'AAAAAAAAAAAAAAAAAAAA')"
echo "-- reads with poly-T (20 T):  $(awk 'NR%4==2' $U | grep -c 'TTTTTTTTTTTTTTTTTTTT')"
echo "-- distinct sequences:         $(awk 'NR%4==2' $U | sort -u | wc -l)"
echo "-- reads in sequences seen >1x: $(awk 'NR%4==2' $U | sort | uniq -c | awk '$1>1 {s+=$1} END {print s+0}')"
echo "-- 10 most repeated sequences (count sequence) -> BLAST the top ones on NCBI:"
awk 'NR%4==2' $U | sort | uniq -c | sort -nr | head -10
echo "-- GC content of unclassified reads (mean %):"
awk 'NR%4==2 {n+=length($0); g+=gsub(/[GC]/,"")} END {printf "%.1f\n", 100*g/n}' $U
echo "-- 10 random unclassified reads in FASTA (paste into NCBI BLASTn):"
awk 'NR%4==1 {h=$1} NR%4==2 {print ">"substr(h,2)"\n"$0}' $U | \
    awk 'NR%2==1 {h=$0; next} {print h"\t"$0}' | shuf -n 10 --random-source=<(yes) | tr '\t' '\n'

echo
echo "[S03d] done. Files in $T"
