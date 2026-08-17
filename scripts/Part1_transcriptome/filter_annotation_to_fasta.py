#!/usr/bin/env python3
"""
Filter a Trinotate simplified annotation table to keep only the transcripts
present in a (decontaminated) Trinity FASTA, and report gene/transcript counts.

Matching is done on the TRANSCRIPT id (isoform level, e.g. TRINITY_DN839_c0_g1_i103),
because both the FASTA headers and the annotation `transcript_id` column are at
that level. Genes (TRINITY_DN839_c0_g1) are counted separately as a cross-check.

Usage:
  python3 filter_annotation_to_fasta.py \
      --fasta Trinity_eukaryoteHits.rmConta.fasta \
      --annot assemblyMergedFungi_filterEuk_simplified_GOKegg.tsv \
      --out   Rmegarrhizum_ChyKol2008_denovo_transcriptome_annotation.tsv
"""
import argparse, re, sys

def gene_of(tx):
    # TRINITY_DN839_c0_g1_i103 -> TRINITY_DN839_c0_g1
    return re.sub(r'_i\d+$', '', tx)

def fasta_transcript_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                # take the token right after '>', drop any trailing "len=.. path=.." etc.
                ids.append(line[1:].split()[0])
    return ids

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--annot', required=True)
    ap.add_argument('--out',   required=True)
    ap.add_argument('--tx-col', default='transcript_id',
                    help="name of the transcript-id column in the annotation header")
    args = ap.parse_args()

    # --- FASTA ---
    fa_tx_list = fasta_transcript_ids(args.fasta)
    fa_tx = set(fa_tx_list)
    fa_genes = {gene_of(t) for t in fa_tx}
    if len(fa_tx_list) != len(fa_tx):
        print(f"WARNING: FASTA has {len(fa_tx_list)-len(fa_tx)} duplicate header id(s)",
              file=sys.stderr)

    # --- annotation ---
    with open(args.annot) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            tx_i = header.index(args.tx_col)
            g_i  = header.index('X.gene_id')
        except ValueError:
            sys.exit(f"ERROR: expected columns '{args.tx_col}' and 'X.gene_id' in header: {header}")

        ann_tx_before, ann_g_before = set(), set()
        kept_tx, kept_g = set(), set()
        n_rows_before = n_rows_kept = 0

        with open(args.out, 'w') as out:
            out.write('\t'.join(header) + '\n')
            for line in fh:
                if not line.strip():
                    continue
                n_rows_before += 1
                f = line.rstrip('\n').split('\t')
                tx, g = f[tx_i], f[g_i]
                ann_tx_before.add(tx); ann_g_before.add(g)
                if tx in fa_tx:
                    out.write(line if line.endswith('\n') else line + '\n')
                    n_rows_kept += 1
                    kept_tx.add(tx); kept_g.add(g)

    # --- report ---
    orphans   = kept_tx - fa_tx                      # should always be empty (built from fa_tx)
    fa_no_ann = fa_tx - ann_tx_before                # FASTA transcripts with no annotation row
    removed   = ann_tx_before - fa_tx                # annotation transcripts dropped (contaminants etc.)

    print("== FASTA (decontaminated transcriptome) ==")
    print(f"  transcripts (headers) : {len(fa_tx):>7}")
    print(f"  genes (unique)        : {len(fa_genes):>7}")
    print("== Annotation BEFORE filtering ==")
    print(f"  rows                  : {n_rows_before:>7}")
    print(f"  transcripts (unique)  : {len(ann_tx_before):>7}")
    print(f"  genes (unique)        : {len(ann_g_before):>7}")
    print("== Annotation AFTER filtering (written to --out) ==")
    print(f"  rows                  : {n_rows_kept:>7}")
    print(f"  transcripts (unique)  : {len(kept_tx):>7}")
    print(f"  genes (unique)        : {len(kept_g):>7}")
    print("== Sanity checks ==")
    print(f"  annotation transcripts removed as not-in-FASTA : {len(removed):>7}")
    print(f"  FASTA transcripts WITHOUT an annotation row    : {len(fa_no_ann):>7}")
    print(f"  kept transcripts NOT in FASTA (must be 0)      : {len(orphans):>7}")
    print(f"  kept genes are a subset of FASTA genes         : {kept_g.issubset(fa_genes)}")

    ok = (len(orphans) == 0) and kept_g.issubset(fa_genes)
    print("\nRESULT:", "OK" if ok else "PROBLEM — check the flags above")
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
