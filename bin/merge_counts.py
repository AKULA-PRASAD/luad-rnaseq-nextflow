#!/usr/bin/env python3
"""Merge STAR ReadsPerGene.out.tab files into one gene x sample count matrix.

STAR writes three count columns per gene:
  column 2 = unstranded, column 3 = forward-stranded, column 4 = reverse-stranded.

With --strandedness auto, the library type is inferred from the gene-assigned
totals of columns 3 and 4 over all samples:
  * one strand holds >= 80 % of stranded counts -> that strand
  * neither strand above 60 %                  -> unstranded
  * in between                                 -> unstranded, with a warning
The per-sample N_noFeature values of all three columns (the check used in the
original manual analysis) are written to the report as well.

Gene names are taken from the GTF (gene_name attribute) so results are readable.
Standard library only.
"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

COLUMNS = {"unstranded": 1, "forward": 2, "reverse": 3}  # index into tab-split row


def read_tab(path):
    summary, genes = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 4:
                sys.exit(f"{path}: expected 4 columns, got {len(f)}")
            vals = [int(x) for x in f[1:]]
            (summary if f[0].startswith("N_") else genes)[f[0]] = vals
    return summary, genes


def gene_names(gtf):
    opener = gzip.open if str(gtf).endswith(".gz") else open
    names, pat_id, pat_name = {}, re.compile(r'gene_id "([^"]+)"'), re.compile(r'gene_name "([^"]+)"')
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid, gname = pat_id.search(f[8]), pat_name.search(f[8])
            if gid:
                names[gid.group(1)] = gname.group(1) if gname else ""
    return names


def infer(totals):
    fwd, rev = totals["forward"], totals["reverse"]
    if fwd + rev == 0:
        sys.exit("No gene-assigned reads in any sample: check the GTF and the alignment.")
    share = max(fwd, rev) / (fwd + rev)
    if share >= 0.8:
        return ("forward" if fwd > rev else "reverse"), share, ""
    if share <= 0.6:
        return "unstranded", share, ""
    return "unstranded", share, "WARNING: strand evidence is ambiguous; using unstranded counts"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tabs", nargs="+", type=Path)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--strandedness", default="auto", choices=["auto", *COLUMNS])
    ap.add_argument("--out", default="counts_matrix.csv")
    ap.add_argument("--report", default="strandedness.tsv")
    a = ap.parse_args()

    samples = {}
    for t in sorted(a.tabs):
        sid = t.name.replace(".ReadsPerGene.out.tab", "")
        if sid in samples:
            sys.exit(f"Duplicate sample id: {sid}")
        samples[sid] = read_tab(t)

    gene_order = list(next(iter(samples.values()))[1])
    for sid, (_, genes) in samples.items():
        if list(genes) != gene_order:
            sys.exit(f"{sid}: gene list differs from the other samples (different index or GTF?)")

    totals = {k: sum(sum(g[i - 1] for g in genes.values()) for _, genes in samples.values()) for k, i in COLUMNS.items()}
    inferred, share, warn = infer(totals)
    chosen = inferred if a.strandedness == "auto" else a.strandedness
    if warn and a.strandedness == "auto":
        print(warn, file=sys.stderr)

    with open(a.report, "w") as fh:
        fh.write(f"# requested\t{a.strandedness}\n# inferred\t{inferred}\n# majority_strand_share\t{share:.3f}\n# used\t{chosen}\n")
        if warn:
            fh.write(f"# {warn}\n")
        fh.write("sample\tN_noFeature_unstranded\tN_noFeature_forward\tN_noFeature_reverse\tassigned_unstranded\tassigned_forward\tassigned_reverse\n")
        for sid, (summ, genes) in samples.items():
            nf = summ.get("N_noFeature", [0, 0, 0])
            assigned = [sum(g[i] for g in genes.values()) for i in range(3)]
            fh.write("\t".join([sid, *map(str, nf), *map(str, assigned)]) + "\n")

    names = gene_names(a.gtf)
    col = COLUMNS[chosen] - 1
    ids = list(samples)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene_id", "gene_name", *ids])
        for g in gene_order:
            w.writerow([g, names.get(g, ""), *(samples[s][1][g][col] for s in ids)])
    print(f"{len(gene_order)} genes x {len(ids)} samples; strandedness used: {chosen} (inferred {inferred}, share {share:.2f})")


if __name__ == "__main__":
    main()
