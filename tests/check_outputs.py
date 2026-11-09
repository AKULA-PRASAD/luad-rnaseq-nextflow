#!/usr/bin/env python3
"""Check a test-profile run against the planted truth. Exit code 1 on any failure.

Usage: python3 tests/check_outputs.py [results_dir] [samplesheet]
       defaults: results_test, tests/data/samplesheet.csv
       paired:   python3 tests/check_outputs.py results_test_paired tests/data_paired/samplesheet.csv
"""
import csv
import sys
from pathlib import Path

res_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "results_test")
here = Path(__file__).resolve().parent
sheet_file = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "data" / "samplesheet.csv"
truth_file = here / "data" / "truth.tsv"
sheet = list(csv.DictReader(open(sheet_file)))
n_samples = len({r["sample"] for r in sheet})
paired_design = "patient" in sheet[0]
fails = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)


truth = {r["gene_id"]: float(r["planted_log2fc"]) for r in csv.DictReader(open(truth_file), delimiter="\t")}

for f in ["fastqc", "star", "counts/counts_matrix.csv", "counts/strandedness.tsv", "deseq2/deseq2_results.csv",
          "deseq2/significant_genes.csv", "deseq2/volcano.pdf", "deseq2/pca.pdf", "multiqc/multiqc_report.html"]:
    check((res_dir / f).exists(), f"output exists: {f}")
if fails:
    sys.exit(1)

counts = list(csv.DictReader(open(res_dir / "counts/counts_matrix.csv")))
check(len(counts) == len(truth), f"count matrix has {len(counts)} genes (expected {len(truth)})")
check(len(counts[0]) - 2 == n_samples, f"count matrix has {len(counts[0]) - 2} samples (expected {n_samples}, runs merged per sample)")
check(all(r["gene_name"] for r in counts), "every gene has a gene_name from the GTF")
assigned = sum(int(v) for r in counts for k, v in r.items() if k not in ("gene_id", "gene_name"))
check(assigned > 50000, f"{assigned} reads (or read pairs) assigned to genes")

used = [l.split("\t")[1].strip() for l in open(res_dir / "counts/strandedness.tsv") if l.startswith("# used")]
check(used == ["unstranded"], f"strandedness detected as {used} (library is unstranded)")

design = [l.split("\t", 1)[1].strip() for l in open(res_dir / "deseq2/summary.txt") if l.startswith("design")]
want = "~batch + condition" if paired_design else "~condition"
check(design == [want], f"DESeq2 design is {design} (expected {want})")

res = {r["gene_id"]: r for r in csv.DictReader(open(res_dir / "deseq2/deseq2_results.csv"))}
sig = {g for g, r in res.items() if r["padj"] not in ("", "NA") and float(r["padj"]) < 0.05}
planted = {g for g, l in truth.items() if l != 0}
check(planted <= sig, f"all {len(planted)} planted genes significant ({len(planted & sig)} found)")
check(all((float(res[g]["log2FoldChange"]) > 0) == (truth[g] > 0) for g in planted), "planted genes have the correct direction")
check(all(abs(float(res[g]["log2FoldChange"]) - truth[g]) < 1.0 for g in planted), "planted log2FC recovered within 1.0")
false_pos = sig - planted
check(len(false_pos) <= 2, f"{len(false_pos)} false positives among {len(truth) - len(planted)} unchanged genes (allowed <= 2)")

print(f"\n{len(fails)} failure(s)")
sys.exit(1 if fails else 0)
