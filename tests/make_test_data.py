#!/usr/bin/env python3
"""Generate the synthetic test dataset in tests/data/ (deterministic, standard library only).

* genome.fa  : one 52 kb chromosome with 40 two-exon genes on both strands
* genes.gtf  : GENCODE-style annotation (gene / transcript / exon, gene_id + gene_name)
* 6 samples  : 3 Normal + 3 Tumor, 27 bp single-end reads, unstranded library
* truth.tsv  : planted log2 fold changes; GENE01-05 up 8x and GENE06-10 down 8x in Tumor,
               GENE11-40 unchanged (biological noise added via a gamma-Poisson model)

tests/data_paired/ (same genome and genes): matched tumour/normal pairs from 3 patients,
paired-end 27 bp reads, each sample sequenced as 2 runs, plus a per-patient expression
effect on GENE11-20 so that the ~ patient + condition design is exercised.
"""
import gzip
import io
import math
import random
from pathlib import Path

SEED, N_GENES, READ_LEN, READS_PER_SAMPLE, DISPERSION = 20260209, 40, 27, 30000, 0.05
EXON1, INTRON, EXON2, GAP = 400, 300, 300, 300
OUT = Path(__file__).resolve().parent / "data"
COMP = str.maketrans("ACGT", "TGCA")

rng = random.Random(SEED)


def gz_text(path):
    """gzip writer with a fixed timestamp, so regenerated files are byte-identical."""
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", compresslevel=9, mtime=0))


def revcomp(s):
    return s.translate(COMP)[::-1]


def main():
    OUT.mkdir(exist_ok=True)
    genes, pos = [], GAP
    for g in range(1, N_GENES + 1):
        start, strand = pos + 1, "+" if g % 2 else "-"
        e1 = (start, start + EXON1 - 1)
        e2 = (e1[1] + INTRON + 1, e1[1] + INTRON + EXON2)
        genes.append(dict(id=f"ENSGTEST{g:08d}.1", name=f"GENE{g:02d}", strand=strand, exons=[e1, e2]))
        pos = e2[1] + GAP
    chrom = "".join(rng.choice("ACGT") for _ in range(pos + GAP))
    (OUT / "genome.fa").write_text(">chrT\n" + "\n".join(chrom[i:i + 60] for i in range(0, len(chrom), 60)) + "\n")

    with open(OUT / "genes.gtf", "w") as fh:
        for g in genes:
            s, e = g["exons"][0][0], g["exons"][1][1]
            tid = g["id"].replace("ENSG", "ENST")
            base = f'gene_id "{g["id"]}"; gene_name "{g["name"]}";'
            fh.write(f"chrT\ttest\tgene\t{s}\t{e}\t.\t{g['strand']}\t.\t{base} gene_type \"protein_coding\";\n")
            fh.write(f"chrT\ttest\ttranscript\t{s}\t{e}\t.\t{g['strand']}\t.\t{base} transcript_id \"{tid}\";\n")
            for n, (a, b) in enumerate(g["exons"], 1):
                fh.write(f"chrT\ttest\texon\t{a}\t{b}\t.\t{g['strand']}\t.\t{base} transcript_id \"{tid}\"; exon_number {n};\n")

    lfc = {g["name"]: (3.0 if i < 5 else -3.0 if i < 10 else 0.0) for i, g in enumerate(genes)}
    base_expr = {g["name"]: math.exp(rng.uniform(math.log(50), math.log(400))) for g in genes}
    with open(OUT / "truth.tsv", "w") as fh:
        fh.write("gene_id\tgene_name\tplanted_log2fc\n")
        for g in genes:
            fh.write(f"{g['id']}\t{g['name']}\t{lfc[g['name']]}\n")

    # mRNA sequence (sense strand) of each gene
    mrna = {}
    for g in genes:
        seq = "".join(chrom[a - 1:b] for a, b in g["exons"])
        mrna[g["name"]] = seq if g["strand"] == "+" else revcomp(seq)

    sheet = ["sample,condition,fastq_1"]
    for cond in ("Normal", "Tumor"):
        for rep in (1, 2, 3):
            sid = f"{cond[0]}{rep}"
            weights = []
            for g in genes:
                mu = base_expr[g["name"]] * (2 ** lfc[g["name"]] if cond == "Tumor" else 1.0)
                weights.append(rng.gammavariate(1 / DISPERSION, mu * DISPERSION))  # biological noise
            total = sum(weights)
            depth = int(READS_PER_SAMPLE * rng.uniform(0.8, 1.2))
            with gz_text(OUT / f"{sid}.fastq.gz") as fq:
                for r in range(depth):
                    x, acc = rng.random() * total, 0.0
                    for g, w in zip(genes, weights):
                        acc += w
                        if acc >= x:
                            break
                    seq = mrna[g["name"]]
                    p = rng.randrange(len(seq) - READ_LEN + 1)
                    read = seq[p:p + READ_LEN]
                    if rng.random() < 0.5:  # unstranded library: either strand
                        read = revcomp(read)
                    read = "".join(rng.choice("ACGT") if rng.random() < 0.002 else c for c in read)
                    fq.write(f"@{sid}.{r}\n{read}\n+\n{'I' * READ_LEN}\n")
            sheet.append(f"{sid},{cond},{sid}.fastq.gz")
    (OUT / "samplesheet.csv").write_text("\n".join(sheet) + "\n")
    print(f"wrote {OUT}")
    make_paired(genes, mrna, base_expr, lfc)


def pick(rng2, genes, weights, total):
    x, acc = rng2.random() * total, 0.0
    for g, w in zip(genes, weights):
        acc += w
        if acc >= x:
            return g
    return genes[-1]


def make_paired(genes, mrna, base_expr, lfc, frag=150, pairs_per_run=8000):
    out = OUT.parent / "data_paired"
    out.mkdir(exist_ok=True)
    rng2 = random.Random(SEED + 1)
    patient_effect = {p: {g["name"]: (rng2.uniform(-1.0, 1.0) if 10 <= i < 20 else 0.0) for i, g in enumerate(genes)}
                      for p in ("P1", "P2", "P3")}
    sheet = ["sample,patient,condition,run,fastq_1,fastq_2"]
    for p in ("P1", "P2", "P3"):
        for cond in ("Normal", "Tumor"):
            sid = f"{p}{cond[0]}"
            weights = []
            for g in genes:
                mu = base_expr[g["name"]] * 2 ** (patient_effect[p][g["name"]] + (lfc[g["name"]] if cond == "Tumor" else 0.0))
                weights.append(rng2.gammavariate(1 / DISPERSION, mu * DISPERSION))
            total = sum(weights)
            for run in (1, 2):
                rid = f"{sid}_run{run}"
                f1, f2 = out / f"{rid}_1.fastq.gz", out / f"{rid}_2.fastq.gz"
                with gz_text(f1) as q1, gz_text(f2) as q2:
                    for r in range(pairs_per_run):
                        g = pick(rng2, genes, weights, total)
                        seq = mrna[g["name"]]
                        a = rng2.randrange(len(seq) - frag + 1)
                        fr = seq[a:a + frag]
                        if rng2.random() < 0.5:  # unstranded: fragment from either strand
                            fr = revcomp(fr)
                        m1, m2 = fr[:READ_LEN], revcomp(fr)[:READ_LEN]
                        q1.write(f"@{rid}.{r}/1\n{m1}\n+\n{'I' * READ_LEN}\n")
                        q2.write(f"@{rid}.{r}/2\n{m2}\n+\n{'I' * READ_LEN}\n")
                sheet.append(f"{sid},{p},{cond},{rid},{f1.name},{f2.name}")
    (out / "samplesheet.csv").write_text("\n".join(sheet) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
