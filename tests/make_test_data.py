#!/usr/bin/env python3
"""Generate the synthetic test dataset in tests/data/ (deterministic, standard library only).

* genome.fa  : one 52 kb chromosome with 40 two-exon genes on both strands
* genes.gtf  : GENCODE-style annotation (gene / transcript / exon, gene_id + gene_name)
* 6 samples  : 3 Normal + 3 Tumor, 27 bp single-end reads, unstranded library
* truth.tsv  : planted log2 fold changes; GENE01-05 up 8x and GENE06-10 down 8x in Tumor,
               GENE11-40 unchanged (biological noise added via a gamma-Poisson model)
"""
import gzip
import math
import random
from pathlib import Path

SEED, N_GENES, READ_LEN, READS_PER_SAMPLE, DISPERSION = 20260209, 40, 27, 30000, 0.05
EXON1, INTRON, EXON2, GAP = 400, 300, 300, 300
OUT = Path(__file__).resolve().parent / "data"
COMP = str.maketrans("ACGT", "TGCA")

rng = random.Random(SEED)


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
            with gzip.open(OUT / f"{sid}.fastq.gz", "wt", compresslevel=9) as fq:
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


if __name__ == "__main__":
    main()
