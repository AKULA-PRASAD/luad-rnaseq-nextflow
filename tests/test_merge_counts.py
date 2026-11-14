#!/usr/bin/env python3
"""Unit tests for bin/merge_counts.py (standard library only).

Run: python3 -m unittest discover -s tests -p 'test_*.py' -v
"""
import gzip
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "bin" / "merge_counts.py"
spec = importlib.util.spec_from_file_location("merge_counts", SCRIPT)
mc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mc)

SUMMARY = "N_unmapped\t5\t5\t5\nN_multimapping\t3\t3\t3\nN_noFeature\t10\t60\t60\nN_ambiguous\t1\t0\t0\n"
GTF = ('chr1\tt\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "ALPHA";\n'
       'chr1\tt\texon\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "ALPHA";\n'
       'chr1\tt\tgene\t200\t300\t.\t-\t.\tgene_id "G2"; gene_name "BETA";\n')


def write_tab(folder, sample, genes):
    p = Path(folder) / f"{sample}.ReadsPerGene.out.tab"
    p.write_text(SUMMARY + "".join(f"{g}\t{u}\t{f}\t{r}\n" for g, (u, f, r) in genes.items()))
    return p


class InferStrandedness(unittest.TestCase):
    def test_unstranded_when_strands_balanced(self):
        self.assertEqual(mc.infer({"forward": 50, "reverse": 50})[0], "unstranded")

    def test_reverse_when_reverse_dominates(self):
        self.assertEqual(mc.infer({"forward": 10, "reverse": 90})[0], "reverse")

    def test_forward_at_threshold(self):
        self.assertEqual(mc.infer({"forward": 80, "reverse": 20})[0], "forward")

    def test_ambiguous_warns_and_uses_unstranded(self):
        kind, share, warn = mc.infer({"forward": 70, "reverse": 30})
        self.assertEqual(kind, "unstranded")
        self.assertAlmostEqual(share, 0.7)
        self.assertIn("ambiguous", warn)

    def test_no_reads_exits(self):
        with self.assertRaises(SystemExit):
            mc.infer({"forward": 0, "reverse": 0})


class ReadInputs(unittest.TestCase):
    def test_gene_names_from_plain_and_gz_gtf(self):
        with tempfile.TemporaryDirectory() as d:
            plain, gz = Path(d) / "a.gtf", Path(d) / "a.gtf.gz"
            plain.write_text(GTF)
            with gzip.open(gz, "wt") as fh:
                fh.write(GTF)
            for p in (plain, gz):
                self.assertEqual(mc.gene_names(p), {"G1": "ALPHA", "G2": "BETA"})

    def test_read_tab_rejects_wrong_column_count(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad.ReadsPerGene.out.tab"
            p.write_text("G1\t1\t2\n")
            with self.assertRaises(SystemExit):
                mc.read_tab(p)


class EndToEnd(unittest.TestCase):
    def run_script(self, d, tabs, *extra):
        gtf = Path(d) / "genes.gtf"
        gtf.write_text(GTF)
        out, rep = Path(d) / "m.csv", Path(d) / "s.tsv"
        cmd = [sys.executable, str(SCRIPT), "--gtf", str(gtf), "--out", str(out), "--report", str(rep), *extra, *map(str, tabs)]
        return subprocess.run(cmd, capture_output=True, text=True), out, rep

    def test_matrix_uses_unstranded_column_and_names(self):
        with tempfile.TemporaryDirectory() as d:
            tabs = [write_tab(d, "S1", {"G1": (100, 50, 50), "G2": (20, 10, 10)}),
                    write_tab(d, "S2", {"G1": (80, 41, 39), "G2": (30, 15, 15)})]
            proc, out, rep = self.run_script(d, tabs)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(out.read_text().splitlines(),
                             ["gene_id,gene_name,S1,S2", "G1,ALPHA,100,80", "G2,BETA,20,30"])
            self.assertIn("# used\tunstranded", rep.read_text())

    def test_reverse_library_uses_column_4(self):
        with tempfile.TemporaryDirectory() as d:
            tabs = [write_tab(d, "S1", {"G1": (100, 5, 95), "G2": (20, 1, 19)})]
            proc, out, _ = self.run_script(d, tabs)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("G1,ALPHA,95", out.read_text())

    def test_explicit_strandedness_overrides_inference(self):
        with tempfile.TemporaryDirectory() as d:
            tabs = [write_tab(d, "S1", {"G1": (100, 5, 95), "G2": (20, 1, 19)})]
            proc, out, rep = self.run_script(d, tabs, "--strandedness", "unstranded")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("G1,ALPHA,100", out.read_text())
            self.assertIn("# inferred\treverse", rep.read_text())

    def test_samples_with_different_gene_lists_fail(self):
        with tempfile.TemporaryDirectory() as d:
            tabs = [write_tab(d, "S1", {"G1": (1, 1, 1), "G2": (1, 1, 1)}),
                    write_tab(d, "S2", {"G1": (1, 1, 1)})]
            proc, _, _ = self.run_script(d, tabs)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("gene list differs", proc.stderr)


if __name__ == "__main__":
    unittest.main()
