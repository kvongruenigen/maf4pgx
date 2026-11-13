"""Smoke test the documented offline demo and its illustrative outputs."""

import csv
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "examples" / "showcase"


class ShowcaseTest(unittest.TestCase):
    def test_offline_quickstart(self):
        result = subprocess.run(
            [sys.executable, "examples/showcase/run_demo.py"], cwd=ROOT,
            check=True, capture_output=True, text=True,
        )
        plot = Path(result.stdout.splitlines()[-1].removeprefix("Open ").removesuffix(" in a browser."))
        self.addCleanup(self.remove_demo, plot.parent)
        self.assertEqual(plot.read_bytes(), (DEMO / plot.name).read_bytes())

        def read(name):
            with (plot.parent / f"synthetic.{name}.tsv").open(newline="") as source:
                return list(csv.DictReader(source, delimiter="\t"))

        genes = {row["gene"]: row for row in read("gene_counts")}
        self.assertEqual(genes["GENE_A"]["percent_of_samples_in_export"], "50.00")
        self.assertEqual(genes["GENE_A"]["distinct_sample_variant_calls"], "4")
        sample = [row for row in read("sample_calls") if row["sample_barcode"] == "SYNTHETIC-01"]
        self.assertEqual([row["distinct_calls"] for row in sample], ["1"])

        recurrent = read("recurrent_variants")
        self.assertEqual(len(recurrent), 1)
        self.assertEqual(
            (recurrent[0]["genes"], recurrent[0]["samples_with_call"], recurrent[0]["samples_in_export"]),
            ("GENE_A", "2", "8"),
        )
        self.assertEqual(read("gene_pairs"), [{
            "gene_a": "GENE_A", "gene_b": "GENE_B", "both": "1",
            "a_only": "3", "b_only": "2", "neither": "2", "samples_in_export": "8",
        }])
        self.assertEqual(
            [row["status"] for row in read("sample_join")],
            ["matched", "no_progenetix_match", "no_gdc_uuid", "ambiguous_gdc_uuid",
             "ambiguous_progenetix_match", "matched", "matched", "matched"],
        )

    @staticmethod
    def remove_demo(output):
        # Guard against deleting any path outside a demo-created directory.
        if output.name.startswith("maf4pgx-"):
            for path in output.glob("synthetic.*"):
                path.unlink()
            output.rmdir()


if __name__ == "__main__":
    unittest.main()
