#!/usr/bin/env python3
"""Check the offline report against its fictional evidence."""

import json
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read_report(path):
    html = path.read_text(encoding="utf-8")
    return json.loads(re.search(r'id="cohort-data">(.*?)</script>', html).group(1))


class CopyNumberReportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        output = TemporaryDirectory()
        cls.addClassCleanup(output.cleanup)
        cls.html = Path(output.name) / "report.html"
        subprocess.run([sys.executable, "-B", str(ROOT / "examples/showcase/run_report.py"),
                        "--output", str(cls.html)], check=True)
        cls.report = read_report(cls.html)

    def test_preview_matches_rebuild(self):
        self.assertEqual(self.html.read_bytes(),
                         (ROOT / "examples/showcase/synthetic.copy-number.html").read_bytes())

    def test_interval_resolution_and_evidence(self):
        tp53 = self.report["records"]["TP53"]
        apc = self.report["records"]["APC"]
        self.assertEqual(len(self.report["samples"]), 8)
        exact = tp53["SYNTHETIC-01"]["gdc"]
        self.assertEqual(exact["segments"], [{"start": 7661779, "end": 7687538, "mean": -0.15}])
        self.assertEqual(exact["mean"], -0.15)
        self.assertEqual(tp53["SYNTHETIC-04"]["gdc"]["reason"], "Partial locus coverage")
        self.assertEqual(tp53["SYNTHETIC-05"]["gdc"]["reason"], "Multiple overlapping segments")
        self.assertEqual(len(tp53["SYNTHETIC-05"]["gdc"]["segments"]), 2)
        self.assertEqual(apc["SYNTHETIC-05"]["gdc"]["reason"], "No overlapping segment")
        self.assertEqual(tp53["SYNTHETIC-07"]["gdc"]["reason"], "Segment file missing from cache")
        self.assertEqual(apc["SYNTHETIC-08"]["gdc"]["reason"], "Ambiguous files")
        self.assertIsNone(tp53["SYNTHETIC-04"]["gdc"]["mean"])
        self.assertIsNone(tp53["SYNTHETIC-05"]["gdc"]["mean"])

    def test_join_and_deduplication(self):
        tp53 = self.report["records"]["TP53"]
        duplicate = tp53["SYNTHETIC-01"]
        self.assertEqual(duplicate["protein_count"], 1)
        self.assertEqual(len(duplicate["variants"]), 1)
        self.assertEqual(len(duplicate["variants"][0]["files"]), 2)
        self.assertEqual(duplicate["pgx"]["intervals"][0]["start"], 7661779)
        self.assertEqual(tp53["SYNTHETIC-03"]["protein_count"], 0)
        self.assertEqual(tp53["SYNTHETIC-04"]["pgx"]["join"], "ambiguous_progenetix_match")
        self.assertEqual(tp53["SYNTHETIC-05"]["pgx"]["join"], "no_progenetix_match")
        self.assertEqual(tp53["SYNTHETIC-06"]["pgx"]["reason"], "Synthetic query failure")
        self.assertEqual(tp53["SYNTHETIC-08"]["pgx"]["join"], "no_gdc_uuid")


if __name__ == "__main__":
    unittest.main()
