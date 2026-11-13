#!/usr/bin/env python3
"""Build the offline copy-number report from fictional inputs."""

import argparse
from pathlib import Path
import subprocess
import sys
from tempfile import mkdtemp

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "copy_number"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Defaults to a new temporary directory")
    args = parser.parse_args()
    output = args.output or Path(mkdtemp(prefix="maf4pgx-")) / "synthetic.copy-number.html"
    subprocess.run([
        sys.executable, "-B", str(HERE.parents[1] / "scripts" / "showcase_report.py"),
        "--variants", str(FIXTURE / "synthetic.variants.tsv"),
        "--pilot-prefix", str(FIXTURE / "synthetic"),
        "--segment-cache", str(FIXTURE / "segments"),
        "--label", "Synthetic cohort · Fictional data", "--synthetic",
        "--output", str(output),
    ], check=True)


if __name__ == "__main__":
    main()
