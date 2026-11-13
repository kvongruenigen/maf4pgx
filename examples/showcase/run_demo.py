#!/usr/bin/env python3
"""Run the public analysis CLIs offline on fictional example data."""

import argparse
import gzip
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parents[1] / "scripts"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="Defaults to a new temporary directory")
    args = parser.parse_args()
    output = args.output_dir or Path(tempfile.mkdtemp(prefix="maf4pgx-"))
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    prefix = output / "synthetic"
    variants = output / "synthetic.variants.tsv.gz"
    with (HERE / "synthetic.variants.tsv").open("rb") as source, gzip.open(variants, "wb") as target:
        shutil.copyfileobj(source, target)

    def run(script, *arguments):
        subprocess.run([sys.executable, "-B", str(SCRIPTS / script), *map(str, arguments)], check=True)

    for script in ("oncoprint.py", "mutation_summary.py"):
        run(script, variants, "--top", "2", "--out-prefix", prefix)
    run("audit_sample_join.py", variants, HERE / "synthetic.gdc_samples.tsv",
        HERE / "synthetic.progenetix_samples.tsv", "--output", output / "synthetic.sample_join.tsv")

    print(f"Open {output / 'synthetic.oncoprint.svg'} in a browser.")


if __name__ == "__main__":
    main()
