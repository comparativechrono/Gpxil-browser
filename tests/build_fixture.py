"""Build the browser's data for the toy genome (tests/fixture.py) into tests/.fixture/.

    python tests/build_fixture.py

Afterwards tests/.fixture/site/data/grch38/ holds what the pipeline makes for
the real genome, small enough to check exactly; the Python and JavaScript
tests both read it.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIX = os.path.join(HERE, ".fixture")
INPUTS = os.path.join(FIX, "inputs")
DATA = os.path.join(FIX, "site", "data", "grch38")


def build() -> str:
    shutil.rmtree(FIX, ignore_errors=True)
    subprocess.run([sys.executable, os.path.join(HERE, "fixture.py"), "--out", INPUTS], check=True)
    subprocess.run(
        [
            sys.executable,
            "-I",
            os.path.join(ROOT, "pipeline", "build.py"),
            "--out",
            DATA,
            "--cache",
            os.path.join(FIX, "cache"),
            "--fasta",
            os.path.join(INPUTS, "toy.fa"),
            "--gtf",
            os.path.join(INPUTS, "toy.gtf"),
            "--hgnc",
            os.path.join(INPUTS, "hgnc.txt"),
            "--bands",
            os.path.join(INPUTS, "bands.json"),
        ],
        check=True,
    )
    return DATA


if __name__ == "__main__":
    print(build())
