"""Build the data the Gpxil browser serves.

    python -I pipeline/build.py --out site/data/grch38            # downloads everything
    python -I pipeline/build.py --out site/data/grch38 --fasta local.fa.gz --gtf local.gtf.gz ...

Inputs default to Ensembl release 116 (GRCh38 soft-masked primary assembly
and gene set), the HGNC complete set, and Ensembl's cytobands. Every input
may be a local path or a URL; downloads are cached in --cache.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import gpxiltools  # noqa: E402

from common import fetch, log, write_json  # noqa: E402
from genes import build_bands, build_genes  # noqa: E402
from sequence import build_sequence  # noqa: E402

ENSEMBL_RELEASE = 116
ENSEMBL = f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}"
DEFAULTS = {
    "fasta": f"{ENSEMBL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz",
    "gtf": f"{ENSEMBL}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz",
    "hgnc": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "bands": "https://rest.ensembl.org/info/assembly/homo_sapiens?bands=1;content-type=application/json",
}
DATA_VERSION = 2  # bump when the layout of the data files changes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output directory (served as data/grch38/)")
    ap.add_argument("--cache", default=".cache", help="download cache directory")
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k}", default=v, help=f"path or URL (default: {v})")
    ap.add_argument("--only", help="comma-separated contigs, for a quick test build")
    ap.add_argument("--skip-sequence", action="store_true", help="reuse gpxil/ and overview/ from an earlier build")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    # the FASTA (about 1 GB) is only needed when the sequence is encoded
    paths = {k: fetch(getattr(a, k), a.cache) for k in DEFAULTS if not (k == "fasta" and a.skip_sequence)}
    only = set(a.only.split(",")) if a.only else None
    manifest_path = os.path.join(a.out, "manifest.json")

    if a.skip_sequence:
        import json

        with open(manifest_path) as fh:
            seq = json.load(fh)["sequence"]
        with open(os.path.join(a.out, "karyotype.json")) as fh:
            seq["karyotype"] = json.load(fh)
    else:
        log("encoding the reference ...")
        seq = build_sequence(paths["fasta"], a.out, only=only, source_name=os.path.basename(DEFAULTS["fasta"]))
    log("gene annotation ...")
    genes = build_genes(paths["gtf"], paths["hgnc"], a.out, [c["name"] for c in seq["contigs"]])
    for c in seq["contigs"]:
        g = genes["per_contig"].get(c["name"])
        c["genes"] = g["file"] if g else None
        c["gene_count"] = g["count"] if g else 0
        c["transcripts"] = g["transcripts"] if g else None
    bands = build_bands(paths["bands"], a.out)
    write_json(os.path.join(a.out, "karyotype.json"), seq["karyotype"])

    manifest = {
        "format": "gpxil-browser-data",
        "data_version": DATA_VERSION,
        "assembly": "GRCh38",
        "created": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "generator": f"gpxil-browser pipeline, gpxiltools {gpxiltools.__version__}",
        "sources": {
            "sequence": f"Ensembl release {ENSEMBL_RELEASE}: {os.path.basename(DEFAULTS['fasta'])}",
            "genes": f"Ensembl release {ENSEMBL_RELEASE}: {os.path.basename(DEFAULTS['gtf'])} (all transcripts)",
            "names": "HGNC complete set (approved symbols, names, aliases, cross-references)",
            "bands": f"Ensembl REST API, info/assembly ({bands['assembly']})",
            "rest": "https://rest.ensembl.org",
        },
        "sequence": {k: v for k, v in seq.items() if k != "karyotype"},
        "karyotype": "karyotype.json",
        "cytobands": bands["file"],
        "genes": {"index": genes["index"], "aliases": genes["aliases"], "count": genes["count"], "transcripts": genes["transcripts"]},
    }
    total = 0
    for root, _, files in os.walk(a.out):
        total += sum(os.path.getsize(os.path.join(root, f)) for f in files)
    manifest["bytes"] = total
    write_json(manifest_path, manifest)
    log(f"wrote {manifest_path}: {len(seq['contigs'])} contigs, {total / 1e6:,.1f} MB of data")
    return 0


if __name__ == "__main__":
    sys.exit(main())
