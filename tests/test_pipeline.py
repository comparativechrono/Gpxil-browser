"""The data pipeline, run on the toy genome and checked against the inputs.

    python -m pytest tests/test_pipeline.py
"""

from __future__ import annotations

import base64
import gzip
import json
import os
import re
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import build_fixture  # noqa: E402
from gpxiltools import open_gpxil  # noqa: E402

BIN_SIZES = (128, 512, 2048, 8192, 32768, 131072)


@pytest.fixture(scope="session")
def data():
    return build_fixture.build()


@pytest.fixture(scope="session")
def manifest(data):
    with open(os.path.join(data, "manifest.json")) as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def fasta(data):
    seqs, name = {}, None
    with open(os.path.join(build_fixture.INPUTS, "toy.fa")) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif line:
                seqs[name].append(line)
    names = {"1": "chr1", "2": "chr2", "MT": "chrM"}
    return {names.get(k, k): "".join(v) for k, v in seqs.items()}


@pytest.fixture(scope="session")
def expected(data):
    with open(os.path.join(build_fixture.INPUTS, "expected.json")) as fh:
        return json.load(fh)


def gz_json(data, rel):
    with gzip.open(os.path.join(data, rel), "rt") as fh:
        return json.load(fh)


def test_manifest(manifest, fasta):
    assert manifest["format"] == "gpxil-browser-data"
    assert manifest["assembly"] == "GRCh38"
    seq = manifest["sequence"]
    assert seq["width"] == 8192
    assert seq["overview_bins"] == list(BIN_SIZES)
    names = [c["name"] for c in seq["contigs"]]
    # chromosomes in karyotype order, then scaffolds longest first
    assert names == ["chr1", "chr2", "chrM", "KI270728.1", "GL000009.2"]
    for c in seq["contigs"]:
        assert c["length"] == len(fasta[c["name"]])
        assert re.fullmatch(r"[0-9a-f]{32}", c["md5"])
        assert c["refget"].startswith("SQ.")
        assert isinstance(c["ifd"], int) and c["ifd"] > 0
        chromosome = c["name"] in ("chr1", "chr2", "chrM")
        assert c["group"] == ("chromosome" if chromosome else "unplaced")
        assert c["file"] == (f"gpxil/{c['name']}.gpxil.tif" if chromosome else "gpxil/unplaced.gpxil.tif")
    chrm = next(c for c in seq["contigs"] if c["name"] == "chrM")
    assert set(chrm["aliases"]) >= {"MT", "M", "chrMT"}
    assert manifest["genes"]["count"] == 4


def test_gpxil_round_trip(data, manifest, fasta):
    for c in manifest["sequence"]["contigs"]:
        with open_gpxil(os.path.join(data, c["file"])) as g:
            contig = g[c["name"]]
            codes = contig.read("ref")
            want = fasta[c["name"]]
            assert len(codes) == len(want)
            assert contig.sequence().upper() == want.upper()
            # soft-masked (lower case) bases carry the 0x10 bit
            masked = [i for i, ch in enumerate(want) if ch.islower()]
            assert all(codes[i] & 0x10 for i in masked)
            assert not any(codes[i] & 0x10 for i in range(len(want)) if want[i].isupper())


def test_unplaced_file_holds_every_scaffold(data):
    with open_gpxil(os.path.join(data, "gpxil/unplaced.gpxil.tif")) as g:
        assert set(g.contigs) == {"KI270728.1", "GL000009.2"}


def expected_bins(seq: str, b: int) -> bytes:
    out = bytearray()
    for i in range(0, len(seq), b):
        part = seq[i : i + b]
        total = len(part)
        up = part.upper()
        acgt = sum(up.count(x) for x in "ACGT")
        gc = up.count("G") + up.count("C")
        mask = sum(ch.islower() for ch in part)
        out.append(round(255 * gc / acgt) if acgt else 0)
        out.append((round(15 * (total - acgt) / total) << 4) | round(15 * mask / total))
    return bytes(out)


def test_overview_bins(data, manifest, fasta):
    for c in manifest["sequence"]["contigs"]:
        with gzip.open(os.path.join(data, c["overview"]["file"])) as fh:
            blob = fh.read()
        for li, b in enumerate(BIN_SIZES):
            n = -(-c["length"] // b)
            off = c["overview"]["offsets"][li]
            assert blob[off : off + 2 * n] == expected_bins(fasta[c["name"]], b), (c["name"], b)


def test_karyotype(data, manifest, fasta):
    with open(os.path.join(data, "karyotype.json")) as fh:
        k = json.load(fh)
    assert k["bin"] == 131072
    assert set(k["bins"]) == {"chr1", "chr2", "chrM"}
    for name, b64 in k["bins"].items():
        assert base64.b64decode(b64) == expected_bins(fasta[name], 131072)


def test_genes(data, manifest, expected):
    chr1 = gz_json(data, "genes/chr1.json.gz")
    assert chr1["contig"] == "chr1"
    rows = {r[0]: r for r in chr1["genes"]}
    assert set(rows) == {"ENSG00000000001", "ENSG00000000002", "ENSG00000000004"}
    starts = [r[2] for r in chr1["genes"]]
    assert starts == sorted(starts)
    for gid, row in rows.items():
        e = expected["genes"][gid]
        _id, name, start, end, strand, biotype, hgnc_name, tx, xref, where = row
        assert (name, start, end, strand) == (e["name"], e["start"], e["end"], e["strand"])
        txid, txname, txbio, mane, ts, te, flat, cds, incomplete, *sec = tx
        exons = [[ts + flat[i], ts + flat[i + 1]] for i in range(0, len(flat), 2)]
        assert exons == e["exons"]
        assert (sec[0] if sec else []) == e["sec"]  # selenocysteine codons, only when there are some
        if e["cds"]:
            assert cds == [e["cds"][0], e["cds"][1], 0]
        else:
            assert cds is None
    pixl1 = rows["ENSG00000000001"]
    assert pixl1[5] == "protein_coding" and pixl1[6] == "pixel protein 1"
    assert pixl1[7][:4] == ["ENST00000000011", "PIXL1-201", "protein_coding", 1]  # the canonical, MANE Select transcript
    assert pixl1[8] == ["90001", "999001", "Q00001", "600001"]
    assert rows["ENSG00000000002"][7][3] == 0  # canonical but not MANE
    assert rows["ENSG00000000004"][5] == "lncRNA"
    mt = gz_json(data, "genes/chrM.json.gz")
    assert [r[1] for r in mt["genes"]] == ["MT-PIX"]
    assert not os.path.exists(os.path.join(data, "genes/chr2.json.gz"))  # no genes there
    files = {c["name"]: c["genes"] for c in manifest["sequence"]["contigs"]}
    assert files["chr1"] == "genes/chr1.json.gz" and files["chr2"] is None


def test_transcripts(data, manifest, expected):
    """Every transcript, one gzip member per gene, found through the gene rows."""
    files = {c["name"]: c["transcripts"] for c in manifest["sequence"]["contigs"]}
    assert files["chr1"] == "transcripts/chr1.bin" and files["chrM"] == "transcripts/chrM.bin" and files["chr2"] is None
    assert manifest["genes"]["transcripts"] == 5  # PIXL1 has two
    with open(os.path.join(data, files["chr1"]), "rb") as fh:
        blob = fh.read()
    rows = {r[0]: r for r in gz_json(data, "genes/chr1.json.gz")["genes"]}
    members = []
    for gid, row in rows.items():
        off, length, n = row[9]
        txs = json.loads(gzip.decompress(blob[off : off + length]))
        assert len(txs) == n
        members.append((off, length))
        g0 = row[2]
        canonical = txs[0]
        assert canonical[0] == row[7][0] and canonical[3] & 1  # canonical first, flagged
        e = expected["genes"][gid]
        assert [[g0 + canonical[6][i], g0 + canonical[6][i + 1]] for i in range(0, len(canonical[6]), 2)] == e["exons"]
        if e["cds"]:
            assert [g0 + canonical[7][0], g0 + canonical[7][1]] == e["cds"]
    # the members tile the file
    members.sort()
    assert members[0][0] == 0 and sum(n for _, n in members) == len(blob)
    pixl1 = json.loads(gzip.decompress(blob[rows["ENSG00000000001"][9][0] :][: rows["ENSG00000000001"][9][1]]))
    assert [t[0] for t in pixl1] == ["ENST00000000011", "ENST00000000012"]
    assert pixl1[0][3] == 1 | 2 | 8  # canonical, MANE Select, basic
    assert pixl1[1][3] == 8 and pixl1[1][7] is None  # basic, non-coding
    assert pixl1[1][6] == [0, 100, 10000, 10300]  # exons relative to the gene start
    assert pixl1[0][9] == [5070]  # the selenocysteine codon, relative too


def test_index_and_aliases(data):
    index = gz_json(data, "genes/index.json.gz")
    assert sorted(r[0] for r in index) == ["MT-PIX", "PIXL1", "PIXL1-AS1", "PIXL2"]  # GHOST's contig isn't in the FASTA
    assert ["PIXL1", "ENSG00000000001", "chr1", 100000, 110300, 1, "protein_coding"] in index
    aliases = dict((a, s) for a, s in gz_json(data, "genes/aliases.json.gz"))
    assert aliases["PIX1"] == ["PIXL1"]  # the withdrawn entry's alias is ignored
    assert aliases["OLDPIX"] == ["PIXL1"]
    assert aliases["PXL"] == ["PIXL1", "PIXL2"]


def test_cytobands(data):
    bands = gz_json(data, "cytobands.json.gz")
    assert set(bands) == {"chr1", "chr2"}
    assert bands["chr1"][0] == [0, 90000, "p12", "gneg"]
    assert [b[2] for b in bands["chr1"]] == ["p12", "p11", "q11", "q21.1", "q21.2"]
    assert bands["chr1"][-1][1] == 400000
