"""Reference sequence -> gpxil files, overview bins and karyotype bins.

Each primary chromosome becomes its own reference-only gpxil file (one
`ref` page, 8,192 bases per row, Deflate, strips of 256 KiB) so the browser
can read any interval with a couple of HTTP range requests. The unplaced
and unlocalised scaffolds share one multi-contig file.

Overview bins summarise the sequence for zoomed-out views: for bins of
128, 512, ... 131,072 bases, two bytes per bin: the GC fraction of the
A/C/G/T bases (0-255), and the fraction of other bases (N) and of
soft-masked bases (4 bits each).
"""

from __future__ import annotations

import base64
import gzip
import json
import os

import numpy as np
import tifffile
from gpxiltools.container import GpxilWriter
from gpxiltools.fasta import iter_fasta, sequence_digests

from common import PRIMARY, aliases_for, browser_name, log, write_bytes

WIDTH = 8192
BIN_SIZES = (128, 512, 2048, 8192, 32768, 131072)
KARYOTYPE_BIN = 131072

_ACGT = np.zeros(256, dtype=bool)
_GC = np.zeros(256, dtype=bool)
for _c in (1, 2, 4, 8):
    _ACGT[_c] = _ACGT[_c | 0x10] = True
for _c in (2, 4):
    _GC[_c] = _GC[_c | 0x10] = True
_MASK = (np.arange(256) & 0x10) != 0


def bin_counts(codes: np.ndarray, b: int = BIN_SIZES[0], chunk: int = 1 << 24) -> dict[str, np.ndarray]:
    """Per-bin counts of all, A/C/G/T, G/C and soft-masked bases (bins of ``b`` bases)."""
    n = codes.shape[0]
    nb = -(-n // b)
    out = {k: np.zeros(nb, dtype=np.int64) for k in ("total", "acgt", "gc", "mask")}
    out["total"][:] = b
    out["total"][-1] = n - (nb - 1) * b
    chunk -= chunk % b
    for s in range(0, n, chunk):
        seg = codes[s : s + chunk]
        m = seg.shape[0]
        if m % b:
            seg = np.concatenate([seg, np.zeros(b - m % b, dtype=np.uint8)])
        i0 = s // b
        i1 = i0 + seg.shape[0] // b
        out["acgt"][i0:i1] = _ACGT[seg].reshape(-1, b).sum(1)
        out["gc"][i0:i1] = _GC[seg].reshape(-1, b).sum(1)
        out["mask"][i0:i1] = _MASK[seg].reshape(-1, b).sum(1)
    return out


def coarsen(counts: dict[str, np.ndarray], factor: int = 4) -> dict[str, np.ndarray]:
    out = {}
    for k, v in counts.items():
        pad = (-v.shape[0]) % factor
        if pad:
            v = np.concatenate([v, np.zeros(pad, dtype=v.dtype)])
        out[k] = v.reshape(-1, factor).sum(1)
    return out


def encode_bins(c: dict[str, np.ndarray]) -> bytes:
    """Two bytes per bin: GC fraction of A/C/G/T (0-255); N fraction << 4 | soft-mask fraction (0-15 each)."""
    acgt = c["acgt"].astype(np.float64)
    total = np.maximum(c["total"], 1).astype(np.float64)
    gc8 = np.where(acgt > 0, np.rint(255 * c["gc"] / np.maximum(acgt, 1)), 0).astype(np.uint8)
    n4 = np.rint(15 * (total - acgt) / total).astype(np.uint8)
    m4 = np.rint(15 * c["mask"] / total).astype(np.uint8)
    out = np.empty(2 * acgt.shape[0], dtype=np.uint8)
    out[0::2] = gc8
    out[1::2] = (n4 << 4) | m4
    return out.tobytes()


def overview(codes: np.ndarray) -> tuple[bytes, list[int], bytes]:
    """All overview levels of one contig, their offsets, and the karyotype level."""
    counts = bin_counts(codes, BIN_SIZES[0])
    parts, offsets, pos, karyo = [], [], 0, b""
    for i, b in enumerate(BIN_SIZES):
        if i:
            counts = coarsen(counts, b // BIN_SIZES[i - 1])
        blob = encode_bins(counts)
        offsets.append(pos)
        parts.append(blob)
        pos += len(blob)
        if b == KARYOTYPE_BIN:
            karyo = blob
    return b"".join(parts), offsets, karyo


def build_sequence(fasta: str, out_dir: str, *, only: set[str] | None = None, source_name: str | None = None) -> dict:
    """Write gpxil/, overview/ and the karyotype bins; return manifest pieces."""
    gdir, odir = os.path.join(out_dir, "gpxil"), os.path.join(out_dir, "overview")
    os.makedirs(gdir, exist_ok=True)
    os.makedirs(odir, exist_ok=True)
    reference = {"name": "GRCh38", "source": source_name or os.path.basename(fasta)}
    contigs: list[dict] = []
    karyotype: dict[str, str] = {}
    unplaced = None
    unplaced_ov: list[bytes] = []
    unplaced_pos = 0
    for rec in iter_fasta(fasta):
        if only is not None and rec.name not in only and browser_name(rec.name) not in only:
            continue
        if rec.length == 0:
            continue
        name = browser_name(rec.name)
        primary = rec.name in PRIMARY
        digests = sequence_digests(rec.codes)
        extra = {"aliases": aliases_for(rec.name), **digests}
        if primary:
            fname = f"gpxil/{name}.gpxil.tif"
            with GpxilWriter(os.path.join(out_dir, fname), width=WIDTH, reference=reference, bigtiff=False) as w:
                w.add_contig(name, {"ref": rec.codes}, **extra)
        else:
            fname = "gpxil/unplaced.gpxil.tif"
            if unplaced is None:
                unplaced = GpxilWriter(os.path.join(out_dir, fname), width=WIDTH, reference=reference, bigtiff=False)
            unplaced.add_contig(name, {"ref": rec.codes}, **extra)
        blob, offsets, karyo = overview(rec.codes)
        if primary:
            ofile = f"overview/{name}.ov.gz"
            write_bytes(os.path.join(out_dir, ofile), gzip.compress(blob, 9, mtime=0))
            base = 0
            karyotype[name] = base64.b64encode(karyo).decode("ascii")
        else:
            ofile = "overview/unplaced.ov.gz"
            base = unplaced_pos
            unplaced_ov.append(blob)
            unplaced_pos += len(blob)
        contigs.append(
            {
                "name": name,
                "ensembl": rec.name,
                "aliases": extra["aliases"],
                "length": rec.length,
                "group": "chromosome" if primary else "unplaced",
                "file": fname,
                "md5": digests["md5"],
                "refget": digests["refget"],
                "overview": {"file": ofile, "offsets": [base + o for o in offsets]},
            }
        )
        log(f"  {name}: {rec.length:,} bp -> {fname}")
    if unplaced is not None:
        unplaced.close()
        write_bytes(os.path.join(out_dir, "overview/unplaced.ov.gz"), gzip.compress(b"".join(unplaced_ov), 9, mtime=0))
    # where each contig's page (IFD) starts, so the browser can read it without walking the IFD chain
    files = sorted({c["file"] for c in contigs})
    for f in files:
        with tifffile.TiffFile(os.path.join(out_dir, f)) as tf:
            meta = json.loads(tf.pages[0].description)
            ifd = {e["name"]: tf.pages[e["pages"]["ref"]].offset for e in meta["contigs"]}
        for c in contigs:
            if c["file"] == f:
                c["ifd"] = int(ifd[c["name"]])
                c["bytes"] = os.path.getsize(os.path.join(out_dir, f)) if c["group"] == "chromosome" else None
    order = {n: i for i, n in enumerate(PRIMARY)}
    contigs.sort(key=lambda c: (0, order[c["ensembl"]]) if c["group"] == "chromosome" else (1, -c["length"]))
    return {
        "width": WIDTH,
        "overview_bins": list(BIN_SIZES),
        "contigs": contigs,
        "karyotype": {"bin": KARYOTYPE_BIN, "bins": karyotype},
    }
