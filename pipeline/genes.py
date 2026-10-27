"""Gene annotation for the browser: Ensembl GTF + HGNC names, aliases and cross-references.

Writes:

- ``genes/index.json.gz``: every gene, for search: [symbol, gene id, contig,
  start (0-based), end, strand, biotype].
- ``genes/aliases.json.gz``: HGNC alias and previous symbols -> approved symbol.
- ``genes/<contig>.json.gz``: the genes of one contig, each with its canonical
  transcript (MANE Select where there is one, otherwise Ensembl canonical),
  its HGNC name and cross-references, and where its other transcripts are.
  A coding transcript's CDS is [start, end, frame, Ensembl protein ID], the
  end including the stop codon.
- ``transcripts/<contig>.bin``: every transcript of every gene on the contig,
  one gzip member per gene, so the browser reads a gene's transcripts with a
  single HTTP range request (each gene row gives offset, length and count).

Coordinates are 0-based and end-exclusive throughout.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import os
import re
from collections import defaultdict

from common import browser_name, log, write_bytes, write_json

_ATTR = re.compile(r'(\w+) "([^"]*)"')
KINDS = {"gene", "transcript", "exon", "CDS", "stop_codon", "Selenocysteine"}
INCOMPLETE = ("cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF")
# transcript flags in the transcripts files
CANONICAL, MANE_SELECT, MANE_CLINICAL, BASIC = 1, 2, 4, 8


def _attrs(text: str) -> dict:
    out: dict = {}
    for k, v in _ATTR.findall(text):
        if k == "tag":
            out.setdefault("tags", set()).add(v)
        else:
            out[k] = v
    return out


def _attr(text: str, key: str) -> str | None:
    """One attribute of a GTF line, without parsing them all (exon and CDS lines are millions)."""
    i = text.find(key + ' "')
    if i < 0:
        return None
    i += len(key) + 2
    return text[i : text.index('"', i)]


def iter_gtf(path: str):
    """Each gene of an Ensembl GTF with all its transcripts, one gene at a time.

    Ensembl writes all of a gene's lines straight after its gene line; this
    relies on that, and stops with an error for a file ordered otherwise.
    """
    gene, txs, seen = None, {}, set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 8)
            if len(f) < 9 or f[2] not in KINDS:
                continue
            kind = f[2]
            start0, end = int(f[3]) - 1, int(f[4])
            if kind == "gene":
                if gene is not None:
                    yield gene, list(txs.values())
                a = _attrs(f[8])
                gene = {
                    "id": a["gene_id"],
                    "name": a.get("gene_name") or a["gene_id"],
                    "contig": browser_name(f[0]),
                    "start": start0,
                    "end": end,
                    "strand": 1 if f[6] == "+" else -1,
                    "biotype": a.get("gene_biotype", ""),
                }
                if gene["id"] in seen:
                    raise ValueError(f"{path}: {gene['id']} has two gene lines")
                seen.add(gene["id"])
                txs = {}
                continue
            gid = _attr(f[8], "gene_id")
            if gene is None or gid != gene["id"]:
                raise ValueError(f"{path}: a line of {gid} is not with the other lines of its gene")
            if kind == "transcript":
                a = _attrs(f[8])
                tags = a.get("tags", set())
                txs[a["transcript_id"]] = {
                    "id": a["transcript_id"],
                    "name": a.get("transcript_name", a["transcript_id"]),
                    "biotype": a.get("transcript_biotype", ""),
                    "flags": (CANONICAL if "Ensembl_canonical" in tags else 0)
                    | (MANE_SELECT if "MANE_Select" in tags else 0)
                    | (MANE_CLINICAL if "MANE_Plus_Clinical" in tags else 0)
                    | (BASIC if "basic" in tags else 0),
                    "start": start0,
                    "end": end,
                    "strand": 1 if f[6] == "+" else -1,
                    "exons": [],
                    "cds": None,  # [lo, hi) of CDS + stop codon
                    "frame": None,
                    "protein": None,  # Ensembl protein ID (ENSP)
                    "sec": [],
                    "incomplete": sorted(t for t in tags if t in INCOMPLETE),
                }
                continue
            t = txs.get(_attr(f[8], "transcript_id"))
            if t is None:
                continue
            if kind == "exon":
                t["exons"].append((start0, end))
            elif kind == "Selenocysteine":
                # a TGA read as Sec: the genomic position of its first base, in transcript orientation
                t["sec"].append(start0 if t["strand"] == 1 else end - 1)
            else:
                lo, hi = t["cds"] or (start0, end)
                t["cds"] = (min(lo, start0), max(hi, end))
                if kind == "CDS":
                    if t["protein"] is None:
                        t["protein"] = _attr(f[8], "protein_id")
                    # frame of the 5'-most CDS piece (in transcript orientation)
                    five_prime = start0 if t["strand"] == 1 else -end
                    if t["frame"] is None or five_prime < t["frame"][0]:
                        t["frame"] = (five_prime, int(f[7]) if f[7] != "." else 0)
    if gene is not None:
        yield gene, list(txs.values())


def _order(t: dict):
    """Canonical first, then MANE Plus Clinical, then coding transcripts, then by name."""
    return (not t["flags"] & CANONICAL, not t["flags"] & MANE_CLINICAL, t["cds"] is None, t["name"])


def canonical_row(t: dict) -> list:
    """The canonical transcript as stored in a gene row (absolute coordinates)."""
    flat = []
    for s, e in sorted(t["exons"]):
        flat += [s - t["start"], e - t["start"]]
    cds = [t["cds"][0], t["cds"][1], t["frame"][1] if t["frame"] else 0, t["protein"] or ""] if t["cds"] else None
    row = [t["id"], t["name"], t["biotype"], int(bool(t["flags"] & MANE_SELECT)), t["start"], t["end"], flat, cds, ",".join(t["incomplete"])]
    if t["sec"]:
        row.append(sorted(t["sec"]))
    return row


def transcript_row(t: dict, g0: int) -> list:
    """A transcript in a transcripts file: coordinates relative to the gene's start."""
    flat = []
    for s, e in sorted(t["exons"]):
        flat += [s - g0, e - g0]
    cds = [t["cds"][0] - g0, t["cds"][1] - g0, t["frame"][1] if t["frame"] else 0, t["protein"] or ""] if t["cds"] else None
    row = [t["id"], t["name"], t["biotype"], t["flags"], t["start"] - g0, t["end"] - g0, flat, cds, ",".join(t["incomplete"])]
    if t["sec"]:
        row.append(sorted(p - g0 for p in t["sec"]))
    return row


def read_hgnc(path: str):
    """HGNC approved name and cross-references by Ensembl gene id, and alias -> symbol."""
    by_ensembl: dict[str, dict] = {}
    aliases: dict[str, set] = defaultdict(set)
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") != "Approved":
                continue
            sym = row["symbol"]
            ens = row.get("ensembl_gene_id") or ""
            if ens:
                by_ensembl[ens] = {
                    "symbol": sym,
                    "name": row.get("name", ""),
                    "hgnc": row["hgnc_id"].replace("HGNC:", ""),
                    "entrez": row.get("entrez_id", ""),
                    "uniprot": (row.get("uniprot_ids") or "").split("|")[0],
                    "omim": (row.get("omim_id") or "").split("|")[0],
                    "location": row.get("location", ""),
                }
            for col in ("alias_symbol", "prev_symbol"):
                for al in (row.get(col) or "").split("|"):
                    al = al.strip()
                    if al and al != sym:
                        aliases[al].add(sym)
    return by_ensembl, aliases


def build_genes(gtf: str, hgnc: str | None, out_dir: str, contigs: list[str]) -> dict:
    log("reading the GTF ...")
    by_ens, aliases = read_hgnc(hgnc) if hgnc else ({}, {})
    known = set(contigs)
    per_contig: dict[str, list] = defaultdict(list)
    bins: dict[str, io.BytesIO] = {}
    index = []
    n_genes = n_tx = 0
    for g, txs in iter_gtf(gtf):
        if g["contig"] not in known:
            continue
        h = by_ens.get(g["id"], {})
        txs = sorted((t for t in txs if t["exons"]), key=_order)
        n_genes += 1
        n_tx += len(txs)
        canonical = next((t for t in txs if t["flags"] & CANONICAL), None)
        where = None
        if txs:
            # all the gene's transcripts, as one gzip member appended to the contig's file
            blob = gzip.compress(json.dumps([transcript_row(t, g["start"]) for t in txs], separators=(",", ":")).encode(), 9, mtime=0)
            buf = bins.setdefault(g["contig"], io.BytesIO())
            where = [buf.tell(), len(blob), len(txs)]
            buf.write(blob)
        xref = [h.get("hgnc", ""), h.get("entrez", ""), h.get("uniprot", ""), h.get("omim", "")] if h else None
        tx = canonical_row(canonical) if canonical else None
        per_contig[g["contig"]].append([g["id"], g["name"], g["start"], g["end"], g["strand"], g["biotype"], h.get("name", ""), tx, xref, where])
        index.append([g["name"], g["id"], g["contig"], g["start"], g["end"], g["strand"], g["biotype"]])
    log(f"  {n_genes:,} genes, {n_tx:,} transcripts")
    sizes = {}
    files = {}
    for c, rows in per_contig.items():
        rows.sort(key=lambda r: (r[2], r[3]))
        f = f"genes/{c}.json.gz"
        sizes[f] = write_json(os.path.join(out_dir, f), {"contig": c, "genes": rows}, gz=True)
        files[c] = {"file": f, "count": len(rows), "transcripts": None}
        if c in bins:
            tf = f"transcripts/{c}.bin"
            sizes[tf] = write_bytes(os.path.join(out_dir, tf), bins[c].getvalue())
            files[c]["transcripts"] = tf
    index.sort(key=lambda r: (r[2], r[3], r[4]))
    sizes["genes/index.json.gz"] = write_json(os.path.join(out_dir, "genes/index.json.gz"), index, gz=True)
    names = {row[0] for row in index}
    alias_rows = sorted([al, sorted(s)] for al, s in aliases.items() if al not in names)
    sizes["genes/aliases.json.gz"] = write_json(os.path.join(out_dir, "genes/aliases.json.gz"), alias_rows, gz=True)
    tx_bytes = sum(v for k, v in sizes.items() if k.startswith("transcripts/"))
    log(f"  gene files: {(sum(sizes.values()) - tx_bytes) / 1e6:.1f} MB for {len(index):,} genes, {len(alias_rows):,} aliases; transcripts: {tx_bytes / 1e6:.1f} MB")
    return {"index": "genes/index.json.gz", "aliases": "genes/aliases.json.gz", "per_contig": files, "count": len(index), "transcripts": n_tx, "bytes": sum(sizes.values())}


def build_bands(assembly_json: str, out_dir: str) -> dict:
    """Cytobands from the Ensembl REST ``info/assembly`` response."""
    import json

    with open(assembly_json) as fh:
        d = json.load(fh)
    bands: dict[str, list] = {}
    for r in d.get("top_level_region", []):
        rows = sorted(([b["start"] - 1, b["end"], b["id"], b["stain"]] for b in r.get("bands", [])), key=lambda x: x[0])
        if rows:
            bands[browser_name(r["name"])] = rows
    n = write_json(os.path.join(out_dir, "cytobands.json.gz"), bands, gz=True)
    log(f"  cytobands: {sum(len(v) for v in bands.values())} bands on {len(bands)} chromosomes")
    return {"file": "cytobands.json.gz", "assembly": d.get("assembly_name"), "bytes": n}
