"""Gene annotation for the browser: Ensembl GTF + HGNC names, aliases and cross-references.

Writes, gzip-compressed:

- ``genes/index.json.gz``: every gene, for search: [symbol, gene id, contig,
  start (0-based), end, strand, biotype].
- ``genes/aliases.json.gz``: HGNC alias and previous symbols -> approved symbol.
- ``genes/<contig>.json.gz``: the genes of one contig with their canonical
  transcript (MANE Select where there is one, otherwise Ensembl canonical):
  exons, CDS, selenocysteine codons, and the HGNC name and cross-references.
  Other transcripts are fetched from the Ensembl REST API when a user asks
  for them.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
from collections import defaultdict

from common import browser_name, log, write_json

_ATTR = re.compile(r'(\w+) "([^"]*)"')


def _attrs(text: str) -> dict:
    out: dict = {}
    for k, v in _ATTR.findall(text):
        if k == "tag":
            out.setdefault("tags", set()).add(v)
        else:
            out[k] = v
    return out


def read_gtf(path: str):
    """Genes, transcripts and the exons/CDS of canonical transcripts from an Ensembl GTF."""
    genes: dict[str, dict] = {}
    txs: dict[str, dict] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 8)
            kind = f[2]
            if kind in ("exon", "CDS", "stop_codon", "Selenocysteine") and "Ensembl_canonical" not in f[8]:
                continue  # exons and CDS only for canonical transcripts
            if kind not in ("gene", "transcript", "exon", "CDS", "stop_codon", "Selenocysteine"):
                continue
            a = _attrs(f[8])
            start0, end = int(f[3]) - 1, int(f[4])
            if kind == "gene":
                genes[a["gene_id"]] = {
                    "id": a["gene_id"],
                    "name": a.get("gene_name") or a["gene_id"],
                    "contig": browser_name(f[0]),
                    "start": start0,
                    "end": end,
                    "strand": 1 if f[6] == "+" else -1,
                    "biotype": a.get("gene_biotype", ""),
                }
            elif kind == "transcript":
                tags = a.get("tags", set())
                if "Ensembl_canonical" not in tags:
                    continue
                txs[a["transcript_id"]] = {
                    "id": a["transcript_id"],
                    "gene": a["gene_id"],
                    "name": a.get("transcript_name", a["transcript_id"]),
                    "biotype": a.get("transcript_biotype", ""),
                    "mane": int("MANE_Select" in tags),
                    "start": start0,
                    "end": end,
                    "strand": 1 if f[6] == "+" else -1,
                    "exons": [],
                    "cds": [],
                    "sec": [],
                    "frame": None,
                    "incomplete": sorted(t for t in tags if t in ("cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF")),
                }
            else:
                t = txs.get(a["transcript_id"])
                if t is None:
                    continue
                if kind == "exon":
                    t["exons"].append((start0, end))
                elif kind == "Selenocysteine":
                    # a TGA read as Sec: the genomic position of its first base, in transcript orientation
                    t["sec"].append(start0 if t["strand"] == 1 else end - 1)
                else:
                    t["cds"].append((start0, end))
                    if kind == "CDS":
                        # frame of the 5'-most CDS piece (in transcript orientation)
                        five_prime = start0 if t["strand"] == 1 else -end
                        if t["frame"] is None or five_prime < t["frame"][0]:
                            t["frame"] = (five_prime, int(f[7]) if f[7] != "." else 0)
    return genes, txs


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
    genes, txs = read_gtf(gtf)
    log(f"  {len(genes):,} genes, {len(txs):,} canonical transcripts")
    by_ens, aliases = read_hgnc(hgnc) if hgnc else ({}, {})
    canonical = {t["gene"]: t for t in txs.values()}
    known = set(contigs)
    per_contig: dict[str, list] = defaultdict(list)
    index = []
    for g in sorted(genes.values(), key=lambda g: (g["contig"], g["start"], g["end"])):
        if g["contig"] not in known:
            continue
        h = by_ens.get(g["id"], {})
        name = g["name"]
        t = canonical.get(g["id"])
        tx = None
        if t is not None and t["exons"]:
            ex = sorted(t["exons"])
            flat = []
            for s, e in ex:
                flat += [s - t["start"], e - t["start"]]
            cds = None
            if t["cds"]:
                cds = [min(s for s, _ in t["cds"]), max(e for _, e in t["cds"]), t["frame"][1] if t["frame"] else 0]
            tx = [t["id"], t["name"], t["biotype"], t["mane"], t["start"], t["end"], flat, cds, ",".join(t["incomplete"])]
            if t["sec"]:
                tx.append(sorted(t["sec"]))
        xref = [h.get("hgnc", ""), h.get("entrez", ""), h.get("uniprot", ""), h.get("omim", "")] if h else None
        per_contig[g["contig"]].append([g["id"], name, g["start"], g["end"], g["strand"], g["biotype"], h.get("name", ""), tx, xref])
        index.append([name, g["id"], g["contig"], g["start"], g["end"], g["strand"], g["biotype"]])
    sizes = {}
    files = {}
    for c, rows in per_contig.items():
        f = f"genes/{c}.json.gz"
        sizes[f] = write_json(os.path.join(out_dir, f), {"contig": c, "genes": rows}, gz=True)
        files[c] = {"file": f, "count": len(rows)}
    sizes["genes/index.json.gz"] = write_json(os.path.join(out_dir, "genes/index.json.gz"), index, gz=True)
    names = {row[0] for row in index}
    alias_rows = sorted([al, sorted(s)] for al, s in aliases.items() if al not in names)
    sizes["genes/aliases.json.gz"] = write_json(os.path.join(out_dir, "genes/aliases.json.gz"), alias_rows, gz=True)
    log(f"  gene files: {sum(sizes.values()) / 1e6:.1f} MB for {len(index):,} genes, {len(alias_rows):,} aliases")
    return {"index": "genes/index.json.gz", "aliases": "genes/aliases.json.gz", "per_contig": files, "count": len(index), "bytes": sum(sizes.values())}


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
