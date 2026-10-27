"""Shared helpers for the data pipeline."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import sys
import urllib.request

# Ensembl names of the primary chromosomes, in karyotype order
PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def browser_name(ensembl_name: str) -> str:
    """'1' -> 'chr1', 'MT' -> 'chrM'; scaffolds keep their Ensembl names."""
    if ensembl_name == "MT":
        return "chrM"
    if ensembl_name in PRIMARY:
        return "chr" + ensembl_name
    return ensembl_name


def aliases_for(ensembl_name: str) -> list[str]:
    name = browser_name(ensembl_name)
    out = [ensembl_name] if ensembl_name != name else []
    if name == "chrM":
        out += ["M", "chrMT"]
    return out


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def fetch(src: str, cache_dir: str) -> str:
    """Return a local path for ``src``: a path is returned as is, a URL is downloaded once into ``cache_dir``."""
    if not src.startswith(("http://", "https://")):
        return src
    os.makedirs(cache_dir, exist_ok=True)
    name = src.rstrip("/").split("/")[-1].split("?")[0] or "download"
    if "?" in src:  # REST queries: name the file by a hash of the URL
        name = hashlib.sha1(src.encode()).hexdigest()[:12] + ".json"
    path = os.path.join(cache_dir, name)
    if not os.path.exists(path):
        log(f"downloading {src}")
        tmp = path + ".part"
        req = urllib.request.Request(src, headers={"User-Agent": "gpxil-browser-pipeline", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh, 1 << 20)
        os.replace(tmp, path)
    return path


def write_json(path: str, obj, *, gz: bool = False) -> int:
    """Write compact JSON (optionally gzip-compressed); returns the number of bytes written."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    data = json.dumps(obj, separators=(",", ":"), ensure_ascii=True).encode()
    if gz:
        data = gzip.compress(data, 9, mtime=0)
    with open(path, "wb") as fh:
        fh.write(data)
    return len(data)


def write_bytes(path: str, data: bytes) -> int:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)
    return len(data)
