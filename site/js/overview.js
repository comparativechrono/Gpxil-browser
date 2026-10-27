// Summary bins for zoomed-out views. For bins of 128, 512, ... bases each
// contig has two bytes per bin: GC fraction of A/C/G/T (0-255), and
// N fraction << 4 | soft-masked fraction (0-15 each). See pipeline/sequence.py.

import { fetchBytes, gunzip } from "./util.js";

export class Overview {
  constructor(manifest, base) {
    this.bins = manifest.sequence.overview_bins;
    this.base = base;
    this.files = new Map(); // url -> Promise<Uint8Array>
  }

  file(url) {
    if (!this.files.has(url)) {
      // gunzip here, unless the server already did (gzip magic: 1f 8b)
      const p = fetchBytes(url, "overview").then((b) => (b[0] === 0x1f && b[1] === 0x8b ? gunzip(b) : b));
      p.catch(() => this.files.delete(url));
      this.files.set(url, p);
    }
    return this.files.get(url);
  }

  /** The level whose bins are the largest not exceeding `bpp` bases (at least the finest). */
  levelFor(bpp) {
    let li = 0;
    for (let i = 0; i < this.bins.length; i++) if (this.bins[i] <= bpp) li = i;
    return li;
  }

  /** Uint8Array of 2 bytes per bin for one contig and level, or a promise of it. */
  async level(contig, li) {
    const all = await this.file(new URL(contig.overview.file, this.base).href);
    const n = Math.ceil(contig.length / this.bins[li]);
    const off = contig.overview.offsets[li];
    return all.subarray(off, off + 2 * n);
  }

  isLoaded(contig) {
    return this.files.has(new URL(contig.overview.file, this.base).href);
  }
}

/** Decodes one bin: [gc 0..1, N 0..1, masked 0..1]. */
export function decodeBin(bytes, i, out = [0, 0, 0]) {
  out[0] = bytes[2 * i] / 255;
  const b = bytes[2 * i + 1];
  out[1] = (b >> 4) / 15;
  out[2] = (b & 15) / 15;
  return out;
}
