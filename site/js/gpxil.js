// The browser's genome: each contig's bases, read from the gpxil files of data/grch38 with gpxiltools' JavaScript
// reader (vendor/gpxil-read.js, a copy; vendor/gpxil-read.source.json says from which commit).
//
// A gpxil plane is a raster of W bases per row stored in strips of whole rows (spec section 8), so any interval of
// a chromosome maps to a run of strips: the browser reads only those bytes, the way cloud-optimised GeoTIFFs work.
// The manifest gives each contig's file and the offset of its page's directory, so the reader goes straight to it.

import "./vendor/gpxil-read.js";
import { fetchRange } from "./util.js";

export const GpxilRead = globalThis.GpxilRead;

/** Decoded strips shared by all files: 160 MB of uncompressed bases at most. */
GpxilRead.setCacheBytes(160 * 1024 * 1024);

// what each read is, as the loading indicator names it
const LABEL = { header: "gpxil header", directory: "gpxil page", tags: "gpxil tags", strips: "DNA" };
const fetchLabelled = (url, start, end, what) => fetchRange(url, start, end, LABEL[what] || "gpxil");

/** Opens each contig's file once and reads bases. */
export class Genome {
  constructor(manifest, base) {
    this.manifest = manifest;
    this.base = base; // URL of the data directory
    this.files = new Map();
    this.contigs = new Map(manifest.sequence.contigs.map((c) => [c.name, c]));
    this.aliases = new Map();
    for (const c of manifest.sequence.contigs) {
      this.aliases.set(c.name.toLowerCase(), c.name);
      for (const a of c.aliases || []) this.aliases.set(a.toLowerCase(), c.name);
    }
  }
  resolve(name) {
    if (!name) return null;
    return this.aliases.get(String(name).toLowerCase()) || null;
  }
  contig(name) {
    return this.contigs.get(this.resolve(name));
  }
  /** the gpxil file holding a contig (scaffolds share one) */
  file(contig) {
    const url = new URL(contig.file, this.base).href;
    if (!this.files.has(url)) this.files.set(url, new GpxilRead.TiffFile(GpxilRead.fromUrl(url, { fetchRange: fetchLabelled })));
    return this.files.get(url);
  }
  /** ref codes of contig positions [start, end) (0-based, end exclusive). */
  async bases(name, start, end) {
    const c = this.contig(name);
    if (!c) throw new Error(`unknown contig ${name}`);
    start = Math.max(0, start);
    end = Math.min(c.length, end);
    const f = this.file(c);
    const page = await f.pageAt(c.ifd);
    return f.readPage(page, start, end);
  }
}
