// A minimal gpxil reader for the web: TIFF (classic or BigTIFF) parsed over
// HTTP range requests, strips inflated with the browser's DecompressionStream.
//
// A gpxil plane is a raster of W bases per row stored in strips of whole rows
// (spec section 8), so any interval of a chromosome maps to a run of strips:
// the browser reads only those bytes, the way cloud-optimised GeoTIFFs work.

import { LRU, fetchRange, inflate } from "./util.js";

const TAG = {
  ImageWidth: 256,
  ImageLength: 257,
  BitsPerSample: 258,
  Compression: 259,
  ImageDescription: 270,
  StripOffsets: 273,
  SamplesPerPixel: 277,
  RowsPerStrip: 278,
  StripByteCounts: 279,
  PageName: 285,
  Predictor: 317,
  TileWidth: 322,
};
const TYPE_SIZE = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8, 16: 8, 17: 8, 18: 8 };
const IFD_PROBE = 16384; // bytes read at an IFD: enough for the tags and, with tifffile's layout, their arrays

/** Decoded strips shared by all files: 160 MB of uncompressed bases at most. */
const stripCache = new LRU(160 * 1024 * 1024, (a) => a.byteLength);
const pending = new Map();

class Reader {
  constructor(bytes, offset, little) {
    this.v = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    this.base = offset; // file offset of bytes[0]
    this.little = little;
  }
  has(off, n) {
    return off >= this.base && off + n <= this.base + this.v.byteLength;
  }
  u16(off) {
    return this.v.getUint16(off - this.base, this.little);
  }
  u32(off) {
    return this.v.getUint32(off - this.base, this.little);
  }
  u64(off) {
    return Number(this.v.getBigUint64(off - this.base, this.little));
  }
}

export class GpxilFile {
  constructor(url) {
    this.url = url;
    this.headerPromise = null;
    this.pages = new Map(); // IFD offset -> page promise
  }

  async header() {
    if (!this.headerPromise) {
      this.headerPromise = (async () => {
        const b = await fetchRange(this.url, 0, 16, "gpxil header");
        const little = b[0] === 0x49;
        if (!(little || b[0] === 0x4d) || b[0] !== b[1]) throw new Error(`${this.url} is not a TIFF file`);
        const r = new Reader(b, 0, little);
        const magic = r.u16(2);
        if (magic === 42) return { little, big: false, firstIfd: r.u32(4) };
        if (magic === 43) return { little, big: true, firstIfd: r.u64(8) };
        throw new Error(`${this.url}: bad TIFF magic ${magic}`);
      })();
    }
    return this.headerPromise;
  }

  /** The plane stored at an IFD (a page); without an offset, the first page. */
  page(ifdOffset) {
    const key = ifdOffset ?? "first";
    if (!this.pages.has(key)) {
      const p = this.readPage(ifdOffset);
      p.catch(() => this.pages.delete(key));
      this.pages.set(key, p);
    }
    return this.pages.get(key);
  }

  async readPage(ifdOffset) {
    const h = await this.header();
    const at = ifdOffset ?? h.firstIfd;
    let bytes = await fetchRange(this.url, at, at + IFD_PROBE, "gpxil page");
    let r = new Reader(bytes, at, h.little);
    const n = h.big ? r.u64(at) : r.u16(at);
    const entrySize = h.big ? 20 : 12;
    const first = at + (h.big ? 8 : 2);
    if (!r.has(first, n * entrySize)) {
      bytes = await fetchRange(this.url, at, first + n * entrySize + 8, "gpxil page");
      r = new Reader(bytes, at, h.little);
    }
    const tags = new Map();
    for (let i = 0; i < n; i++) {
      const e = first + i * entrySize;
      const tag = r.u16(e);
      const type = r.u16(e + 2);
      const count = h.big ? r.u64(e + 4) : r.u32(e + 4);
      const size = (TYPE_SIZE[type] || 1) * count;
      const inline = size <= (h.big ? 8 : 4);
      const valueAt = inline ? e + (h.big ? 12 : 8) : h.big ? r.u64(e + 12) : r.u32(e + 8);
      tags.set(tag, { type, count, valueAt, size });
    }
    const values = async (tag) => {
      const t = tags.get(tag);
      if (!t) return null;
      let rr = r;
      if (!rr.has(t.valueAt, t.size)) {
        const b = await fetchRange(this.url, t.valueAt, t.valueAt + t.size, "gpxil tags");
        rr = new Reader(b, t.valueAt, h.little);
      }
      const out = new Array(t.count);
      const step = TYPE_SIZE[t.type] || 1;
      for (let i = 0; i < t.count; i++) {
        const o = t.valueAt + i * step;
        out[i] = t.type === 3 ? rr.u16(o) : t.type === 4 || t.type === 9 ? rr.u32(o) : t.type === 16 || t.type === 17 ? rr.u64(o) : rr.v.getUint8(o - rr.base);
      }
      return out;
    };
    const one = async (tag, fallback) => {
      const v = await values(tag);
      return v ? v[0] : fallback;
    };
    if (tags.has(TAG.TileWidth)) throw new Error("tiled gpxil pages are not supported by this reader");
    const page = {
      width: await one(TAG.ImageWidth),
      height: await one(TAG.ImageLength),
      bits: await one(TAG.BitsPerSample, 1),
      samples: await one(TAG.SamplesPerPixel, 1),
      compression: await one(TAG.Compression, 1),
      predictor: await one(TAG.Predictor, 1),
      rowsPerStrip: await one(TAG.RowsPerStrip, 2 ** 32 - 1),
      stripOffsets: await values(TAG.StripOffsets),
      stripByteCounts: await values(TAG.StripByteCounts),
      ifd: at,
      describe: async () => {
        const v = await values(TAG.ImageDescription);
        return v ? new TextDecoder().decode(Uint8Array.from(v.slice(0, -1))) : null;
      },
    };
    if (page.bits !== 8 || page.samples !== 1) throw new Error(`unsupported page: ${page.bits} bits x ${page.samples} samples`);
    if (![1, 8, 32946].includes(page.compression)) throw new Error(`unsupported compression ${page.compression}`);
    if (page.predictor !== 1) throw new Error(`unsupported predictor ${page.predictor}`);
    page.rowsPerStrip = Math.min(page.rowsPerStrip, page.height);
    return page;
  }

  /** Decoded strips [s0, s1] of a page, fetched with as few range requests as possible. */
  async strips(page, s0, s1) {
    const out = new Array(s1 - s0 + 1);
    const missing = [];
    for (let s = s0; s <= s1; s++) {
      const key = `${this.url}#${page.ifd}#${s}`;
      const hit = stripCache.get(key);
      if (hit) out[s - s0] = hit;
      else if (pending.has(key)) out[s - s0] = pending.get(key);
      else missing.push(s);
    }
    // group missing strips that are adjacent in the file into one request (<= 8 MB each)
    const groups = [];
    for (const s of missing) {
      const g = groups[groups.length - 1];
      const off = page.stripOffsets[s];
      if (g && off === g.end && g.end - g.start < 8 << 20) {
        g.strips.push(s);
        g.end = off + page.stripByteCounts[s];
      } else groups.push({ start: off, end: off + page.stripByteCounts[s], strips: [s] });
    }
    for (const g of groups) {
      const fetched = fetchRange(this.url, g.start, g.end, "DNA");
      for (const s of g.strips) {
        const key = `${this.url}#${page.ifd}#${s}`;
        const p = fetched
          .then(async (bytes) => {
            const o = page.stripOffsets[s] - g.start;
            const raw = bytes.subarray(o, o + page.stripByteCounts[s]);
            const data = page.compression === 1 ? raw.slice() : await inflate(raw);
            stripCache.set(key, data);
            return data;
          })
          .finally(() => pending.delete(key));
        pending.set(key, p);
        out[s - s0] = p;
      }
    }
    return Promise.all(out);
  }

  /** Codes for 1-D positions [start, end) of a plane stored at width W. */
  async read(page, start, end) {
    const W = page.width;
    end = Math.min(end, W * page.height);
    if (end <= start) return new Uint8Array(0);
    const r0 = Math.floor(start / W);
    const r1 = Math.floor((end - 1) / W);
    const rps = page.rowsPerStrip;
    const s0 = Math.floor(r0 / rps);
    const s1 = Math.floor(r1 / rps);
    const strips = await this.strips(page, s0, s1);
    const out = new Uint8Array(end - start);
    let pos = 0;
    for (let s = s0; s <= s1; s++) {
      const data = strips[s - s0];
      const stripStart = s * rps * W; // 1-D position of the strip's first pixel
      const a = Math.max(start, stripStart) - stripStart;
      const b = Math.min(end, stripStart + data.length) - stripStart;
      if (b > a) {
        out.set(data.subarray(a, b), pos);
        pos += b - a;
      }
    }
    return out;
  }
}

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
  file(contig) {
    const url = new URL(contig.file, this.base).href;
    if (!this.files.has(url)) this.files.set(url, new GpxilFile(url));
    return this.files.get(url);
  }
  /** ref codes of contig positions [start, end) (0-based, end exclusive). */
  async bases(name, start, end) {
    const c = this.contig(name);
    if (!c) throw new Error(`unknown contig ${name}`);
    start = Math.max(0, start);
    end = Math.min(c.length, end);
    const f = this.file(c);
    const page = await f.page(c.ifd);
    return f.read(page, start, end);
  }
}
