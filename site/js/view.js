// The pixel viewer: one chromosome as rows of pixels, one base per pixel.
//
// Layout. In "auto" wrap the rows are as wide as the screen, so the view is
// always full: zoomed out, a pixel is a run of `bpp` bases (a summary bin);
// zoomed in, a base is an s x s block, and rows hold fewer bases. Zooming
// reflows the rows, like text, around the base under the cursor. In "fixed"
// wrap the row width is chosen by the user (e.g. a repeat period) and the
// view scrolls sideways when rows are wider than the screen.
//
// Data. Zoomed in, pixels come from the gpxil file itself (strips read by
// range requests); zoomed out, from the precomputed summary bins.

import { buildLut, binColour, boost, IUPAC, inkOn } from "./palettes.js";
import { decodeBin } from "./overview.js";
import { REG_COL } from "./live.js";
import { Emitter, clamp, fmt, fmtBp } from "./util.js";

// zoom ladder, in screen pixels per base: whole bases per pixel (262,144 ... 3, 2,
// in steps of about 1.4x), then whole pixels per base (1 ... 64)
export const SCALES = [];
for (let k = 18; k >= 1; k--) {
  SCALES.push(1 / 2 ** k);
  if (k >= 2) SCALES.push(1 / (3 * 2 ** (k - 2)));
}
SCALES.push(1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64);
const IS_ACGT = new Uint8Array(256);
const IS_GC = new Uint8Array(256);
for (const c of [1, 2, 4, 8]) IS_ACGT[c] = IS_ACGT[c | 16] = 1;
for (const c of [2, 4]) IS_GC[c] = IS_GC[c | 16] = 1;
const BASE_LIMIT = 24e6; // read bases (not bins) when the view shows at most this many
export const MAX_WRAP = 16384; // widest fixed wrap: keeps a screenful of rows under BASE_LIMIT
const LETTER_SCALE = 12; // draw base letters from this many pixels per base
const GUTTER_STEP = 22; // minimum pixels between ruler labels
// which regulatory element wins a row in the gutter's live track
const LIVE_RANK = { severe: 9, promoter: 5, CTCF_binding_site: 4, enhancer: 3, open_chromatin_region: 2, TF_binding_site: 1 };

export class PixelView extends Emitter {
  constructor({ canvas, gutter, genome, overview, genes }) {
    super();
    this.canvas = canvas;
    this.gutter = gutter;
    this.genome = genome;
    this.overview = overview;
    this.genes = genes;
    this.ctx = canvas.getContext("2d", { alpha: false });
    this.frame = document.createElement("canvas");
    this.overlay = document.createElement("canvas");
    this.state = { contig: null, zi: 0, scrollRow: 0, scrollX: 0, wrap: "auto", width: 171 };
    this.opts = { palette: "classic", repeats: "fade", letters: true, grid: true, genesTrack: true, exons: true, cpg: false, contrast: true };
    this.lut = buildLut(this.opts.palette);
    this.data = null; // {contig, start, end, codes}
    this.request = null;
    this.highlight = null; // {gene, tx, mode}
    this.selection = null; // [start, end)
    this.hover = null;
    this.pulse = 0;
    this.dirtyFrame = true;
    this.raf = 0;
    this.css = { w: 1, h: 1, dpr: 1 };
    this.genesInView = [];
    this.live = null; // live Ensembl tracks, set by the app
    this.rowGC = null;
    this.resizeObserver = new ResizeObserver(() => this.resize());
    this.resizeObserver.observe(canvas.parentElement);
    this.bindEvents();
    this.resize();
  }

  // ------------------------------------------------------------------ setup
  setOptions(o) {
    const before = { ...this.opts };
    Object.assign(this.opts, o);
    if (before.palette !== this.opts.palette || before.repeats !== this.opts.repeats) {
      this.lutVersion = (this.lutVersion || 0) + 1;
      this.lut = buildLut(this.opts.palette);
      if (this.opts.repeats !== "fade") this.adjustRepeats();
      this.lut.mean = null;
      this.lut.mean = binColour(this.lut, 0.41, 0, this.opts.repeats === "fade" ? 0.5 : 0);
    }
    this.dirtyFrame = true;
    this.schedule();
  }

  adjustRepeats() {
    // "off": soft-masked bases look like the others; "hide": they fade into the background
    const plain = buildLut(this.opts.palette);
    const p = plain.palette;
    for (let code = 0x11; code < 0x20; code++) {
      const src = this.opts.repeats === "off" ? code & 15 : code;
      let [r, g, b] = [plain.rgb[3 * src], plain.rgb[3 * src + 1], plain.rgb[3 * src + 2]];
      if (this.opts.repeats === "hide") [r, g, b] = [r, g, b].map((v, i) => Math.round(v * 0.12 + [8, 8, 18][i] * 0.88));
      this.lut.rgb.set([r, g, b], 3 * code);
      this.lut.packed[code] = this.lut.little ? (255 << 24) | (b << 16) | (g << 8) | r : (r << 24) | (g << 16) | (b << 8) | 255;
    }
    void p;
  }

  resize() {
    const box = this.canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    const gw = this.gutter ? this.gutter.getBoundingClientRect().width : 0;
    const w = Math.max(50, Math.floor(box.width));
    const h = Math.max(50, Math.floor(box.height));
    if (w === this.css.w && h === this.css.h && dpr === this.css.dpr) return;
    const anchor = this.state.contig ? this.centreBase() : null;
    this.css = { w, h, dpr };
    for (const c of [this.canvas]) {
      c.width = Math.round(w * dpr);
      c.height = Math.round(h * dpr);
      c.style.width = `${w}px`;
      c.style.height = `${h}px`;
    }
    if (this.gutter) {
      this.gutter.width = Math.round(gw * dpr);
      this.gutter.height = Math.round(h * dpr);
      this.gutter.style.height = `${h}px`;
    }
    this.overlay.width = this.canvas.width;
    this.overlay.height = this.canvas.height;
    if (anchor !== null) this.centreOn(anchor);
    this.dirtyFrame = true;
    this.schedule();
  }

  get contig() {
    return this.state.contig ? this.genome.contig(this.state.contig) : null;
  }

  // ------------------------------------------------------------------ layout
  layout(zi = this.state.zi) {
    const s = SCALES[zi];
    const { w: cw } = this.css;
    if (this.state.wrap === "fixed" && s >= 1) {
      const w = this.state.width;
      return { s, bpp: 1, w, rowH: s, cellW: s, cols: w, rowPx: w * s };
    }
    if (s >= 1) {
      const w = Math.max(1, Math.floor(cw / s));
      return { s, bpp: 1, w, rowH: s, cellW: s, cols: w, rowPx: w * s };
    }
    const bpp = Math.round(1 / s);
    return { s, bpp, w: cw * bpp, rowH: 1, cellW: 1, cols: cw, rowPx: cw };
  }

  totalRows(L = this.layout()) {
    const c = this.contig;
    return c ? Math.ceil(c.length / L.w) : 0;
  }

  visibleRows(L = this.layout()) {
    return this.css.h / L.rowH;
  }

  /**
   * Zoom index that shows [start, end) whole, as large as possible. By default
   * with a margin; `exact` restores a window this view described (a link), so
   * reloading or sharing a view keeps its zoom.
   */
  fitZoom(start, end, { padding = 0.85, exact = false } = {}) {
    const len = Math.max(1, end - start);
    const fixed = this.state.wrap === "fixed";
    const first = fixed ? SCALES.indexOf(1) : 0; // fixed rows are whole bases: never zoomed out past 1 px per base
    let best = first;
    for (let zi = first; zi < SCALES.length; zi++) {
      const L = this.layout(zi);
      const rows = Math.ceil(len / L.w);
      // a described window is whole rows, rounded out at both ends
      const fits = exact ? rows <= Math.ceil(this.visibleRows(L)) + 1 : (rows + 1) * L.rowH <= this.css.h * padding && (!fixed || L.rowPx <= this.css.w);
      if (fits) best = zi;
    }
    return best;
  }

  maxScroll(L = this.layout()) {
    return Math.max(0, this.totalRows(L) - this.visibleRows(L) * 0.5);
  }

  centreBase() {
    const L = this.layout();
    const row = this.state.scrollRow + this.visibleRows(L) / 2;
    return clamp(Math.floor(row) * L.w + Math.floor(L.w / 2), 0, this.contig.length - 1);
  }

  centreOn(pos, L = this.layout()) {
    const row = pos / L.w;
    this.state.scrollRow = clamp(row - this.visibleRows(L) / 2, 0, this.maxScroll(L));
    if (this.state.wrap === "fixed") {
      const x = (pos % L.w) * L.cellW;
      this.state.scrollX = clamp(x - this.css.w / 2, 0, Math.max(0, L.rowPx - this.css.w));
    }
  }

  /** Visible 1-D window [start, end). */
  window(L = this.layout()) {
    const c = this.contig;
    if (!c) return [0, 0];
    const r0 = Math.max(0, Math.floor(this.state.scrollRow));
    const r1 = Math.ceil(this.state.scrollRow + this.visibleRows(L));
    return [clamp(r0 * L.w, 0, c.length), clamp(r1 * L.w, 0, c.length)];
  }

  // ------------------------------------------------------------------ navigation
  show(contigName, start = null, end = null, { zoom = null, exact = false } = {}) {
    const c = this.genome.contig(contigName);
    if (!c) throw new Error(`unknown sequence ${contigName}`);
    const changed = this.state.contig !== c.name;
    this.state.contig = c.name;
    if (changed) {
      this.data = null;
      this.selection = null;
      this.state.scrollX = 0;
    }
    if (start === null) {
      start = 0;
      end = c.length;
    }
    start = clamp(start, 0, c.length - 1);
    end = clamp(end ?? start + 1, start + 1, c.length);
    this.state.zi = zoom ?? this.fitZoom(start, end, { exact });
    this.centreOn(Math.floor((start + end) / 2));
    this.dirtyFrame = true;
    this.emit("navigate", this.describe());
    this.schedule();
  }

  zoomAt(steps, x = this.css.w / 2, y = this.css.h / 2) {
    if (!this.contig) return;
    const pos = this.baseAt(x, y, true) ?? this.centreBase();
    let zi = clamp(this.state.zi + steps, 0, SCALES.length - 1);
    if (this.state.wrap === "fixed") zi = Math.max(zi, SCALES.indexOf(1));
    // don't zoom out further than needed to see the whole contig
    const fit = this.fitZoom(0, this.contig.length, { padding: 1 });
    zi = Math.max(zi, Math.min(fit, this.state.zi));
    if (zi === this.state.zi) return;
    this.state.zi = zi;
    const L = this.layout();
    this.state.scrollRow = clamp(pos / L.w - y / L.rowH, 0, this.maxScroll(L));
    if (this.state.wrap === "fixed") this.state.scrollX = clamp((pos % L.w) * L.cellW - x, 0, Math.max(0, L.rowPx - this.css.w));
    this.dirtyFrame = true;
    this.emit("navigate", this.describe());
    this.schedule();
  }

  scrollRows(d) {
    const L = this.layout();
    this.state.scrollRow = clamp(this.state.scrollRow + d, 0, this.maxScroll(L));
    this.dirtyFrame = true;
    this.emit("navigate", this.describe());
    this.schedule();
  }

  scrollPx(dx, dy) {
    const L = this.layout();
    if (dy) this.state.scrollRow = clamp(this.state.scrollRow + dy / L.rowH, 0, this.maxScroll(L));
    if (dx && this.state.wrap === "fixed") this.state.scrollX = clamp(this.state.scrollX + dx, 0, Math.max(0, L.rowPx - this.css.w));
    this.dirtyFrame = true;
    this.emit("navigate", this.describe());
    this.schedule();
  }

  setWrap(mode, width = this.state.width) {
    const pos = this.contig ? this.centreBase() : 0;
    this.state.wrap = mode;
    this.state.width = clamp(Math.round(width) || 171, 1, MAX_WRAP);
    if (mode === "fixed") this.state.zi = Math.max(this.state.zi, SCALES.indexOf(1));
    this.state.scrollX = 0;
    if (this.contig) this.centreOn(pos);
    this.dirtyFrame = true;
    this.emit("navigate", this.describe());
    this.schedule();
  }

  describe() {
    const c = this.contig;
    if (!c) return null;
    const L = this.layout();
    const [start, end] = this.window(L);
    return { contig: c.name, start, end, length: c.length, bpp: L.bpp, s: L.s, w: L.w, wrap: this.state.wrap };
  }

  // ------------------------------------------------------------------ geometry
  /** 0-based position under a screen point (CSS px), or null. */
  baseAt(x, y, clampToContig = false) {
    const c = this.contig;
    if (!c) return null;
    const L = this.layout();
    const row = Math.floor(this.state.scrollRow + y / L.rowH);
    const col = Math.floor((x + (this.state.wrap === "fixed" ? this.state.scrollX : 0)) / L.cellW) * L.bpp;
    if (col < 0 || col >= L.w || row < 0) return clampToContig ? clamp(row * L.w, 0, c.length - 1) : null;
    const pos = row * L.w + col;
    if (pos >= c.length) return clampToContig ? c.length - 1 : null;
    return pos;
  }

  /** Screen rectangles [x, y, w, h] covering positions [start, end). */
  rects(start, end, L = this.layout()) {
    const out = [];
    if (end <= start) return out;
    const { w, rowH, cellW, bpp } = L;
    const sx = this.state.wrap === "fixed" ? this.state.scrollX : 0;
    const top = this.state.scrollRow;
    const vis = this.visibleRows(L);
    let ra = Math.floor(start / w);
    let rb = Math.floor((end - 1) / w);
    if (rb < top - 1 || ra > top + vis + 1) return out;
    const xA = Math.floor((start - ra * w) / bpp) * cellW - sx;
    const xB = Math.ceil((end - rb * w) / bpp) * cellW - sx;
    const fullW = Math.ceil(w / bpp) * cellW;
    const y = (r) => (r - top) * rowH;
    if (ra === rb) return [[xA, y(ra), Math.max(xB - xA, 1), rowH]];
    out.push([xA, y(ra), fullW - sx - xA, rowH]);
    const m0 = Math.max(ra + 1, Math.floor(top) - 1);
    const m1 = Math.min(rb - 1, Math.ceil(top + vis) + 1);
    if (m1 >= m0) out.push([-sx, y(m0), fullW, (m1 - m0 + 1) * rowH]);
    out.push([-sx, y(rb), xB + sx, rowH]);
    return out;
  }

  // ------------------------------------------------------------------ data
  wantBases(L, start, end) {
    return L.s >= 1 || end - start <= BASE_LIMIT;
  }

  ensureBases(start, end) {
    const c = this.contig;
    const d = this.data;
    if (d && d.contig === c.name && d.start <= start && d.end >= end) return d;
    const r = this.request;
    if (r && r.contig === c.name && r.start <= start && r.end >= end) return null;
    // read a margin around the window so small scrolls need no new request
    const margin = Math.min(Math.floor((end - start) * 0.5), 4e6);
    const a = clamp(start - margin, 0, c.length);
    const b = clamp(end + margin, 0, c.length);
    const req = { contig: c.name, start: a, end: b };
    this.request = req;
    this.genome
      .bases(c.name, a, b)
      .then((codes) => {
        if (this.request === req) this.request = null;
        if (this.state.contig !== req.contig) return;
        this.data = { ...req, codes };
        this.dirtyFrame = true;
        this.schedule();
        this.emit("data");
      })
      .catch((e) => {
        if (this.request === req) this.request = null;
        this.emit("error", e);
      });
    return null;
  }

  ensureBins(li) {
    const c = this.contig;
    const key = `${c.name}:${li}`;
    if (this.bins && this.bins.key === key) return this.bins;
    if (!this.binsPending || this.binsPending !== key) {
      this.binsPending = key;
      this.overview
        .level(c, li)
        .then((bytes) => {
          this.binsPending = null;
          if (this.state.contig !== c.name) return;
          this.bins = { key, li, size: this.overview.bins[li], bytes };
          this.dirtyFrame = true;
          this.schedule();
        })
        .catch((e) => {
          this.binsPending = null;
          this.emit("error", e);
        });
    }
    return this.bins && this.bins.key.startsWith(`${c.name}:`) ? this.bins : null;
  }

  /** Bases of [start, end) if they are loaded (for info panels and letters). */
  codesFor(start, end) {
    const d = this.data;
    if (!d || d.contig !== this.state.contig || d.start > start || d.end < end) return null;
    return d.codes.subarray(start - d.start, end - d.start);
  }

  // ------------------------------------------------------------------ drawing
  schedule() {
    if (!this.raf) this.raf = requestAnimationFrame(() => this.draw());
  }

  draw() {
    this.raf = 0;
    if (this.paused) return; // covered by the transcript lanes (Expanded mode)
    const c = this.contig;
    const { ctx } = this;
    const { dpr, w: cw, h: ch } = this.css;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    if (!c) {
      ctx.fillStyle = "#07070f";
      ctx.fillRect(0, 0, cw, ch);
      return;
    }
    const L = this.layout();
    if (this.dirtyFrame) {
      this.buildFrame(L);
      this.dirtyFrame = false;
    }
    ctx.fillStyle = "#07070f";
    ctx.fillRect(0, 0, cw, ch);
    ctx.imageSmoothingEnabled = false;
    const f = this.frameInfo;
    if (f) ctx.drawImage(this.frame, 0, 0, f.fw, f.fh, f.dx, f.dy, f.fw * f.sx, f.fh * f.sy);
    this.drawOverlays(L);
    this.drawGutter(L);
    if (this.pulse > 0) {
      this.pulse = Math.max(0, this.pulse - 0.03);
      this.schedule();
    }
  }

  buildFrame(L) {
    const c = this.contig;
    const top = this.state.scrollRow;
    const r0 = Math.max(0, Math.floor(top));
    const r1 = Math.min(this.totalRows(L), Math.ceil(top + this.visibleRows(L)) + 1);
    const start = r0 * L.w;
    const end = Math.min(c.length, r1 * L.w);
    this.frameInfo = null;
    this.frameSource = null;
    if (end <= start) return;
    let source = null;
    if (this.wantBases(L, start, end)) {
      const d = this.ensureBases(start, end);
      if (d) source = { kind: "bases", d };
      else if (L.s < 1) {
        const b = this.ensureBins(this.overview.levelFor(L.bpp));
        if (b) source = { kind: "bins", b };
      } else if (this.data && this.data.contig === c.name) source = { kind: "bases", d: this.data, partial: true };
    } else {
      const b = this.ensureBins(this.overview.levelFor(L.bpp));
      if (b) source = { kind: "bins", b };
    }
    this.frameSource = source;
    if (!source) return;
    const sx = this.state.wrap === "fixed" ? this.state.scrollX : 0;
    const c0 = L.s >= 1 ? Math.max(0, Math.floor(sx / L.cellW)) : 0;
    const c1 = L.s >= 1 ? Math.min(L.w, Math.ceil((sx + this.css.w) / L.cellW) + 1) : L.cols;
    const fw = Math.max(1, c1 - c0);
    const fh = Math.max(1, r1 - r0);
    if (this.frame.width < fw || this.frame.height < fh) {
      this.frame.width = Math.max(this.frame.width, fw);
      this.frame.height = Math.max(this.frame.height, fh);
    }
    const fctx = this.frame.getContext("2d");
    const img = fctx.createImageData(fw, fh);
    const px = new Uint32Array(img.data.buffer);
    const bg = this.lut.little ? 0xff0f0707 : 0x07070fff;
    px.fill(bg);
    let gcs = null;
    if (source.kind === "bases") gcs = this.fillFromBases(px, fw, fh, r0, c0, L, source.d);
    else this.fillFromBins(px, fw, fh, r0, L, source.b);
    this.rowGC = { r0, values: gcs || this.computeRowGC(source, r0, r1, L) };
    fctx.putImageData(img, 0, 0);
    this.frameInfo = {
      fw,
      fh,
      dx: c0 * L.cellW - sx,
      dy: (r0 - top) * L.rowH,
      sx: L.cellW,
      sy: L.rowH,
    };
  }

  fillFromBases(px, fw, fh, r0, c0, L, d) {
    const { packed } = this.lut;
    const n = this.contig.length;
    if (L.bpp === 1) {
      for (let r = 0; r < fh; r++) {
        const rowStart = (r0 + r) * L.w;
        for (let c = 0; c < fw; c++) {
          const pos = rowStart + c0 + c;
          if (c0 + c >= L.w || pos >= n) break;
          if (pos < d.start || pos >= d.end) continue;
          px[r * fw + c] = packed[d.codes[pos - d.start]];
        }
      }
      return null;
    }
    // several bases per pixel: the mean colour of the A/C/G/T bases, stretched away from
    // the genome average, then mixed with the N colour by the fraction of other bases.
    // That means reading every base on screen (millions), so finished rows are kept
    // and scrolling only computes the rows it uncovers.
    const cache = this.rowCache(L);
    const gcs = new Float32Array(fh).fill(NaN);
    for (let r = 0; r < fh; r++) {
      const R = r0 + r;
      const rowStart = R * L.w;
      if (rowStart >= n) break;
      let row = cache.rows.get(R);
      if (!row) {
        row = this.binRow(rowStart, Math.min(n, rowStart + L.w), L, d, fw);
        if (row.complete) cache.rows.set(R, row);
      }
      px.set(row.px.subarray(0, fw), r * fw);
      gcs[r] = row.gc;
    }
    // keep the cache to a few screens around the view
    if (cache.rows.size > 3 * fh + 64) for (const R of cache.rows.keys()) if (R < r0 - fh || R > r0 + 2 * fh) cache.rows.delete(R);
    return gcs;
  }

  rowCache(L) {
    const key = `${this.state.contig}|${L.w}|${L.bpp}|${this.lutVersion || 0}|${this.contrastFor(L.bpp)}`;
    if (!this.rowsCache || this.rowsCache.key !== key) this.rowsCache = { key, rows: new Map() };
    return this.rowsCache;
  }

  /** One zoomed-out row of pixels from bases [a, b), and its GC fraction. */
  binRow(a, b, L, d, fw) {
    const { rgb, little } = this.lut;
    const k = L.bpp;
    const kc = this.contrastFor(k);
    const N = this.lut.palette.rgb.N;
    const px = new Uint32Array(fw);
    const c3 = [0, 0, 0];
    let gc = 0;
    let acgt = 0;
    for (let c = 0; c < fw; c++) {
      const p0 = a + c * k;
      if (p0 >= b) break;
      const s = Math.max(p0, d.start);
      const e = Math.min(p0 + k, b, d.end);
      if (e <= s) continue;
      let R = 0;
      let G = 0;
      let B = 0;
      let nb = 0;
      for (let p = s; p < e; p++) {
        const code = d.codes[p - d.start];
        const low = code & 15;
        if (low === 1 || low === 2 || low === 4 || low === 8) {
          const o = code * 3;
          R += rgb[o];
          G += rgb[o + 1];
          B += rgb[o + 2];
          nb++;
          if (low === 2 || low === 4) gc++;
        }
      }
      acgt += nb;
      if (nb) {
        c3[0] = R / nb;
        c3[1] = G / nb;
        c3[2] = B / nb;
        boost(this.lut, c3, kc);
      }
      const f = nb / (e - s);
      const rr = (c3[0] * f + N[0] * (1 - f)) | 0;
      const gg = (c3[1] * f + N[1] * (1 - f)) | 0;
      const bb = (c3[2] * f + N[2] * (1 - f)) | 0;
      px[c] = little ? (255 << 24) | (bb << 16) | (gg << 8) | rr : (rr << 24) | (gg << 16) | (bb << 8) | 255;
    }
    return { px, gc: acgt ? gc / acgt : NaN, complete: d.start <= a && d.end >= b };
  }

  computeRowGC(source, r0, r1, L) {
    const out = new Float32Array(r1 - r0).fill(NaN);
    const n = this.contig.length;
    if (source.kind === "bases") {
      const d = source.d;
      for (let r = r0; r < r1; r++) {
        const a = Math.max(r * L.w, d.start);
        const b = Math.min((r + 1) * L.w, n, d.end);
        let gc = 0;
        let acgt = 0;
        for (let p = a; p < b; p++) {
          const c = d.codes[p - d.start];
          acgt += IS_ACGT[c];
          gc += IS_GC[c];
        }
        if (acgt) out[r - r0] = gc / acgt;
      }
    } else {
      const b = source.b;
      const nb = b.bytes.length / 2;
      for (let r = r0; r < r1; r++) {
        const i0 = Math.floor((r * L.w) / b.size);
        const i1 = Math.min(nb - 1, Math.floor((Math.min(n, (r + 1) * L.w) - 1) / b.size));
        let gc = 0;
        let wsum = 0;
        for (let i = i0; i <= i1; i++) {
          const wgt = 1 - (b.bytes[2 * i + 1] >> 4) / 15;
          gc += (b.bytes[2 * i] / 255) * wgt;
          wsum += wgt;
        }
        if (wsum > 0.05) out[r - r0] = gc / wsum;
      }
    }
    return out;
  }

  contrastFor(bpp) {
    return this.opts.contrast && bpp > 1 ? 1 + 1.1 * Math.min(1, Math.log2(bpp) / 7) : 1;
  }

  fillFromBins(px, fw, fh, r0, L, b) {
    const n = this.contig.length;
    const size = b.size;
    const nb = b.bytes.length / 2;
    const little = this.lut.little;
    // colours of the bins in view, computed once per frame
    const first = Math.floor((r0 * L.w) / size);
    const last = Math.min(nb - 1, Math.floor(Math.min(n - 1, (r0 + fh) * L.w - 1) / size));
    const cols = new Float32Array(3 * Math.max(0, last - first + 1));
    const tmp = [0, 0, 0];
    const col = [0, 0, 0];
    const kc = this.contrastFor(Math.max(L.bpp, size));
    for (let i = first; i <= last; i++) {
      decodeBin(b.bytes, i, tmp);
      binColour(this.lut, tmp[0], tmp[1], tmp[2], col, kc);
      const o = 3 * (i - first);
      cols[o] = col[0];
      cols[o + 1] = col[1];
      cols[o + 2] = col[2];
    }
    for (let r = 0; r < fh; r++) {
      const rowStart = (r0 + r) * L.w;
      for (let c = 0; c < fw; c++) {
        const p0 = rowStart + c * L.bpp;
        if (p0 >= n) break;
        const p1 = Math.min(n, p0 + L.bpp);
        const i0 = Math.floor(p0 / size);
        const i1 = Math.min(last, Math.floor((p1 - 1) / size));
        let R = 0;
        let G = 0;
        let B = 0;
        for (let i = i0; i <= i1; i++) {
          const o = 3 * (i - first);
          R += cols[o];
          G += cols[o + 1];
          B += cols[o + 2];
        }
        const m = i1 - i0 + 1;
        const rr = (R / m) | 0;
        const gg = (G / m) | 0;
        const bb = (B / m) | 0;
        px[r * fw + c] = little ? (255 << 24) | (bb << 16) | (gg << 8) | rr : (rr << 24) | (gg << 16) | (bb << 8) | 255;
      }
    }
  }

  drawOverlays(L) {
    const { ctx } = this;
    const c = this.contig;
    const [ws, we] = this.window(L);
    // genes in view; their labels are drawn last, above everything else
    this.genesInView = this.genes ? this.genes.overlapping(c.name, ws, we) : [];
    this.labels = [];
    const h = this.highlight;
    if (h && h.gene && h.gene.contig === c.name) this.drawHighlight(L, h);
    else if (this.opts.genesTrack && we - ws <= 40e6) this.drawGeneOutlines(L);
    if (this.opts.cpg && L.s >= 3) this.drawCpG(L, ws, we);
    if (this.live) this.drawLive(L, ws, we);
    if (this.opts.letters && L.s >= LETTER_SCALE) this.drawLetters(L);
    if (this.opts.grid && L.s >= 8) this.drawGrid(L);
    for (const args of this.labels) label(ctx, ...args);
    if (this.selection) {
      ctx.fillStyle = "rgba(32, 227, 255, 0.28)";
      ctx.strokeStyle = "rgba(32, 227, 255, 0.95)";
      ctx.lineWidth = 1;
      for (const [x, y, w, hh] of this.rects(this.selection[0], this.selection[1], L)) {
        ctx.fillRect(x, y, w, hh);
        ctx.strokeRect(x + 0.5, y + 0.5, w - 1, hh - 1);
      }
    }
    if (this.hover !== null && this.hover !== undefined && L.s >= 2) {
      const [x, y, w, hh] = this.rects(this.hover, this.hover + 1, L)[0] || [];
      if (x !== undefined) {
        ctx.strokeStyle = "#fff";
        ctx.lineWidth = 2;
        ctx.strokeRect(x + 1, y + 1, w - 2, hh - 2);
      }
    }
    if (!this.frameSource) {
      ctx.fillStyle = "rgba(7,7,15,0.6)";
      ctx.fillRect(0, 0, this.css.w, this.css.h);
    }
  }

  drawGeneOutlines(L) {
    const { ctx } = this;
    const [ws, we] = this.window(L);
    const span = we - ws;
    const named = (g) => !/^ENSG\d/.test(g.name);
    const rank = (g) => (g.biotype === "protein_coding" ? 0 : 2) + (named(g) ? 0 : 4) + (g.biotype === "lncRNA" ? 1 : 0);
    const shown = this.genesInView.filter((g) => (span > 5e6 ? g.biotype === "protein_coding" && named(g) : span > 1e6 ? named(g) : true));
    const strong = L.rowH >= 3;
    for (const g of shown) {
      const rs = this.rects(g.start, g.end, L);
      const coding = g.biotype === "protein_coding";
      if (strong) {
        // a light line on a dark edge, visible over any pixels
        ctx.lineWidth = 3;
        ctx.strokeStyle = "rgba(0,0,0,0.45)";
        outline(ctx, rs, 1.5);
        ctx.lineWidth = 1;
        ctx.strokeStyle = coding ? "rgba(255,255,255,0.8)" : "rgba(170,170,255,0.6)";
        outline(ctx, rs, 1.5);
      } else {
        // rows one pixel high: a fine line, so the outline doesn't hide whole rows
        ctx.lineWidth = 1;
        ctx.strokeStyle = coding ? "rgba(255,255,255,0.55)" : "rgba(170,170,255,0.35)";
        outline(ctx, rs, 0.5);
      }
      if (this.opts.exons && g.tx && L.s >= 1 / 8 && span <= 5e6) {
        ctx.fillStyle = "rgba(255,255,255,0.18)";
        for (const [a, b] of g.tx.exons) for (const [x, y, w, h] of this.rects(a, b, L)) ctx.fillRect(x, y, w, h);
      }
    }
    // labels: most important first, never on top of each other
    const placed = [];
    let n = 0;
    for (const g of [...shown].sort((a, b) => rank(a) - rank(b) || b.end - b.start - (a.end - a.start))) {
      if (n >= 90) break;
      const rs = this.rects(g.start, g.end, L);
      if (!rs.length) continue;
      const [x0, y0] = rs[0];
      const x = Math.max(2, Math.min(x0, this.css.w - 60));
      const y = Math.max(18, y0 - 1);
      if (y > this.css.h + 16) continue;
      const w = labelWidth(ctx, g.name, g.strand);
      if (placed.some((r) => x < r[0] + r[2] + 4 && x + w > r[0] - 4 && y - 16 < r[1] && y > r[1] - 16)) continue;
      this.labels.push([g.name, x, y, g.strand]);
      placed.push([x, y, w]);
      n++;
    }
  }

  drawHighlight(L, h) {
    const o = this.overlay;
    const octx = o.getContext("2d");
    const { dpr } = this.css;
    octx.setTransform(dpr, 0, 0, dpr, 0, 0);
    octx.globalCompositeOperation = "source-over";
    octx.clearRect(0, 0, this.css.w, this.css.h);
    octx.fillStyle = "rgba(4,4,12,0.78)";
    octx.fillRect(0, 0, this.css.w, this.css.h);
    const g = h.gene;
    const geneRects = this.rects(g.start, g.end, L);
    const tx = h.mode === "gene" ? null : h.tx;
    octx.globalCompositeOperation = "destination-out";
    octx.fillStyle = "#000";
    for (const r of geneRects) octx.fillRect(...r);
    if (tx) {
      octx.globalCompositeOperation = "source-over";
      octx.fillStyle = "rgba(4,4,12,0.66)"; // introns
      for (const r of geneRects) octx.fillRect(...r);
      octx.globalCompositeOperation = "destination-out";
      for (const [a, b] of tx.exons) for (const r of this.rects(a, b, L)) octx.fillRect(...r);
      if (tx.cds) {
        octx.globalCompositeOperation = "source-over";
        octx.fillStyle = "rgba(4,4,12,0.28)"; // UTRs: exons outside the CDS
        for (const [a, b] of tx.exons) {
          if (a < tx.cds[0]) for (const r of this.rects(a, Math.min(b, tx.cds[0]), L)) octx.fillRect(...r);
          if (b > tx.cds[1]) for (const r of this.rects(Math.max(a, tx.cds[1]), b, L)) octx.fillRect(...r);
        }
      }
    }
    octx.globalCompositeOperation = "source-over";
    this.ctx.drawImage(o, 0, 0, o.width, o.height, 0, 0, this.css.w, this.css.h);
    const { ctx } = this;
    const glow = 0.55 + 0.45 * Math.sin((1 - this.pulse) * Math.PI * 6) * this.pulse;
    ctx.lineWidth = 2;
    ctx.strokeStyle = `rgba(255, 62, 165, ${glow.toFixed(3)})`;
    outline(ctx, geneRects, 1);
    if (tx) {
      // exons: boxed when big enough, bracketed above and below when only a pixel or two tall
      ctx.lineWidth = 1;
      ctx.strokeStyle = "rgba(255, 225, 77, 0.85)";
      ctx.fillStyle = "rgba(255, 225, 77, 0.9)";
      for (const [a, b] of tx.exons)
        for (const r of this.rects(a, b, L)) {
          if (r[2] >= 3 && r[3] >= 4) ctx.strokeRect(r[0] + 0.5, r[1] + 0.5, r[2] - 1, r[3] - 1);
          else {
            ctx.fillRect(r[0], r[1] - 1, Math.max(1, r[2]), 1);
            ctx.fillRect(r[0], r[1] + r[3], Math.max(1, r[2]), 1);
          }
        }
    }
    if (geneRects.length) {
      const [x, y] = geneRects[0];
      this.labels.push([g.name, Math.max(2, x), Math.max(2, y - 2), g.strand, true]);
    }
  }

  drawLive(L, ws, we) {
    const { ctx } = this;
    const c = this.contig.name;
    const reg = this.live.items("regulatory", c, ws, we);
    if (reg.length) {
      if (L.rowH < 6) {
        // zoomed out: a highlighter tint over the stretch each element covers
        ctx.lineWidth = 1;
        for (const f of reg) {
          const col = REG_COL[f.type] || "rgba(200,200,255,0.9)";
          const rs = this.rects(f.start, f.end, L);
          ctx.globalAlpha = 0.45;
          ctx.fillStyle = col;
          for (const [x, y, w, h] of rs) ctx.fillRect(x, y, w, h);
          ctx.globalAlpha = 1;
          ctx.strokeStyle = col;
          outline(ctx, rs, 0.5);
        }
      } else {
        // zoomed in: dashed outlines on a dark edge, so they show on any base colour
        const lw = L.rowH >= 6 ? 2 : 1;
        for (const f of reg) {
          const rs = this.rects(f.start, f.end, L);
          ctx.setLineDash([]);
          ctx.lineWidth = lw + 2;
          ctx.strokeStyle = "rgba(0,0,0,0.8)";
          outline(ctx, rs, lw / 2 + 1);
          ctx.setLineDash([6, 4]);
          ctx.lineWidth = lw;
          ctx.strokeStyle = REG_COL[f.type] || "rgba(200,200,255,0.9)";
          outline(ctx, rs, lw / 2 + 1);
        }
        ctx.setLineDash([]);
      }
    }
    // variants: a flag in the corner of the first base, and a line under the rest of
    // longer ones; red = clinically significant. dbSNP has a variant every few bases,
    // so when zoomed out only the red ones are shown.
    const vars = this.live.items("variation", c, ws, we);
    for (const v of vars) {
      if (L.rowH < 3 && !v.severe) continue;
      const rs = this.rects(v.start, v.end, L);
      if (!rs.length) continue;
      const fill = v.severe ? "#ff2a2a" : "#ffffff";
      const [x, y, , h] = rs[0];
      if (L.rowH >= 6) {
        if (v.end - v.start > 1) {
          const t = Math.max(2, Math.round(h / 10));
          for (const [rx, ry, rw, rh] of rs) {
            ctx.fillStyle = "#000";
            ctx.fillRect(rx, ry + rh - t - 2, rw, t + 2);
            ctx.fillStyle = fill;
            ctx.fillRect(rx + 1, ry + rh - t - 1, rw - 2, t);
          }
        }
        const k = Math.max(4, Math.round(h * (v.severe ? 0.5 : 0.36)));
        ctx.beginPath();
        ctx.moveTo(x, y);
        ctx.lineTo(x + k, y);
        ctx.lineTo(x, y + k);
        ctx.closePath();
        ctx.fillStyle = fill;
        ctx.fill();
        ctx.lineWidth = 1;
        ctx.strokeStyle = "#000";
        ctx.stroke();
      } else {
        const k = Math.max(1, Math.floor(h / 2));
        ctx.fillStyle = "#000";
        ctx.fillRect(x, y, k + 1, k + 1);
        ctx.fillStyle = fill;
        ctx.fillRect(x, y, k, k);
      }
    }
  }


  drawLetters(L) {
    const { ctx } = this;
    const top = this.state.scrollRow;
    const r0 = Math.max(0, Math.floor(top));
    const r1 = Math.min(this.totalRows(L), Math.ceil(top + this.visibleRows(L)));
    const sx = this.state.wrap === "fixed" ? this.state.scrollX : 0;
    const c0 = Math.max(0, Math.floor(sx / L.cellW));
    const c1 = Math.min(L.w, Math.ceil((sx + this.css.w) / L.cellW));
    const codes = this.codesFor(r0 * L.w, Math.min(this.contig.length, r1 * L.w));
    if (!codes) return;
    ctx.font = `${Math.round(L.s * 0.62)}px "VT323", ui-monospace, monospace`;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    const { rgb } = this.lut;
    for (let r = r0; r < r1; r++) {
      for (let c = c0; c < c1; c++) {
        const pos = r * L.w + c;
        if (pos >= this.contig.length) break;
        const code = codes[pos - r0 * L.w];
        let ch = IUPAC[code & 15];
        if (code & 0x10) ch = ch.toLowerCase();
        ctx.fillStyle = inkOn([rgb[3 * code], rgb[3 * code + 1], rgb[3 * code + 2]]);
        ctx.fillText(ch, c * L.cellW - sx + L.cellW / 2, (r - top) * L.rowH + L.rowH / 2 + 1);
      }
    }
  }

  drawGrid(L) {
    const { ctx } = this;
    const sx = this.state.wrap === "fixed" ? this.state.scrollX : 0;
    ctx.strokeStyle = "rgba(0,0,0,0.35)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    const top = this.state.scrollRow;
    const rowsVisible = this.visibleRows(L);
    for (let i = Math.floor(top); i <= top + rowsVisible + 1; i++) {
      const y = Math.round((i - top) * L.rowH) + 0.5;
      ctx.moveTo(0, y);
      ctx.lineTo(Math.min(this.css.w, L.rowPx - sx), y);
    }
    for (let c = 0; c <= L.w; c++) {
      const x = Math.round(c * L.cellW - sx) + 0.5;
      if (x < 0) continue;
      if (x > this.css.w) break;
      ctx.moveTo(x, 0);
      ctx.lineTo(x, this.css.h);
    }
    ctx.stroke();
  }

  drawCpG(L, ws, we) {
    const codes = this.codesFor(ws, we);
    if (!codes) return;
    const { ctx } = this;
    ctx.fillStyle = "rgba(255,255,255,0.9)";
    const d = Math.max(1, Math.floor(L.cellW / 3));
    for (let i = 0; i + 1 < codes.length; i++) {
      if ((codes[i] & 15) === 2 && (codes[i + 1] & 15) === 4) {
        for (const p of [ws + i, ws + i + 1]) {
          const r = this.rects(p, p + 1, L)[0];
          if (r) ctx.fillRect(r[0] + (r[2] - d) / 2, r[1] + (r[3] - d) / 2, d, d);
        }
      }
    }
  }

  drawGutter(L) {
    const g = this.gutter;
    if (!g) return;
    const gctx = g.getContext("2d");
    const { dpr, h } = this.css;
    const gw = g.width / dpr;
    gctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    gctx.fillStyle = "#0b0b18";
    gctx.fillRect(0, 0, gw, h);
    const c = this.contig;
    if (!c) return;
    const top = this.state.scrollRow;
    const total = this.totalRows(L);
    const r0 = Math.max(0, Math.floor(top));
    const r1 = Math.min(total, Math.ceil(top + this.visibleRows(L)));
    // tracks at the right edge: GC content of each row, genes (exons bright), and
    // live data when switched on (regulatory elements, clinically significant variants)
    const live = this.live && (this.live.enabled("regulatory") || this.live.enabled("variation"));
    const tx = gw - (live ? 18 : 12);
    const rowY = (r) => (r - top) * L.rowH;
    const rh = Math.max(1, L.rowH);
    if (this.rowGC) {
      const { r0: g0, values } = this.rowGC;
      for (let r = r0; r < r1; r++) {
        const v = values[r - g0];
        if (v === undefined || Number.isNaN(v)) continue;
        const t = clamp((v - 0.3) / 0.3, 0, 1); // 30% (blue) .. 60% (orange) GC
        gctx.fillStyle = `rgb(${Math.round(40 + 215 * t)},${Math.round(90 + 40 * t)},${Math.round(220 - 180 * t)})`;
        gctx.fillRect(tx, rowY(r), 5, rh);
      }
    }
    const mark = new Uint8Array(Math.max(0, r1 - r0));
    for (const gene of this.genesInView) {
      const coding = gene.biotype === "protein_coding";
      const a = Math.max(r0, Math.floor(gene.start / L.w));
      const b = Math.min(r1 - 1, Math.floor((gene.end - 1) / L.w));
      for (let r = a; r <= b; r++) mark[r - r0] = Math.max(mark[r - r0], coding ? 2 : 1);
      if (gene.tx && coding)
        for (const [ea, eb] of gene.tx.exons) {
          const x = Math.max(r0, Math.floor(ea / L.w));
          const y = Math.min(r1 - 1, Math.floor((eb - 1) / L.w));
          for (let r = x; r <= y; r++) mark[r - r0] = 3;
        }
    }
    const geneCol = ["", "#3c3c78", "#7a6ac8", "#ffe14d"];
    for (let i = 0; i < mark.length; i++) {
      if (!mark[i]) continue;
      gctx.fillStyle = geneCol[mark[i]];
      gctx.fillRect(tx + 6, rowY(r0 + i), 5, rh);
    }
    if (live) {
      const ws = r0 * L.w;
      const we = r1 * L.w;
      const lm = new Array(Math.max(0, r1 - r0)).fill(null);
      const rank = (t) => LIVE_RANK[t] || 0;
      for (const f of this.live.items("regulatory", c.name, ws, we)) {
        const a = Math.max(r0, Math.floor(f.start / L.w));
        const b = Math.min(r1 - 1, Math.floor((f.end - 1) / L.w));
        for (let r = a; r <= b; r++) if (rank(f.type) >= rank(lm[r - r0])) lm[r - r0] = f.type;
      }
      for (const v of this.live.items("variation", c.name, ws, we)) {
        const r = Math.floor(v.start / L.w);
        if (v.severe && r >= r0 && r < r1) lm[r - r0] = "severe";
      }
      for (let i = 0; i < lm.length; i++) {
        if (!lm[i]) continue;
        gctx.fillStyle = lm[i] === "severe" ? "#ff2a2a" : REG_COL[lm[i]] || "rgba(200,200,255,0.9)";
        gctx.fillRect(tx + 12, rowY(r0 + i), 5, rh);
      }
    }
    // row labels: the position of each labelled row's first base
    const step = Math.max(1, Math.ceil(GUTTER_STEP / L.rowH));
    const nice = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000];
    const k = nice.find((v) => v >= step) || step;
    gctx.font = '15px "VT323", ui-monospace, monospace';
    gctx.textAlign = "right";
    gctx.textBaseline = "top";
    const room = tx - 10;
    const first = Math.max(0, Math.ceil(top / k) * k);
    for (let r = first; r < Math.min(total, top + this.visibleRows(L)); r += k) {
      const y = rowY(r);
      const pos = r * L.w + 1;
      gctx.fillStyle = "#3a3a66";
      gctx.fillRect(tx - 6, Math.round(y), 5, 1);
      gctx.fillStyle = "#9aa0d0";
      let text = fmt(pos);
      if (L.w >= 1e5 || gctx.measureText(text).width > room) text = pos >= 1e6 ? `${(pos / 1e6).toFixed(L.w >= 1e5 ? 1 : 3)}M` : `${(pos / 1e3).toFixed(1)}k`;
      gctx.fillText(text, tx - 8, Math.max(0, y - 1));
    }
  }

  // ------------------------------------------------------------------ input
  bindEvents() {
    const cv = this.canvas;
    let drag = null;
    const pt = (e) => {
      const b = cv.getBoundingClientRect();
      return [e.clientX - b.left, e.clientY - b.top];
    };
    cv.addEventListener(
      "wheel",
      (e) => {
        e.preventDefault();
        const [x, y] = pt(e);
        if (e.ctrlKey || e.metaKey) {
          this.wheelZoom = (this.wheelZoom || 0) + e.deltaY;
          if (Math.abs(this.wheelZoom) >= 40) {
            this.zoomAt(this.wheelZoom > 0 ? -1 : 1, x, y);
            this.wheelZoom = 0;
          }
          return;
        }
        let dx = e.shiftKey ? e.deltaY : e.deltaX;
        let dy = e.shiftKey ? 0 : e.deltaY;
        if (e.deltaMode === 1) {
          dx *= 16;
          dy *= 16;
        }
        this.scrollPx(dx, dy);
      },
      { passive: false }
    );
    const down = new Set(); // pointers touching the canvas
    let multi = false; // more than one at once in this gesture: a pinch, never a click
    // a selection covers whole pixels: zoomed out, a pixel is several bases
    const span = (a, b) => [Math.min(a, b), Math.min(this.contig.length, Math.max(a, b) + this.layout().bpp)];
    cv.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      down.add(e.pointerId);
      if (down.size > 1) {
        multi = true;
        if (drag) drag.moved = true;
        return;
      }
      multi = false;
      cv.setPointerCapture(e.pointerId);
      const [x, y] = pt(e);
      drag = { x, y, sx: x, sy: y, moved: false, select: e.shiftKey, pos: this.baseAt(x, y, true), id: e.pointerId };
      if (drag.select) this.selection = span(drag.pos, drag.pos);
    });
    cv.addEventListener("pointermove", (e) => {
      const [x, y] = pt(e);
      if (drag && drag.id === e.pointerId && !this.pinch) {
        const dx = x - drag.x;
        const dy = y - drag.y;
        if (Math.abs(x - drag.sx) + Math.abs(y - drag.sy) > 3) drag.moved = true;
        if (drag.select) {
          const p = this.baseAt(x, y, true);
          this.selection = span(drag.pos, p);
          this.emit("select", this.selection);
          this.schedule();
        } else if (drag.moved) {
          this.scrollPx(-dx, -dy);
        }
        drag.x = x;
        drag.y = y;
        return;
      }
      const pos = this.baseAt(x, y);
      if (pos !== this.hover) {
        this.hover = pos;
        this.emit("hover", pos, x, y);
        this.schedule();
      }
    });
    const end = (e) => {
      down.delete(e.pointerId);
      if (!drag || drag.id !== e.pointerId) return;
      const [x, y] = pt(e);
      if (!drag.moved && !drag.select && !multi) this.emit("click", this.baseAt(x, y), x, y);
      if (drag.select) this.emit("select", this.selection);
      drag = null;
    };
    cv.addEventListener("pointerup", end);
    cv.addEventListener("pointercancel", (e) => {
      down.delete(e.pointerId);
      if (drag && drag.id === e.pointerId) drag = null;
    });
    cv.addEventListener("pointerleave", () => {
      if (!drag) {
        this.hover = null;
        this.emit("hover", null);
        this.schedule();
      }
    });
    cv.addEventListener("dblclick", (e) => {
      const [x, y] = pt(e);
      this.zoomAt(e.shiftKey ? -2 : 2, x, y);
    });
    // pinch to zoom on touch screens
    const touches = new Map();
    cv.addEventListener("touchstart", (e) => {
      for (const t of e.changedTouches) touches.set(t.identifier, [t.clientX, t.clientY]);
      if (touches.size === 2) {
        const [a, b] = [...touches.values()];
        this.pinch = { d: Math.hypot(a[0] - b[0], a[1] - b[1]) };
      }
    });
    cv.addEventListener(
      "touchmove",
      (e) => {
        for (const t of e.changedTouches) touches.set(t.identifier, [t.clientX, t.clientY]);
        if (this.pinch && touches.size === 2) {
          e.preventDefault();
          const [a, b] = [...touches.values()];
          const d = Math.hypot(a[0] - b[0], a[1] - b[1]);
          const box = cv.getBoundingClientRect();
          const cx = (a[0] + b[0]) / 2 - box.left;
          const cy = (a[1] + b[1]) / 2 - box.top;
          if (d / this.pinch.d > 1.35) {
            this.zoomAt(1, cx, cy);
            this.pinch.d = d;
          } else if (d / this.pinch.d < 0.74) {
            this.zoomAt(-1, cx, cy);
            this.pinch.d = d;
          }
        }
      },
      { passive: false }
    );
    const tend = (e) => {
      for (const t of e.changedTouches) touches.delete(t.identifier);
      if (touches.size < 2) this.pinch = null;
    };
    cv.addEventListener("touchend", tend);
    cv.addEventListener("touchcancel", tend);
  }

  setHighlight(h) {
    this.highlight = h;
    this.pulse = h ? 1 : 0;
    this.schedule();
  }

  setSelection(sel) {
    this.selection = sel;
    this.schedule();
  }

  /** PNG of what is on screen. */
  snapshot() {
    return new Promise((res) => this.canvas.toBlob(res, "image/png"));
  }

  statusText() {
    const L = this.layout();
    const zoom = L.s >= 1 ? `${L.s} px per base` : `${fmt(L.bpp)} bases per pixel`;
    return `${zoom} · ${fmtBp(L.w)} per row`;
  }
}

function outline(ctx, rects, inset = 0.5) {
  // outline of the shape an interval makes: one box, or the shape of a text selection
  if (!rects.length) return;
  ctx.beginPath();
  const box = ([x, y, w, h]) => ctx.rect(x + inset, y + inset, Math.max(1, w - 2 * inset), Math.max(1, h - 2 * inset));
  if (rects.length === 1) box(rects[0]);
  else {
    const [xA, yA, wA, hA] = rects[0];
    const [x0, yB, wB, hB] = rects[rects.length - 1];
    const R = xA + wA;
    const xB = x0 + wB;
    if (rects.length === 2 && xB <= xA) {
      box(rects[0]);
      box(rects[1]);
    } else {
      const i = inset;
      ctx.moveTo(xA + i, yA + i);
      ctx.lineTo(R - i, yA + i);
      ctx.lineTo(R - i, yB + i);
      ctx.lineTo(xB - i, yB + i);
      ctx.lineTo(xB - i, yB + hB - i);
      ctx.lineTo(x0 + i, yB + hB - i);
      ctx.lineTo(x0 + i, yA + hA - i);
      ctx.lineTo(xA + i, yA + hA - i);
      ctx.closePath();
    }
  }
  ctx.stroke();
}

function labelText(text, strand) {
  return `${strand < 0 ? "◀ " : ""}${text}${strand > 0 ? " ▶" : ""}`;
}

function labelWidth(ctx, text, strand, strong = false) {
  ctx.font = strong ? '18px "VT323", ui-monospace, monospace' : '15px "VT323", ui-monospace, monospace';
  return ctx.measureText(labelText(text, strand)).width + 6;
}

function label(ctx, text, x, y, strand, strong = false) {
  const w = labelWidth(ctx, text, strand, strong);
  ctx.textAlign = "left";
  ctx.textBaseline = "bottom";
  ctx.fillStyle = strong ? "rgba(255, 62, 165, 0.92)" : "rgba(7, 7, 15, 0.8)";
  ctx.fillRect(x, y - 16, w, 16);
  ctx.fillStyle = strong ? "#fff" : "#e8e8ff";
  ctx.fillText(labelText(text, strand), x + 3, y);
}
