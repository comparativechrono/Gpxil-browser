// The karyotype: every chromosome as a column of pixels (one pixel per
// 131,072 bases, coloured by the active palette from the summary bins),
// with its cytobands alongside, like a banded ideogram made of the genome.

import { buildLut, binColour } from "./palettes.js";
import { Emitter, clamp, fmtBp, prefersReducedMotion } from "./util.js";

export const STAIN = {
  gneg: [236, 236, 246],
  gpos25: [190, 190, 206],
  gpos50: [140, 140, 160],
  gpos75: [90, 90, 112],
  gpos100: [40, 40, 56],
  acen: [255, 62, 165],
  gvar: [120, 170, 255],
  stalk: [160, 120, 220],
};

function decodeB64(s) {
  const bin = atob(s);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

export function bandAt(bands, pos) {
  if (!bands) return null;
  let lo = 0;
  let hi = bands.length - 1;
  while (lo <= hi) {
    const m = (lo + hi) >> 1;
    if (bands[m][1] <= pos) lo = m + 1;
    else if (bands[m][0] > pos) hi = m - 1;
    else return bands[m];
  }
  return null;
}

export class Karyotype extends Emitter {
  constructor({ canvas, genome, karyotype, cytobands, genes, overview }) {
    super();
    this.canvas = canvas;
    this.genome = genome;
    this.overview = overview;
    this.ring = null; // chrM is far smaller than one karyotype pixel: drawn as a ring of its own bins
    this.data = karyotype; // {bin, bins: {chr: base64}}
    this.cytobands = cytobands;
    this.genes = genes;
    this.palette = "classic";
    this.chroms = genome.manifest.sequence.contigs.filter((c) => c.group === "chromosome");
    this.bins = new Map(this.chroms.map((c) => [c.name, decodeB64(this.data.bins[c.name] || "")]));
    this.density = null;
    this.hover = null;
    this.reveal = 1;
    this.boxes = [];
    this.bindEvents();
  }

  setPalette(p) {
    this.palette = p;
    this.draw();
  }

  async loadDensity() {
    // genes per row, from the search index (counted once)
    const idx = await this.genes.index();
    const d = new Map();
    for (const c of this.chroms) d.set(c.name, new Uint16Array(Math.ceil(c.length / this.data.bin) + 1));
    for (const r of idx.rows) {
      const a = d.get(r[2]);
      if (a && r[6] === "protein_coding") a[Math.floor((r[3] + r[4]) / 2 / this.data.bin)]++;
    }
    this.density = d;
    this.draw();
  }

  show() {
    this.reveal = prefersReducedMotion() ? 1 : 0;
    const t0 = performance.now();
    const step = () => {
      this.reveal = Math.min(1, (performance.now() - t0) / 900);
      this.draw();
      if (this.reveal < 1) requestAnimationFrame(step);
    };
    if (this.reveal < 1) requestAnimationFrame(step);
    else this.draw();
    if (!this.density) this.loadDensity().catch(() => {});
    if (!this.ring) this.loadRings().catch(() => {});
  }

  async loadRings() {
    // contigs shorter than one karyotype bin (chrM) get a ring of 128-base bins instead
    const ring = new Map();
    for (const c of this.chroms) if (c.length < this.data.bin && this.overview) ring.set(c.name, await this.overview.level(c, 0));
    this.ring = ring;
    this.draw();
  }

  layout() {
    const box = this.canvas.parentElement.getBoundingClientRect();
    const W = Math.max(300, box.width);
    const H = Math.max(200, box.height);
    const n = this.chroms.length;
    const longest = Math.max(...this.chroms.map((c) => c.length));
    const binsLongest = Math.ceil(longest / this.data.bin);
    // one or two lines of chromosomes, pixel columns 4-16 bins wide
    let best = null;
    for (const lines of [1, 2]) {
      const perLine = Math.ceil(n / lines);
      const slot = W / perLine;
      const avail = (H - 34 * lines - 24) / lines;
      for (const per of [16, 12, 8, 6, 4, 3, 2]) {
        const rows = Math.ceil(binsLongest / per);
        const px = Math.floor(avail / rows);
        if (px < 1) continue;
        const width = per * px + 8; // + band stripe
        if (width > slot - 8) continue;
        const score = rows * px * per * px; // area of the longest chromosome
        if (!best || score > best.score) best = { lines, perLine, slot, per, px, rows, score, lineH: avail + 34 };
      }
    }
    return best || { lines: 2, perLine: Math.ceil(n / 2), slot: W / Math.ceil(n / 2), per: 2, px: 1, lineH: (H - 24) / 2 };
  }

  draw() {
    const cv = this.canvas;
    const box = cv.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    cv.width = Math.round(box.width * dpr);
    cv.height = Math.round(box.height * dpr);
    cv.style.width = `${box.width}px`;
    cv.style.height = `${box.height}px`;
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.fillStyle = "#07070f";
    ctx.fillRect(0, 0, box.width, box.height);
    const L = this.layout();
    this.L = L;
    const lut = buildLut(this.palette);
    const col = [0, 0, 0];
    this.boxes = [];
    ctx.font = '20px "VT323", ui-monospace, monospace';
    ctx.textAlign = "center";
    this.chroms.forEach((c, k) => {
      const line = Math.floor(k / L.perLine);
      const i = k % L.perLine;
      const x0 = Math.round(i * L.slot + (L.slot - (L.per * L.px + 8)) / 2);
      const y0 = Math.round(line * L.lineH + 12);
      if (c.length < this.data.bin) return this.drawRing(ctx, c, x0, y0, L, lut);
      const bins = this.bins.get(c.name);
      const nb = bins.length / 2;
      const rows = Math.max(3, Math.ceil(nb / L.per));
      const shown = Math.ceil(rows * this.reveal);
      const img = ctx.createImageData(L.per, rows);
      for (let b = 0; b < nb; b++) {
        const r = Math.floor(b / L.per);
        if (r >= shown) break;
        const byte = bins[2 * b + 1];
        binColour(lut, bins[2 * b] / 255, (byte >> 4) / 15, (byte & 15) / 15, col, 2.6); // stretched so isochores show
        const o = 4 * (r * L.per + (b % L.per));
        img.data[o] = col[0];
        img.data[o + 1] = col[1];
        img.data[o + 2] = col[2];
        img.data[o + 3] = 255;
      }
      const tmp = document.createElement("canvas");
      tmp.width = L.per;
      tmp.height = rows;
      tmp.getContext("2d").putImageData(img, 0, 0);
      const cx = x0 + 8;
      ctx.drawImage(tmp, cx, y0, L.per * L.px, rows * L.px);
      // cytoband stripe
      const bands = this.cytobands[c.name];
      if (bands) {
        for (let r = 0; r < shown; r++) {
          const pos = (r + 0.5) * L.per * this.data.bin;
          const b = bandAt(bands, Math.min(pos, c.length - 1));
          if (!b) continue;
          const s = STAIN[b[3]] || STAIN.gneg;
          ctx.fillStyle = `rgb(${s[0]},${s[1]},${s[2]})`;
          ctx.fillRect(x0 + 1, y0 + r * L.px, 5, L.px);
        }
      }
      // gene density: a thin heat line on the right
      if (this.density && this.density.get(c.name)) {
        const d = this.density.get(c.name);
        for (let r = 0; r < shown; r++) {
          let sum = 0;
          for (let b = r * L.per; b < Math.min(nb, (r + 1) * L.per); b++) sum += d[b];
          if (!sum) continue;
          const a = clamp(sum / (L.per * 1.6), 0.08, 1);
          ctx.fillStyle = `rgba(255, 225, 77, ${a.toFixed(2)})`;
          ctx.fillRect(cx + L.per * L.px + 1, y0 + r * L.px, 2, L.px);
        }
      }
      const h = rows * L.px;
      const label = c.name.replace(/^chr/, "");
      ctx.fillStyle = this.hover && this.hover.contig === c.name ? "#ff3ea5" : "#c8c8f0";
      ctx.fillText(label, cx + (L.per * L.px) / 2, y0 + h + 20);
      this.boxes.push({ contig: c, x: cx, y: y0, w: L.per * L.px, h, labelY: y0 + h, bandX: x0 });
      if (this.hover && this.hover.contig === c.name && this.hover.pos !== null) {
        const r = Math.floor(this.hover.pos / this.data.bin / L.per);
        ctx.strokeStyle = "#ff3ea5";
        ctx.lineWidth = 2;
        ctx.strokeRect(cx - 3, y0 + r * L.px - 2, L.per * L.px + 6, L.px + 4);
      }
    });
  }

  /** A circular genome (chrM) as a ring of pixels, starting at the top and running clockwise. */
  drawRing(ctx, c, x0, y0, L, lut) {
    const cell = Math.max(1, Math.min(3, L.px));
    const D = Math.max(26, L.per * L.px);
    const R = D / 2;
    const cx = x0 + 8 + (L.per * L.px) / 2;
    const cy = y0 + R;
    const thick = Math.max(4, Math.round(D * 0.2));
    const bins = this.ring && this.ring.get(c.name);
    const col = [0, 0, 0];
    const n = bins ? bins.length / 2 : 0;
    const reveal = this.reveal;
    for (let yy = -R; yy < R; yy += cell) {
      for (let xx = -R; xx < R; xx += cell) {
        const px = xx + cell / 2;
        const py = yy + cell / 2;
        const d = Math.hypot(px, py);
        if (d > R || d < R - thick) continue;
        let a = Math.atan2(px, -py) / (2 * Math.PI);
        if (a < 0) a += 1;
        if (a > reveal) continue;
        if (n) {
          const i = Math.min(n - 1, Math.floor(a * n));
          const byte = bins[2 * i + 1];
          binColour(lut, bins[2 * i] / 255, (byte >> 4) / 15, (byte & 15) / 15, col, 2.6);
          ctx.fillStyle = `rgb(${col[0] | 0},${col[1] | 0},${col[2] | 0})`;
        } else ctx.fillStyle = "#3a3a66";
        ctx.fillRect(Math.round(cx + xx), Math.round(cy + yy), cell, cell);
      }
    }
    const hovered = this.hover && this.hover.contig === c.name;
    ctx.fillStyle = hovered ? "#ff3ea5" : "#c8c8f0";
    ctx.fillText(c.name.replace(/^chr/, ""), cx, y0 + D + 20);
    this.boxes.push({ contig: c, x: cx - R, y: y0, w: D, h: D, labelY: y0 + D, bandX: cx - R, ring: { cx, cy, R } });
    if (hovered && this.hover.pos !== null) {
      const a = (this.hover.pos / c.length) * 2 * Math.PI;
      const r = R - thick / 2;
      ctx.strokeStyle = "#ff3ea5";
      ctx.lineWidth = 2;
      ctx.strokeRect(Math.round(cx + Math.sin(a) * r) - 4, Math.round(cy - Math.cos(a) * r) - 4, 8, 8);
    }
  }

  hit(x, y) {
    for (const b of this.boxes) {
      if (b.ring) {
        if (x < b.x - 6 || x > b.x + b.w + 6 || y < b.y - 6 || y > b.labelY + 26) continue;
        const { cx, cy } = b.ring;
        let a = Math.atan2(x - cx, cy - y) / (2 * Math.PI);
        if (a < 0) a += 1;
        const onLabel = y > b.labelY + 2;
        return { contig: b.contig.name, pos: onLabel ? null : clamp(Math.floor(a * b.contig.length), 0, b.contig.length - 1), length: b.contig.length, circular: true };
      }
      if (x >= b.bandX - 4 && x <= b.x + b.w + 6 && y >= b.y - 4 && y <= b.labelY + 26) {
        const onLabel = y > b.labelY + 2;
        const r = clamp(Math.floor((y - b.y) / this.L.px), 0, Math.ceil(b.contig.length / this.data.bin / this.L.per) - 1);
        const pos = onLabel ? null : clamp((r + 0.5) * this.L.per * this.data.bin, 0, b.contig.length - 1);
        return { contig: b.contig.name, pos, length: b.contig.length };
      }
    }
    return null;
  }

  bindEvents() {
    const cv = this.canvas;
    const pt = (e) => {
      const b = cv.getBoundingClientRect();
      return [e.clientX - b.left, e.clientY - b.top];
    };
    cv.addEventListener("pointermove", (e) => {
      const h = this.hit(...pt(e));
      const key = h ? `${h.contig}:${h.pos}` : null;
      if (key !== (this.hover ? `${this.hover.contig}:${this.hover.pos}` : null)) {
        this.hover = h;
        cv.style.cursor = h ? "pointer" : "default";
        this.draw();
        if (h) {
          const band = h.pos !== null ? bandAt(this.cytobands[h.contig], h.pos) : null;
          const where = h.pos !== null ? ` · ${h.circular ? fmtBp(h.pos + 1) : fmtBp(h.pos)}${band ? ` · ${h.contig.replace("chr", "")}${band[2]}` : ""}` : ` · ${fmtBp(h.length)}`;
          this.emit("hover", `${h.contig}${where}${h.circular ? " · circular" : ""}`, e.clientX, e.clientY);
        } else this.emit("hover", null);
      }
    });
    cv.addEventListener("pointerleave", () => {
      this.hover = null;
      this.draw();
      this.emit("hover", null);
    });
    cv.addEventListener("click", (e) => {
      const h = this.hit(...pt(e));
      if (h) this.emit("open", h.contig, h.pos);
    });
    new ResizeObserver(() => {
      if (this.canvas.offsetParent) this.draw();
    }).observe(cv.parentElement);
  }
}
