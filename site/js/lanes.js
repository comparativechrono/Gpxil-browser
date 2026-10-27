// Expanded mode: every transcript of a gene on its own lane, drawn with the
// gene's own pixels, like the expanded gene track of a classic browser.
// Coding exons fill the lane, UTRs and non-coding exons half of it, introns
// are a thin line, so alternative exons, starts and ends line up down the
// lanes. The top lane is the DNA itself. Click a lane to open that transcript.

import { Emitter, clamp, fmt } from "./util.js";

const RULER = 24; // height of the position ruler
const PAD = 8; // space between the ruler and the first lane

export class TranscriptLanes extends Emitter {
  constructor({ canvas, genome, lut, gutterWidth }) {
    super();
    this.canvas = canvas;
    this.genome = genome;
    this.lutOf = lut; // () => the viewer's colours (palette and repeats)
    this.gutterWidth = gutterWidth; // () => width of the label column, matching the viewer's gutter
    this.ctx = canvas.getContext("2d");
    this.off = document.createElement("canvas");
    this.active = false;
    this.gene = null;
    this.txs = [];
    this.selected = null;
    this.a = 0; // window [a, b) in genome coordinates
    this.b = 1;
    this.scrollY = 0;
    this.bases = null; // {start, end, codes} around the gene
    this.hover = null; // {lane, pos}
    this.css = { w: 1, h: 1, dpr: 1 };
    this.raf = 0;
    this.bind();
    new ResizeObserver(() => this.active && this.resize()).observe(canvas.parentElement);
  }

  // ------------------------------------------------------------------ state
  show(gene, txs, selectedId) {
    const fresh = !this.gene || this.gene.id !== gene.id;
    this.gene = gene;
    this.txs = txs;
    this.selected = selectedId;
    if (fresh) {
      this.scrollY = 0;
      this.hover = null;
      this.bases = null;
      this.loadBases();
    }
    if (!this.active) {
      this.active = true;
      this.canvas.hidden = false;
      this.resize(true);
    }
    if (fresh) this.reset();
    else {
      this.scrollY = clamp(this.scrollY, 0, this.maxScroll());
      this.schedule();
      this.emit("navigate", this.describe());
    }
  }

  hide() {
    if (!this.active) return;
    this.active = false;
    this.canvas.hidden = true;
    this.emit("hover", null);
  }

  /** Bases for the gene and the room around it that zooming out can show. */
  bounds() {
    const g = this.gene;
    const c = this.genome.contig(g.contig);
    const len = g.end - g.start;
    return [Math.max(0, g.start - len), Math.min(c.length, g.end + len)];
  }

  async loadBases() {
    const g = this.gene;
    const [lo, hi] = this.bounds();
    try {
      const codes = await this.genome.bases(g.contig, lo, hi);
      if (this.gene !== g) return;
      this.bases = { start: lo, end: hi, codes };
      this.schedule();
    } catch (e) {
      this.emit("error", e);
    }
  }

  /** The whole gene, with a little room on either side. */
  reset() {
    const g = this.gene;
    const m = Math.max(10, Math.round((g.end - g.start) * 0.03));
    this.setWindow(g.start - m, g.end + m);
  }

  setWindow(a, b) {
    const [lo, hi] = this.bounds();
    const W = this.lanesWidth();
    const minLen = Math.max(20, W / 64); // at most 64 pixels per base
    let len = clamp(b - a, minLen, hi - lo);
    a = clamp(a, lo, hi - len);
    this.a = a;
    this.b = a + len;
    this.schedule();
    this.emit("navigate", this.describe());
  }

  describe() {
    const bpp = (this.b - this.a) / this.lanesWidth();
    return { contig: this.gene.contig, start: Math.floor(this.a), end: Math.ceil(this.b), bpp, transcripts: this.txs.length };
  }

  zoom(steps, x = null) {
    const W = this.lanesWidth();
    const fx = x === null ? 0.5 : clamp((x - this.labelW) / W, 0, 1);
    const at = this.a + (this.b - this.a) * fx;
    const len = (this.b - this.a) / 1.5 ** steps;
    this.setWindow(at - len * fx, at - len * fx + len);
  }

  pan(px) {
    const bpp = (this.b - this.a) / this.lanesWidth();
    this.setWindow(this.a + px * bpp, this.b + px * bpp);
  }

  scroll(px) {
    this.scrollY = clamp(this.scrollY + px, 0, this.maxScroll());
    this.schedule();
  }

  centreOn(pos) {
    const half = (this.b - this.a) / 2;
    this.setWindow(pos - half, pos + half);
  }

  // ------------------------------------------------------------------ layout
  get labelW() {
    return this.gutterWidth();
  }

  lanesWidth() {
    return Math.max(10, this.css.w - this.labelW);
  }

  laneH() {
    return this.css.w < 600 ? 14 : 18;
  }

  gap() {
    return this.css.w < 600 ? 5 : 6;
  }

  /** Top of lane i (0 = the DNA, 1.. = transcripts), in canvas CSS pixels. */
  laneTop(i) {
    return RULER + PAD + i * (this.laneH() + this.gap()) - this.scrollY;
  }

  maxScroll() {
    const content = PAD + (this.txs.length + 1) * (this.laneH() + this.gap()) + 12;
    return Math.max(0, content - (this.css.h - RULER));
  }

  laneAt(y) {
    if (y < RULER) return null;
    const i = Math.floor((y - RULER - PAD + this.scrollY) / (this.laneH() + this.gap()));
    if (i < 0 || i > this.txs.length) return null;
    const inLane = y - this.laneTop(i) < this.laneH() + this.gap() / 2;
    return inLane ? i : null;
  }

  posAt(x) {
    if (x < this.labelW) return null;
    return Math.floor(this.a + ((x - this.labelW) / this.lanesWidth()) * (this.b - this.a));
  }

  resize(force = false) {
    const box = this.canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    const w = Math.max(50, Math.floor(box.width));
    const h = Math.max(50, Math.floor(box.height));
    if (!force && w === this.css.w && h === this.css.h && dpr === this.css.dpr) return;
    const keep = this.gene ? [this.a, this.b] : null;
    this.css = { w, h, dpr };
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.canvas.style.width = `${w}px`;
    this.canvas.style.height = `${h}px`;
    this.scrollY = clamp(this.scrollY, 0, this.maxScroll());
    if (keep) this.setWindow(...keep);
    this.schedule();
  }

  // ------------------------------------------------------------------ drawing
  schedule() {
    if (this.active && !this.raf) this.raf = requestAnimationFrame(() => this.draw());
  }

  /** What each column of a transcript's lane shows: 0 nothing, 1 intron, 2 UTR or non-coding exon, 3 coding exon. */
  classes(t, W) {
    const out = new Uint8Array(W);
    const bpp = (this.b - this.a) / W;
    const col = (p) => (p - this.a) / bpp;
    const fill = (s, e, v) => {
      const x0 = clamp(Math.floor(col(s)), 0, W);
      const x1 = clamp(Math.max(Math.ceil(col(e)), Math.floor(col(s)) + 1), 0, W);
      for (let x = x0; x < x1; x++) if (out[x] < v) out[x] = v;
    };
    fill(t.start, t.end, 1);
    for (const [s, e] of t.exons) {
      fill(s, e, 2);
      if (t.cds) {
        const cs = Math.max(s, t.cds[0]);
        const ce = Math.min(e, t.cds[1]);
        if (ce > cs) fill(cs, ce, 3);
      }
    }
    return out;
  }

  /** The mean colour of the bases under each column (kept until the window or the colours change). */
  columnColours(W) {
    const lut = this.lutOf();
    const d = this.bases;
    const key = `${this.a}|${this.b}|${W}|${d ? d.start : "-"}`;
    if (this.cols && this.cols.key === key && this.cols.lut === lut) return this.cols.out;
    const out = new Float32Array(3 * W).fill(70);
    if (d) this.cols = { key, lut, out };
    if (!d) return out;
    const bpp = (this.b - this.a) / W;
    for (let x = 0; x < W; x++) {
      const p0 = Math.floor(this.a + x * bpp);
      const p1 = Math.max(p0 + 1, Math.floor(this.a + (x + 1) * bpp));
      let R = 0;
      let G = 0;
      let B = 0;
      let n = 0;
      for (let p = Math.max(p0, d.start); p < Math.min(p1, d.end); p++) {
        const o = d.codes[p - d.start] * 3;
        R += lut.rgb[o];
        G += lut.rgb[o + 1];
        B += lut.rgb[o + 2];
        n++;
      }
      if (n) {
        out[3 * x] = R / n;
        out[3 * x + 1] = G / n;
        out[3 * x + 2] = B / n;
      }
    }
    return out;
  }

  draw() {
    this.raf = 0;
    if (!this.active || !this.gene) return;
    const { ctx } = this;
    const { w, h, dpr } = this.css;
    const L = this.labelW;
    const W = this.lanesWidth();
    const H = Math.max(1, h - RULER);
    const LH = this.laneH();
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.fillStyle = "#07070f";
    ctx.fillRect(0, 0, w, h);
    // the lanes, as pixels, in an offscreen image of the lanes area
    const cols = this.columnColours(W);
    const off = this.off;
    if (off.width !== W || off.height !== H) {
      off.width = W;
      off.height = H;
    }
    const octx = off.getContext("2d");
    const img = octx.createImageData(W, H);
    const px = new Uint32Array(img.data.buffer);
    const little = new Uint8Array(new Uint32Array([1]).buffer)[0] === 1;
    const pack = (r, g, b) => (little ? (255 << 24) | (b << 16) | (g << 8) | r : (r << 24) | (g << 16) | (b << 8) | 255);
    px.fill(pack(7, 7, 15));
    const shade = (x, k) => pack((7 + (cols[3 * x] - 7) * k) | 0, (7 + (cols[3 * x + 1] - 7) * k) | 0, (15 + (cols[3 * x + 2] - 15) * k) | 0);
    const rows = (y0, y1, x, c) => {
      for (let y = Math.max(0, y0); y < Math.min(H, y1); y++) px[y * W + x] = c;
    };
    const first = Math.max(0, Math.floor((this.scrollY - PAD) / (LH + this.gap())));
    for (let i = first; i <= this.txs.length; i++) {
      const top = this.laneTop(i) - RULER;
      if (top >= H) break;
      if (top + LH < 0) continue;
      if (i === 0) {
        for (let x = 0; x < W; x++) rows(top, top + LH, x, shade(x, 0.95)); // the DNA
        continue;
      }
      const t = this.txs[i - 1];
      const cls = this.classes(t, W);
      const mid = top + (LH >> 1);
      const coding = !!t.cds;
      for (let x = 0; x < W; x++) {
        const c = cls[x];
        if (c === 3) rows(top, top + LH, x, shade(x, 1));
        else if (c === 2) rows(top + Math.round(LH * 0.25), top + Math.round(LH * 0.75), x, shade(x, coding ? 0.8 : 0.95));
        else if (c === 1) rows(mid - 1, mid + 1, x, shade(x, 0.42));
      }
    }
    octx.putImageData(img, 0, 0);
    ctx.drawImage(off, 0, 0, W, H, L, RULER, W, H);
    this.drawDecorations(L, W, H, LH);
    this.drawRuler(L, W);
    this.drawLabels(L, LH);
  }

  drawDecorations(L, W, H, LH) {
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.rect(L, RULER, W, H);
    ctx.clip();
    const bpp = (this.b - this.a) / W;
    for (let i = 1; i <= this.txs.length; i++) {
      const top = this.laneTop(i);
      if (top > RULER + H || top + LH < RULER) continue;
      const t = this.txs[i - 1];
      // strand chevrons along the introns
      const cls = this.classes(t, W);
      const mid = top + (LH >> 1);
      ctx.strokeStyle = "rgba(220, 220, 255, 0.55)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (let x = 24; x < W - 4; x += 36) {
        if (cls[x] !== 1 || cls[x - 3] !== 1 || cls[x + 3] !== 1) continue;
        const d = t.strand > 0 ? 2 : -2;
        ctx.moveTo(L + x - d + 0.5, mid - 2.5);
        ctx.lineTo(L + x + d + 0.5, mid + 0.5);
        ctx.lineTo(L + x - d + 0.5, mid + 3.5);
      }
      ctx.stroke();
      // the transcript in use
      if (t.id === this.selected) {
        ctx.strokeStyle = "#ffe14d";
        ctx.lineWidth = 1;
        ctx.strokeRect(L + 0.5, top - 2.5, W - 1, LH + 5);
      } else if (this.hover && this.hover.lane === i) {
        ctx.strokeStyle = "rgba(32, 227, 255, 0.7)";
        ctx.strokeRect(L + 0.5, top - 2.5, W - 1, LH + 5);
      }
    }
    // the position under the pointer, down all lanes
    if (this.hover && this.hover.pos !== null) {
      const x = L + Math.floor((this.hover.pos - this.a) / bpp);
      ctx.fillStyle = "rgba(255, 255, 255, 0.35)";
      ctx.fillRect(x, RULER, Math.max(1, Math.round(1 / bpp)), H);
    }
    ctx.restore();
  }

  drawRuler(L, W) {
    const { ctx } = this;
    ctx.fillStyle = "#0b0b18";
    ctx.fillRect(0, 0, this.css.w, RULER);
    ctx.fillStyle = "#3a3a66";
    ctx.fillRect(L, RULER - 1, W, 1);
    const span = this.b - this.a;
    const nice = [1, 2, 5];
    let step = 1;
    for (let k = 0; step < span / Math.max(2, W / 110); k++) step = nice[k % 3] * 10 ** Math.floor(k / 3);
    ctx.font = '15px "VT323", ui-monospace, monospace';
    ctx.textBaseline = "top";
    ctx.textAlign = "left";
    for (let p = Math.ceil(this.a / step) * step; p < this.b; p += step) {
      const x = L + ((p - this.a) / span) * W;
      ctx.fillStyle = "#3a3a66";
      ctx.fillRect(Math.round(x), RULER - 7, 1, 6);
      ctx.fillStyle = "#9aa0d0";
      if (x + 70 < L + W) ctx.fillText(fmt(p + 1), Math.round(x) + 3, 4);
    }
    ctx.fillStyle = "#9aa0d0";
    ctx.textAlign = "right";
    ctx.fillText(this.gene.contig, L - 8, 4);
  }

  drawLabels(L, LH) {
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, RULER, L - 2, this.css.h - RULER);
    ctx.clip();
    ctx.fillStyle = "#0b0b18";
    ctx.fillRect(0, RULER, L - 2, this.css.h - RULER);
    ctx.fillStyle = "#1c1c38";
    ctx.fillRect(L - 2, RULER, 2, this.css.h - RULER);
    ctx.font = '15px "VT323", ui-monospace, monospace';
    ctx.textBaseline = "middle";
    ctx.textAlign = "right";
    const fit = (text, max) => {
      if (ctx.measureText(text).width <= max) return text;
      while (text.length > 2 && ctx.measureText(`${text}…`).width > max) text = text.slice(0, -1);
      return `${text}…`;
    };
    for (let i = 0; i <= this.txs.length; i++) {
      const top = this.laneTop(i);
      if (top > this.css.h) break;
      if (top + LH < RULER) continue;
      const y = top + LH / 2 + 1;
      if (i === 0) {
        ctx.fillStyle = "#ff3ea5";
        ctx.fillText(fit(`${this.gene.name} ${this.gene.strand > 0 ? "▶" : "◀"}`, L - 12), L - 8, y);
        continue;
      }
      const t = this.txs[i - 1];
      ctx.fillStyle = t.id === this.selected ? "#ffe14d" : t.mane ? "#ffe14d" : t.canonical ? "#20e3ff" : t.cds ? "#c8c8f0" : "#8a8ab0";
      ctx.fillText(fit(`${t.mane ? "★ " : ""}${t.name}`, L - 12), L - 8, y);
    }
    ctx.restore();
  }

  /** PNG of the lanes. */
  snapshot() {
    return new Promise((res) => this.canvas.toBlob(res, "image/png"));
  }

  // ------------------------------------------------------------------ input
  info(lane, pos) {
    const g = this.gene;
    if (lane === 0) return `${g.name} · the DNA · ${g.contig}:${fmt(pos + 1)}`;
    const t = this.txs[lane - 1];
    const tags = [t.mane && "MANE Select", t.maneClinical && "MANE Plus Clinical", t.canonical && !t.mane && "canonical"].filter(Boolean).join(", ");
    let where = "outside the transcript";
    if (pos !== null && pos >= t.start && pos < t.end) {
      const i = t.exons.findIndex(([s, e]) => pos >= s && pos < e);
      if (i < 0) where = "intron";
      else {
        const n = t.strand > 0 ? i + 1 : t.exons.length - i;
        const coding = t.cds && pos >= t.cds[0] && pos < t.cds[1];
        where = `exon ${n} of ${t.exons.length}${t.cds ? (coding ? " · coding" : " · UTR") : ""}`;
      }
    }
    const at = pos === null ? "" : ` · ${g.contig}:${fmt(pos + 1)} · ${where}`;
    return `${t.name} · ${t.biotype.replace(/_/g, " ")}${tags ? ` · ${tags}` : ""}${at} · click to open`;
  }

  bind() {
    const cv = this.canvas;
    const pt = (e) => {
      const b = cv.getBoundingClientRect();
      return [e.clientX - b.left, e.clientY - b.top];
    };
    const pointers = new Map();
    let drag = null;
    let pinch = null;
    cv.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      cv.setPointerCapture(e.pointerId);
      pointers.set(e.pointerId, pt(e));
      if (pointers.size === 2) {
        const [p, q] = [...pointers.values()];
        pinch = { d: Math.hypot(p[0] - q[0], p[1] - q[1]) };
        if (drag) drag.moved = true;
        return;
      }
      const [x, y] = pt(e);
      drag = { x, y, sx: x, sy: y, moved: false, id: e.pointerId };
    });
    cv.addEventListener("pointermove", (e) => {
      const [x, y] = pt(e);
      if (pointers.has(e.pointerId)) pointers.set(e.pointerId, [x, y]);
      if (pinch && pointers.size === 2) {
        const [p, q] = [...pointers.values()];
        const d = Math.hypot(p[0] - q[0], p[1] - q[1]);
        if (d / pinch.d > 1.3 || d / pinch.d < 0.77) {
          this.zoom(d > pinch.d ? 1 : -1, (p[0] + q[0]) / 2);
          pinch.d = d;
        }
        return;
      }
      if (drag && drag.id === e.pointerId) {
        if (Math.abs(x - drag.sx) + Math.abs(y - drag.sy) > 4) drag.moved = true;
        if (drag.moved) {
          this.pan(drag.x - x);
          this.scroll(drag.y - y);
        }
        drag.x = x;
        drag.y = y;
        return;
      }
      const lane = this.laneAt(y);
      const pos = this.posAt(x);
      const hover = lane === null ? null : { lane, pos };
      if (JSON.stringify(hover) !== JSON.stringify(this.hover)) {
        this.hover = hover;
        this.schedule();
        this.emit("hover", hover ? this.info(lane, pos) : null, e.clientX, e.clientY, pos);
      }
      cv.style.cursor = lane > 0 ? "pointer" : "default";
    });
    const end = (e) => {
      pointers.delete(e.pointerId);
      if (pointers.size < 2) pinch = null;
      if (!drag || drag.id !== e.pointerId) return;
      const [, y] = pt(e);
      if (!drag.moved) {
        const lane = this.laneAt(y);
        if (lane > 0) this.emit("pick", this.txs[lane - 1].id);
      }
      drag = null;
    };
    cv.addEventListener("pointerup", end);
    cv.addEventListener("pointercancel", (e) => {
      pointers.delete(e.pointerId);
      pinch = null;
      drag = null;
    });
    cv.addEventListener("pointerleave", () => {
      if (drag) return;
      this.hover = null;
      this.schedule();
      this.emit("hover", null);
    });
    cv.addEventListener(
      "wheel",
      (e) => {
        e.preventDefault();
        const [x] = pt(e);
        if (e.ctrlKey || e.metaKey) {
          this.wheelZoom = (this.wheelZoom || 0) + e.deltaY;
          if (Math.abs(this.wheelZoom) >= 40) {
            this.zoom(this.wheelZoom > 0 ? -1 : 1, x);
            this.wheelZoom = 0;
          }
          return;
        }
        const k = e.deltaMode === 1 ? 16 : 1;
        if (e.shiftKey || Math.abs(e.deltaX) > Math.abs(e.deltaY)) this.pan((e.shiftKey ? e.deltaY : e.deltaX) * k);
        else this.scroll(e.deltaY * k);
      },
      { passive: false }
    );
    cv.addEventListener("dblclick", (e) => {
      const [x] = pt(e);
      this.zoom(e.shiftKey ? -2 : 2, x);
    });
  }
}

