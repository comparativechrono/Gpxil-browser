// The chromosome bar above the viewer: cytobands, the chromosome's own pixel
// summary, and a bracket showing the part on screen. Click or drag to move.

import { buildLut, binColour } from "./palettes.js";
import { STAIN, bandAt } from "./karyotype.js";
import { Emitter, clamp, fmtBp } from "./util.js";

export class Ideogram extends Emitter {
  constructor({ canvas, cytobands, karyotype }) {
    super();
    this.canvas = canvas;
    this.cytobands = cytobands;
    this.karyotype = karyotype;
    this.contig = null;
    this.window = null;
    this.palette = "classic";
    this.markers = [];
    this.bind();
    new ResizeObserver(() => this.draw()).observe(canvas.parentElement);
  }

  set(contig, window, markers = []) {
    this.contig = contig;
    this.window = window;
    this.markers = markers;
    this.draw();
  }

  draw() {
    const cv = this.canvas;
    const box = cv.getBoundingClientRect();
    if (!box.width) return;
    const dpr = window.devicePixelRatio || 1;
    cv.width = Math.round(box.width * dpr);
    cv.height = Math.round(box.height * dpr);
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = "#0b0b18";
    ctx.fillRect(0, 0, box.width, box.height);
    const c = this.contig;
    if (!c) return;
    const pad = 8;
    const W = box.width - 2 * pad;
    const x = (pos) => pad + (pos / c.length) * W;
    this.geom = { pad, W, len: c.length };
    // cytobands
    const bands = this.cytobands[c.name];
    const yb = 6;
    const hb = 12;
    if (bands) {
      for (const [s, e, , stain] of bands) {
        const col = STAIN[stain] || STAIN.gneg;
        ctx.fillStyle = `rgb(${col[0]},${col[1]},${col[2]})`;
        const x0 = Math.floor(x(s));
        ctx.fillRect(x0, yb, Math.max(1, Math.ceil(x(e)) - x0), hb);
      }
    } else {
      ctx.fillStyle = "#3a3a66";
      ctx.fillRect(pad, yb, W, hb);
    }
    // the chromosome's own pixels (summary bins), one column per screen pixel
    const b64 = this.karyotype.bins[c.name];
    if (b64) {
      if (!this.binCache || this.binCache.name !== c.name) {
        const bin = atob(b64);
        const bytes = new Uint8Array(bin.length);
        for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
        this.binCache = { name: c.name, bytes };
      }
      const bytes = this.binCache.bytes;
      const nb = bytes.length / 2;
      const lut = buildLut(this.palette);
      const col = [0, 0, 0];
      for (let px = 0; px < W; px++) {
        const i = Math.min(nb - 1, Math.floor((px / W) * nb));
        const byte = bytes[2 * i + 1];
        binColour(lut, bytes[2 * i] / 255, (byte >> 4) / 15, (byte & 15) / 15, col, 2.6);
        ctx.fillStyle = `rgb(${col[0] | 0},${col[1] | 0},${col[2] | 0})`;
        ctx.fillRect(pad + px, yb + hb + 3, 1, 8);
      }
    }
    // markers (selected gene, selection)
    for (const m of this.markers) {
      ctx.fillStyle = m.colour;
      const x0 = x(m.start);
      ctx.fillRect(Math.floor(x0), 2, Math.max(2, Math.ceil(x(m.end) - x0)), box.height - 4);
    }
    // visible window
    if (this.window) {
      const [s, e] = this.window;
      const x0 = Math.floor(x(s));
      const x1 = Math.max(x0 + 3, Math.ceil(x(e)));
      ctx.strokeStyle = "#20e3ff";
      ctx.lineWidth = 2;
      ctx.strokeRect(x0 - 1, 2, x1 - x0 + 2, box.height - 4);
      ctx.fillStyle = "rgba(32,227,255,0.15)";
      ctx.fillRect(x0, 2, x1 - x0, box.height - 4);
    }
  }

  posAt(clientX) {
    const b = this.canvas.getBoundingClientRect();
    if (!this.geom) return null;
    return clamp(Math.floor(((clientX - b.left - this.geom.pad) / this.geom.W) * this.geom.len), 0, this.geom.len - 1);
  }

  bind() {
    const cv = this.canvas;
    let dragging = false;
    cv.addEventListener("pointerdown", (e) => {
      dragging = true;
      cv.setPointerCapture(e.pointerId);
      const p = this.posAt(e.clientX);
      if (p !== null) this.emit("seek", p);
    });
    cv.addEventListener("pointermove", (e) => {
      const p = this.posAt(e.clientX);
      if (p === null || !this.contig) return;
      if (dragging) this.emit("seek", p);
      const band = bandAt(this.cytobands[this.contig.name], p);
      this.emit("hover", `${this.contig.name}:${fmtBp(p)}${band ? ` · ${this.contig.name.replace("chr", "")}${band[2]}` : ""}`, e.clientX, e.clientY);
    });
    const stop = () => (dragging = false);
    cv.addEventListener("pointerup", stop);
    cv.addEventListener("pointercancel", stop);
    cv.addEventListener("pointerleave", () => this.emit("hover", null));
  }
}
