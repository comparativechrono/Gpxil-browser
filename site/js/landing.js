// "Enter the pixelverse": the landing screen. Everything is drawn on a
// quarter-resolution canvas and scaled up without smoothing, so the stars,
// the double helix and the title are made of real chunky pixels.

import { prefersReducedMotion, storageGet, storageSet } from "./util.js";

const BASE_COL = { A: [255, 62, 62], C: [60, 110, 255], G: [57, 255, 120], T: [255, 225, 60] };
const PAIR = { A: "T", T: "A", C: "G", G: "C" };
const SEQ = "GATTACACGTGCATGCCGTAACGGATCCATGGCTAGCAAGCTTGCGGCCGCTCTAGAGGATCC";

export class Landing {
  constructor(root) {
    this.root = root;
    this.canvas = root.querySelector("canvas");
    this.ctx = this.canvas.getContext("2d");
    this.low = document.createElement("canvas");
    this.lctx = this.low.getContext("2d");
    this.t0 = performance.now();
    this.mouse = [0.5, 0.5];
    this.exiting = 0;
    this.particles = [];
    this.stars = [];
    this.running = false;
    this.reduced = prefersReducedMotion();
    this.onResize = () => this.resize();
    this.onMove = (e) => (this.mouse = [e.clientX / innerWidth, e.clientY / innerHeight]);
  }

  async start() {
    this.root.hidden = false;
    this.resize();
    addEventListener("resize", this.onResize);
    addEventListener("pointermove", this.onMove);
    try {
      await document.fonts.load('16px "Press Start 2P"');
    } catch {
      /* fall back to monospace */
    }
    if (this.stopped) return; // entered before the font arrived
    this.makeTitle();
    this.running = true;
    const tick = (t) => {
      if (!this.running) return;
      this.frame(t);
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
    const sub = this.root.querySelector(".tagline");
    if (sub) typewriter(sub, sub.dataset.text, this.reduced ? 0 : 28, this.reduced ? 0 : 900);
    const skip = this.root.querySelector("#skip-intro");
    if (skip) {
      skip.checked = !!storageGet("gpxil.skipIntro", false);
      skip.addEventListener("change", () => storageSet("gpxil.skipIntro", skip.checked));
    }
  }

  resize() {
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(innerWidth * dpr);
    this.canvas.height = Math.round(innerHeight * dpr);
    this.px = innerWidth < 600 ? 3 : 4;
    this.low.width = Math.ceil(innerWidth / this.px);
    this.low.height = Math.ceil(innerHeight / this.px);
    this.stars = Array.from({ length: Math.round((this.low.width * this.low.height) / 140) }, () => this.newStar(true));
    if (this.titlePoints) this.makeTitle();
  }

  newStar(anywhere = false) {
    const keys = "ACGT";
    return { x: (Math.random() - 0.5) * 2, y: (Math.random() - 0.5) * 2, z: anywhere ? Math.random() : 1, base: keys[(Math.random() * 4) | 0] };
  }

  makeTitle() {
    // rasterise the title at low resolution and turn its pixels into particle targets
    const W = this.low.width;
    const H = this.low.height;
    const lines = W < 150 ? ["ENTER", "THE", "PIXELVERSE"] : ["ENTER THE", "PIXELVERSE"];
    const size = Math.max(8, Math.floor(Math.min(W / (lines.reduce((m, l) => Math.max(m, l.length), 0) + 2), H / 9)));
    const c = document.createElement("canvas");
    c.width = W;
    c.height = H;
    const x = c.getContext("2d");
    x.fillStyle = "#fff";
    x.textAlign = "center";
    x.textBaseline = "top";
    x.font = `${size}px "Press Start 2P", ui-monospace, monospace`;
    const top = Math.round(H * 0.16);
    lines.forEach((l, i) => x.fillText(l, W / 2, top + i * Math.round(size * 1.5)));
    const d = x.getImageData(0, 0, W, H).data;
    const pts = [];
    for (let yy = 0; yy < H; yy++) for (let xx = 0; xx < W; xx++) if (d[4 * (yy * W + xx) + 3] > 128) pts.push([xx, yy]);
    this.titlePoints = pts;
    this.titleBottom = top + lines.length * Math.round(size * 1.5);
    this.particles = pts.map(([tx, ty], i) => ({
      tx,
      ty,
      x: Math.random() * W,
      y: Math.random() * H,
      delay: this.reduced ? 0 : Math.random() * 0.9 + (i / pts.length) * 0.4,
      base: "ACGT"[(tx + ty * 3) % 4],
    }));
  }

  frame(now) {
    const t = (now - this.t0) / 1000;
    const { lctx: g, low } = this;
    const W = low.width;
    const H = low.height;
    const warp = this.exiting ? Math.min(1, (now - this.exiting) / 700) : 0;
    g.fillStyle = `rgba(7,7,15,${this.reduced ? 1 : 0.55 - 0.4 * warp})`;
    g.fillRect(0, 0, W, H);
    // starfield of bases flying past
    const speed = this.reduced ? 0 : 0.0035 + warp * 0.06;
    for (const s of this.stars) {
      s.z -= speed;
      if (s.z <= 0.02) Object.assign(s, this.newStar());
      const k = 0.5 / s.z;
      const sx = Math.round(W / 2 + s.x * k * W * 0.5 + (this.mouse[0] - 0.5) * 10 * (1 - s.z));
      const sy = Math.round(H / 2 + s.y * k * H * 0.5 + (this.mouse[1] - 0.5) * 10 * (1 - s.z));
      if (sx < 0 || sy < 0 || sx >= W || sy >= H) continue;
      const c = BASE_COL[s.base];
      const a = Math.min(1, (1 - s.z) * 1.4);
      g.fillStyle = `rgba(${c[0]},${c[1]},${c[2]},${a.toFixed(2)})`;
      const sz = s.z < 0.25 ? 2 : 1;
      g.fillRect(sx, sy, sz, sz + (warp > 0.2 ? Math.round(warp * 6) : 0));
    }
    // the double helix, turning
    this.drawHelix(g, W, H, t, warp);
    // the title, assembling from flying bases
    for (const p of this.particles) {
      const k = this.reduced ? 1 : clamp01((t - p.delay) / 1.1);
      const e = 1 - (1 - k) ** 3;
      let x = p.x + (p.tx - p.x) * e;
      let y = p.y + (p.ty - p.y) * e;
      if (warp) {
        x += (p.tx - W / 2) * warp * 3;
        y += (p.ty - H / 2) * warp * 3;
      }
      const shimmer = k >= 1 && !this.reduced ? (Math.sin(t * 3 + p.tx * 0.15) + 1) / 2 : 0;
      const c = k >= 1 ? mix([255, 255, 255], BASE_COL[p.base], 0.25 + 0.5 * shimmer) : BASE_COL[p.base];
      g.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
      g.fillRect(Math.round(x), Math.round(y), 1, 1);
    }
    // scale up, no smoothing
    const ctx = this.ctx;
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(low, 0, 0, W, H, 0, 0, this.canvas.width, this.canvas.height);
    if (warp) {
      ctx.fillStyle = `rgba(255,255,255,${Math.max(0, warp - 0.6) * 2})`;
      ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
    }
  }

  drawHelix(g, W, H, t, warp) {
    const top = (this.titleBottom || H * 0.4) + 6;
    const ui = this.root.querySelector(".landing-ui");
    const bottom = ui ? Math.min(H - 6, Math.floor(ui.getBoundingClientRect().top / this.px) - 4) : H - 6;
    if (bottom - top < 20) return;
    const vertical = W < H * 0.9;
    const cx = W / 2 + (this.mouse[0] - 0.5) * 12;
    const len = vertical ? bottom - top : W * 0.84;
    const R = vertical ? Math.min(W * 0.18, 26) : Math.min((bottom - top) * 0.38, 26);
    const steps = Math.floor(len / 3);
    const spin = this.reduced ? 0.6 : t * 1.6 + warp * 20;
    const pts = [];
    for (let i = 0; i < steps; i++) {
      const a = i * 0.32 + spin;
      const u = i * 3;
      const base = SEQ[i % SEQ.length];
      pts.push({ i, a, u, base, z: Math.sin(a) });
    }
    // far strand first, then rungs, then near strand
    const place = (u, off) => (vertical ? [Math.round(cx + off), Math.round(top + u)] : [Math.round(W * 0.08 + u), Math.round((top + bottom) / 2 + off)]);
    for (const p of pts) {
      const off = Math.cos(p.a) * R;
      // base pair rung, drawn every step
      const [x1, y1] = place(p.u, off);
      const [x2, y2] = place(p.u, -off);
      const n = Math.max(1, Math.round(Math.hypot(x2 - x1, y2 - y1)));
      const shade = 0.45 + 0.55 * ((p.z + 1) / 2);
      for (let k = 0; k <= n; k++) {
        const b = k < n / 2 ? p.base : PAIR[p.base];
        const c = BASE_COL[b];
        g.fillStyle = `rgb(${(c[0] * shade) | 0},${(c[1] * shade) | 0},${(c[2] * shade) | 0})`;
        g.fillRect(Math.round(x1 + ((x2 - x1) * k) / n), Math.round(y1 + ((y2 - y1) * k) / n), 1, 1);
      }
    }
    for (const p of pts) {
      for (const sgn of [1, -1]) {
        const off = sgn * Math.cos(p.a) * R;
        const depth = sgn * p.z;
        const [x, y] = place(p.u, off);
        const v = Math.round(150 + 105 * ((depth + 1) / 2));
        g.fillStyle = `rgb(${v},${v},${Math.min(255, v + 30)})`;
        g.fillRect(x - 1, y - 1, depth > 0 ? 3 : 2, depth > 0 ? 3 : 2);
      }
    }
  }

  /** Play the warp-out and resolve when the landing has gone. */
  exit() {
    return new Promise((resolve) => {
      if (this.reduced || !this.running) {
        this.stop();
        return resolve();
      }
      this.exiting = performance.now();
      this.root.classList.add("leaving");
      setTimeout(() => {
        this.stop();
        resolve();
      }, 820);
    });
  }

  stop() {
    this.stopped = true;
    this.running = false;
    this.root.hidden = true;
    removeEventListener("resize", this.onResize);
    removeEventListener("pointermove", this.onMove);
  }
}

function clamp01(x) {
  return Math.max(0, Math.min(1, x));
}

function mix(a, b, k) {
  return [0, 1, 2].map((i) => Math.round(a[i] + (b[i] - a[i]) * k));
}

function typewriter(el, text, ms, delay) {
  el.textContent = "";
  if (!ms) {
    el.textContent = text;
    return;
  }
  let i = 0;
  setTimeout(function step() {
    el.textContent = text.slice(0, ++i);
    if (i < text.length) setTimeout(step, ms);
  }, delay);
}
