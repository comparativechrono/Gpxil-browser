// The DNA loader: a little pixel double helix that spins while data is on its
// way. It appears only when something takes a while (over 350 ms).

import { activity, busyWith, prefersReducedMotion } from "./util.js";

const COL = { A: "#ff3e3e", C: "#3c6eff", G: "#39ff78", T: "#ffe13c" };
const PAIR = { A: "T", T: "A", C: "G", G: "C" };
const SEQ = "ACGTTGCAGCTA";
// what the loader says depends on what it is waiting for
const LINES = {
  local: ["UNZIPPING DNA", "FETCHING BASES", "COUNTING Gs AND Cs", "READING STRIPS", "PIXELATING"],
  ensembl: ["ASKING ENSEMBL", "PHONING HINXTON", "CONSULTING ENSEMBL"],
  mygene: ["ASKING MYGENE.INFO", "READING THE NCBI NOTES"],
};

function themeOf(labels) {
  if (labels.some((l) => l.startsWith("Ensembl"))) return "ensembl";
  if (labels.some((l) => l.startsWith("NCBI"))) return "mygene";
  return "local";
}

export class DnaLoader {
  constructor(root) {
    this.root = root;
    this.canvas = root.querySelector("canvas");
    this.text = root.querySelector(".loader-text");
    this.dots = root.querySelector(".loader-dots");
    this.what = root.querySelector(".loader-what");
    this.ctx = this.canvas.getContext("2d");
    this.canvas.width = 20;
    this.canvas.height = 40;
    this.visible = false;
    this.timer = null;
    this.raf = 0;
    this.count = 0;
    activity.on("change", (n) => this.set(n));
  }

  set(n) {
    this.count = n;
    if (n > 0 && !this.visible && !this.timer) this.timer = setTimeout(() => this.show(), 350);
    if (n === 0) {
      clearTimeout(this.timer);
      this.timer = null;
      if (this.visible) this.hide();
    }
  }

  show() {
    this.timer = null;
    if (this.count === 0) return;
    this.visible = true;
    this.root.hidden = false;
    this.root.classList.add("on");
    this.t0 = performance.now();
    this.theme = null;
    const tick = (t) => {
      if (!this.visible) return;
      this.draw(t);
      this.raf = requestAnimationFrame(tick);
    };
    this.raf = requestAnimationFrame(tick);
  }

  hide() {
    this.visible = false;
    cancelAnimationFrame(this.raf);
    this.root.classList.remove("on");
    this.root.hidden = true;
  }

  draw(now) {
    const t = (now - this.t0) / 1000;
    const g = this.ctx;
    g.fillStyle = "#07070f";
    g.fillRect(0, 0, 20, 40);
    const spin = prefersReducedMotion() ? 0 : t * 5;
    for (let y = 1; y < 39; y += 2) {
      const i = (y - 1) / 2;
      const a = i * 0.55 + spin;
      const off = Math.cos(a) * 7;
      const x1 = Math.round(10 + off);
      const x2 = Math.round(10 - off);
      const b = SEQ[i % SEQ.length];
      const lo = Math.min(x1, x2);
      const hi = Math.max(x1, x2);
      for (let x = lo + 1; x < hi; x++) {
        g.fillStyle = COL[x < 10 === x1 < x2 ? b : PAIR[b]];
        g.fillRect(x, y, 1, 1);
      }
      const near = Math.sin(a) > 0;
      g.fillStyle = near ? "#ffffff" : "#8a8ab0";
      g.fillRect(x1, y, 1, 1);
      g.fillStyle = near ? "#8a8ab0" : "#ffffff";
      g.fillRect(x2, y, 1, 1);
    }
    const labels = busyWith();
    const theme = themeOf(labels);
    if (theme !== this.theme) {
      this.theme = theme;
      this.line = LINES[theme][(Math.random() * LINES[theme].length) | 0];
    }
    // the words change rarely (screen readers announce them); the dots, hidden from them, every frame
    if (this.text.textContent !== this.line) this.text.textContent = this.line;
    if (this.dots) this.dots.textContent = ".".repeat(1 + (Math.floor(t * 3) % 3)).padEnd(3, " ");
    if (this.what) {
      let what = [...new Set(labels.map((l) => l.replace(/^gpxil /, "")))].slice(0, 3).join(" · ");
      if (theme !== "local" && t > 4) what += theme === "ensembl" ? " · Ensembl can take a few seconds" : " · still waiting";
      if (this.what.textContent !== what) this.what.textContent = what;
    }
  }
}
