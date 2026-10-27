// Small shared helpers: formatting, fetching (ranges, gzip, deflate), caches, events.

export const fmt = (n) => Math.round(n).toLocaleString("en-US");

/** 1234567 -> "1.23 Mb", 4200 -> "4.2 kb", 87 -> "87 bp" */
export function fmtBp(n) {
  const a = Math.abs(n);
  if (a >= 1e9) return `${(n / 1e9).toFixed(2)} Gb`;
  if (a >= 1e6) return `${(n / 1e6).toFixed(a >= 1e8 ? 0 : a >= 1e7 ? 1 : 2)} Mb`;
  if (a >= 1e3) return `${(n / 1e3).toFixed(a >= 1e5 ? 0 : 1)} kb`;
  return `${Math.round(n)} bp`;
}

export const clamp = (x, lo, hi) => Math.max(lo, Math.min(hi, x));

export function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
}

/** Least-recently-used cache with a size budget. */
export class LRU {
  constructor(budget, sizeOf = () => 1) {
    this.budget = budget;
    this.sizeOf = sizeOf;
    this.map = new Map();
    this.size = 0;
  }
  get(k) {
    if (!this.map.has(k)) return undefined;
    const v = this.map.get(k);
    this.map.delete(k);
    this.map.set(k, v);
    return v;
  }
  has(k) {
    return this.map.has(k);
  }
  set(k, v) {
    if (this.map.has(k)) {
      this.size -= this.sizeOf(this.map.get(k));
      this.map.delete(k);
    }
    this.map.set(k, v);
    this.size += this.sizeOf(v);
    while (this.size > this.budget && this.map.size > 1) {
      const [oldK, oldV] = this.map.entries().next().value;
      this.map.delete(oldK);
      this.size -= this.sizeOf(oldV);
    }
  }
}

/** Tiny event emitter. */
export class Emitter {
  constructor() {
    this.handlers = new Map();
  }
  on(name, fn) {
    if (!this.handlers.has(name)) this.handlers.set(name, new Set());
    this.handlers.get(name).add(fn);
    return () => this.handlers.get(name).delete(fn);
  }
  emit(name, ...args) {
    for (const fn of this.handlers.get(name) || []) fn(...args);
  }
}

/** Activity tracking so the DNA loader can show while anything is in flight. */
export const activity = new Emitter();
let inFlight = 0;
const pending = new Map(); // label -> number in flight

/** Labels of everything loading right now, e.g. ["DNA", "Ensembl variants"]. */
export function busyWith() {
  return [...pending.keys()].filter(Boolean);
}

/** Await a promise, counting it as "loading" (the DNA loader shows when this takes a while). */
export async function track(promise, label = "") {
  inFlight++;
  pending.set(label, (pending.get(label) || 0) + 1);
  activity.emit("change", inFlight, label);
  try {
    return await promise;
  } finally {
    inFlight--;
    const n = pending.get(label) - 1;
    if (n > 0) pending.set(label, n);
    else pending.delete(label);
    activity.emit("change", inFlight, label);
  }
}

async function streamToBytes(stream) {
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** zlib (TIFF Deflate) -> bytes, with the browser's own decompressor. */
export async function inflate(bytes) {
  return streamToBytes(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("deflate")));
}

export async function gunzip(bytes) {
  return streamToBytes(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip")));
}

// Each fetch is tracked until its body has arrived, not just its headers.

export function fetchBytes(url, label) {
  return track(
    (async () => {
      const r = await fetch(url);
      if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
      return new Uint8Array(await r.arrayBuffer());
    })(),
    label
  );
}

/** Fetch and parse a .json.gz (decompressed here, so it works on any static host). */
export async function fetchJsonGz(url, label) {
  const bytes = await fetchBytes(url, label);
  const plain = bytes[0] === 0x1f && bytes[1] === 0x8b ? await gunzip(bytes) : bytes;
  return JSON.parse(new TextDecoder().decode(plain));
}

export function fetchJson(url, label, options) {
  return track(
    (async () => {
      const r = await fetch(url, options);
      if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
      return r.json();
    })(),
    label
  );
}

/** Bytes [start, end) of a file, by an HTTP range request. */
export function fetchRange(url, start, end, label) {
  return track(
    (async () => {
      const r = await fetch(url, { headers: { Range: `bytes=${start}-${end - 1}` } });
      if (r.status === 206) return new Uint8Array(await r.arrayBuffer());
      if (r.status === 200) {
        // the server ignored the range: still correct, but slow for big files
        console.warn(`${url}: server does not support range requests`);
        const all = new Uint8Array(await r.arrayBuffer());
        return all.slice(start, end);
      }
      throw new Error(`${url}: HTTP ${r.status}`);
    })(),
    label
  );
}

export function debounce(fn, ms) {
  let t = null;
  return (...args) => {
    clearTimeout(t);
    t = setTimeout(() => fn(...args), ms);
  };
}

export const prefersReducedMotion = () => window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

export function storageGet(key, fallback) {
  try {
    const v = localStorage.getItem(key);
    return v === null ? fallback : JSON.parse(v);
  } catch {
    return fallback;
  }
}

export function storageSet(key, value) {
  try {
    localStorage.setItem(key, JSON.stringify(value));
  } catch {
    /* private mode or storage disabled: settings just aren't remembered */
  }
}

/** Upper bound: first index i with arr[i] > x (arr sorted ascending). */
export function upperBound(arr, x) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (arr[mid] <= x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

export function lowerBound(arr, x) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (arr[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}
