// Palettes turn gpxil codes into colours at display time (spec section 12).
// The four gpxil palettes match gpxiltools exactly; "neon" and "gc" are
// extra views for this browser.

export const PALETTES = {
  classic: {
    label: "Classic",
    about: "A red, C blue, G green, T yellow: the colours of the original gpxil notebooks.",
    rgb: { A: [255, 0, 0], C: [0, 0, 255], G: [0, 255, 0], T: [255, 255, 0], N: [128, 128, 128], gap: [255, 255, 255] },
    mix: "mean",
    fade: 0.45,
    fadeTo: [255, 255, 255],
  },
  gpxi: {
    label: "GPxI greys",
    about: "Albrecht-Buehler's 2009 genome pixel images: A black, C dark grey, T light grey, G white.",
    rgb: { A: [0, 0, 0], C: [85, 85, 85], G: [255, 255, 255], T: [170, 170, 170], N: [170, 200, 235], gap: [255, 250, 230] },
    mix: "mean",
    fade: 0,
    fadeTo: [255, 255, 255],
  },
  "okabe-ito": {
    label: "Okabe-Ito",
    about: "Colour-blind safe, in the IGV order: A green, C blue, G orange, T vermillion.",
    rgb: { A: [0, 158, 115], C: [0, 114, 178], G: [230, 159, 0], T: [213, 94, 0], N: [153, 153, 153], gap: [255, 255, 255] },
    mix: "mean",
    fade: 0.45,
    fadeTo: [255, 255, 255],
  },
  process: {
    label: "Process inks",
    about: "One printing ink per base: A cyan, C magenta, G yellow, T black.",
    rgb: { A: [0, 174, 239], C: [236, 0, 140], G: [255, 242, 0], T: [35, 31, 32], N: [160, 160, 160], gap: [255, 255, 255] },
    mix: "overprint",
    fade: 0.45,
    fadeTo: [255, 255, 255],
  },
  neon: {
    label: "Neon",
    about: "The pixelverse: bright bases on a dark sky. Repeats sink into the dark.",
    rgb: { A: [255, 62, 165], C: [32, 227, 255], G: [57, 255, 136], T: [255, 225, 77], N: [46, 46, 72], gap: [7, 7, 15] },
    mix: "mean",
    fade: 0.55,
    fadeTo: [10, 10, 24],
  },
  gc: {
    label: "GC heat",
    about: "G and C hot, A and T cold: isochores and CpG islands stand out.",
    rgb: { A: [36, 84, 196], C: [255, 122, 36], G: [255, 122, 36], T: [36, 84, 196], N: [90, 90, 90], gap: [255, 255, 255] },
    mix: "mean",
    fade: 0.3,
    fadeTo: [255, 255, 255],
  },
};

export const BASES = ["A", "C", "G", "T"];
export const SOFTMASK = 0x10;
export const IUPAC = "-ACMGRSVTWYHKDBN"; // index = 4-bit code

function mix(cols, how) {
  if (how === "overprint") return [0, 1, 2].map((i) => Math.round(cols.reduce((p, c) => p * (c[i] / 255), 1) * 255));
  return [0, 1, 2].map((i) => Math.round(cols.reduce((s, c) => s + c[i], 0) / cols.length));
}

function fade(c, t, to) {
  return c.map((v, i) => Math.round(v + (to[i] - v) * t));
}

/** True for the name of a palette (and not, say, "toString", which every object has). */
export const hasPalette = (name) => typeof name === "string" && Object.hasOwn(PALETTES, name);

/** RGB for every ref code (256 x 3) and a packed RGBA lookup for ImageData. */
export function buildLut(name) {
  const p = hasPalette(name) ? PALETTES[name] : PALETTES.classic;
  const rgb = new Uint8Array(256 * 3);
  for (let code = 0; code < 256; code++) {
    const low = code & 15;
    let c;
    if (code > 31 || code === SOFTMASK) c = p.rgb.gap;
    else if (low === 0) c = p.rgb.gap;
    else if (low === 15) c = p.rgb.N;
    else c = mix(BASES.filter((b, i) => low & (1 << i)).map((b) => p.rgb[b]), p.mix);
    if (code & SOFTMASK && low !== 0) c = fade(c, p.fade, p.fadeTo);
    rgb.set(c, code * 3);
  }
  const packed = new Uint32Array(256);
  const little = new Uint8Array(new Uint32Array([1]).buffer)[0] === 1;
  for (let i = 0; i < 256; i++) {
    const [r, g, b] = [rgb[3 * i], rgb[3 * i + 1], rgb[3 * i + 2]];
    packed[i] = little ? (255 << 24) | (b << 16) | (g << 8) | r : (r << 24) | (g << 16) | (b << 8) | 255;
  }
  const lut = { name, palette: p, rgb, packed, little };
  lut.mean = binColour(lut, 0.41, 0, 0.5); // the genome's average A/C/G/T colour, for contrast stretching
  return lut;
}

/** Stretch a colour away from the genome's average A/C/G/T colour (zoomed-out pixels average towards it). */
export function boost(lut, c, k) {
  if (k === 1) return c;
  const m = lut.mean;
  for (let i = 0; i < 3; i++) c[i] = Math.max(0, Math.min(255, m[i] + (c[i] - m[i]) * k));
  return c;
}

/**
 * Colour of a summary bin from its GC fraction (of A/C/G/T), N fraction and
 * soft-masked fraction (0..1 each). With k > 1 the A/C/G/T part is stretched
 * away from the genome's average colour, so isochores stand out when many
 * bases are averaged; N stays the palette's N colour.
 */
export function binColour(lut, gc, n, masked, out = [0, 0, 0], k = 1) {
  const p = lut.palette;
  for (let i = 0; i < 3; i++) {
    let v = 0.5 * (1 - gc) * (p.rgb.A[i] + p.rgb.T[i]) + 0.5 * gc * (p.rgb.C[i] + p.rgb.G[i]);
    v += (p.fadeTo[i] - v) * p.fade * masked;
    out[i] = v;
  }
  if (k !== 1 && lut.mean) boost(lut, out, k);
  for (let i = 0; i < 3; i++) out[i] = out[i] * (1 - n) + p.rgb.N[i] * n;
  return out;
}

/** Packs r, g, b into the Uint32 layout of ImageData. */
export function pack(lut, r, g, b) {
  return lut.little ? (255 << 24) | (b << 16) | (g << 8) | r : (r << 24) | (g << 16) | (b << 8) | 255;
}

/** Colours for amino acids (Taylor 1997 scheme), "*" for stop. */
export const AMINO = {
  A: [204, 255, 0], R: [0, 0, 255], N: [204, 0, 255], D: [255, 0, 0], C: [255, 255, 0],
  Q: [255, 0, 204], E: [255, 0, 102], G: [255, 153, 0], H: [0, 102, 255], I: [102, 255, 0],
  L: [51, 255, 0], K: [102, 0, 255], M: [0, 255, 0], F: [0, 255, 102], P: [255, 204, 0],
  S: [255, 51, 0], T: [255, 102, 0], W: [0, 204, 255], Y: [0, 255, 204], V: [153, 255, 0],
  U: [255, 120, 200], "*": [40, 40, 48], X: [110, 110, 120],
};

/** Readable text colour (black or white) on a background colour. */
export function inkOn([r, g, b]) {
  return 0.299 * r + 0.587 * g + 0.114 * b > 140 ? "#000" : "#fff";
}
