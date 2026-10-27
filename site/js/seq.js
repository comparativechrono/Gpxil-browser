// Sequence helpers: codes to letters, reverse complement, splicing, translation.

import { IUPAC } from "./palettes.js";

export function lettersOf(codes, keepCase = false) {
  let s = "";
  for (let i = 0; i < codes.length; i++) {
    const c = codes[i];
    const ch = IUPAC[c & 15];
    s += keepCase && c & 0x10 ? ch.toLowerCase() : ch;
  }
  return s;
}

const COMP = { A: "T", C: "G", G: "C", T: "A", N: "N", R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K", B: "V", V: "B", D: "H", H: "D", "-": "-" };
export function revcomp(s) {
  let out = "";
  for (let i = s.length - 1; i >= 0; i--) {
    const ch = s[i];
    const up = ch.toUpperCase();
    const c = COMP[up] || "N";
    out += ch === up ? c : c.toLowerCase();
  }
  return out;
}

/** mRNA (5'->3') of a transcript from the codes of [regionStart, ...). Returns {seq, map} with genomic positions. */
export function spliced(codes, regionStart, tx) {
  let seq = "";
  const map = [];
  for (const [a, b] of tx.exons) {
    const part = lettersOf(codes.subarray(a - regionStart, b - regionStart));
    seq += part;
    for (let p = a; p < b; p++) map.push(p);
  }
  if (tx.strand < 0) {
    seq = revcomp(seq);
    map.reverse();
  }
  return { seq, map };
}

const BASES = "TCAG";
const STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const VERT_MITO = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";

export function translate(cds, mito = false) {
  const table = mito ? VERT_MITO : STANDARD;
  let aa = "";
  for (let i = 0; i + 3 <= cds.length; i += 3) {
    const a = BASES.indexOf(cds[i].toUpperCase());
    const b = BASES.indexOf(cds[i + 1].toUpperCase());
    const c = BASES.indexOf(cds[i + 2].toUpperCase());
    aa += a < 0 || b < 0 || c < 0 ? "X" : table[a * 16 + b * 4 + c];
  }
  return aa;
}

/**
 * Selenoproteins read some TGA codons as selenocysteine (U). `sec` holds the
 * genomic position of the first base of each such codon; `map` and `start` are
 * the spliced mRNA's positions and where translation starts in it.
 */
export function markSelenocysteine(aa, map, start, sec) {
  if (!sec || !sec.length) return aa;
  const at = new Set(sec);
  const out = aa.split("");
  for (let i = 0; i < out.length; i++) if (at.has(map[start + 3 * i])) out[i] = "U";
  return out.join("");
}

/** CDS of a transcript as an mRNA slice: [from, to) indexes into the spliced sequence. */
export function cdsRange(tx, map) {
  if (!tx.cds) return null;
  const [cs, ce] = tx.cds;
  let from = -1;
  let to = -1;
  for (let i = 0; i < map.length; i++) {
    const p = map[i];
    if (p >= cs && p < ce) {
      if (from < 0) from = i;
      to = i + 1;
    }
  }
  return from < 0 ? null : [from, to];
}

export function gcFraction(codes) {
  let gc = 0;
  let acgt = 0;
  for (let i = 0; i < codes.length; i++) {
    const b = codes[i] & 15;
    if (b === 1 || b === 2 || b === 4 || b === 8) {
      acgt++;
      if (b === 2 || b === 4) gc++;
    }
  }
  return acgt ? gc / acgt : 0;
}

export function wrapFasta(name, seq, width = 60) {
  let out = `>${name}\n`;
  for (let i = 0; i < seq.length; i += width) out += `${seq.slice(i, i + width)}\n`;
  return out;
}
