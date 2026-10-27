// pxil files: one sequence as a picture, one pixel per residue. A pxil file
// is a PNG whose palette indices are the residues, with the alphabet and other
// details in an iTXt chunk called "pxil" (https://github.com/comparativechrono/pxil).
// The browser writes them for genes, transcripts and proteins, and zips them
// when several are wanted at once.

export const PXIL_VERSION = "0.1";
export const DNA = "-ACMGRSVTWYHKDBN-acmgrsvtwyhkdbn"; // the gpxil reference codes; +16 = soft-masked
export const PROTEIN = "-ACDEFGHIKLMNPQRSTVWYUOBZJX*";

const CRC = new Uint32Array(256);
for (let n = 0; n < 256; n++) {
  let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  CRC[n] = c >>> 0;
}

export function crc32(bytes, crc = 0) {
  let c = ~crc >>> 0;
  for (let i = 0; i < bytes.length; i++) c = CRC[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return ~c >>> 0;
}

const ascii = (s) => Uint8Array.from(s, (ch) => ch.charCodeAt(0) & 0xff);

function concat(parts) {
  const out = new Uint8Array(parts.reduce((n, p) => n + p.length, 0));
  let o = 0;
  for (const p of parts) {
    out.set(p, o);
    o += p.length;
  }
  return out;
}

async function zlib(bytes) {
  return new Uint8Array(await new Response(new Blob([bytes]).stream().pipeThrough(new CompressionStream("deflate"))).arrayBuffer());
}

function chunk(type, data) {
  const out = new Uint8Array(12 + data.length);
  const v = new DataView(out.buffer);
  v.setUint32(0, data.length);
  out.set(ascii(type), 4);
  out.set(data, 8);
  v.setUint32(8 + data.length, crc32(out.subarray(4, 8 + data.length)));
  return out;
}

/** Indices of a sequence's characters in an alphabet (the first index of each character). */
export function indicesOf(seq, alphabet) {
  const at = new Map();
  [...alphabet].forEach((ch, i) => at.has(ch) || at.set(ch, i));
  const out = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) {
    const k = at.get(seq[i]);
    if (k === undefined) throw new Error(`“${seq[i]}” is not in the alphabet ${alphabet}`);
    out[i] = k;
  }
  return out;
}

export async function sha256Hex(text) {
  const d = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text));
  return [...new Uint8Array(d)].map((b) => b.toString(16).padStart(2, "0")).join("");
}

/**
 * A pxil file (the bytes of a PNG): `indices` into `alphabet`, `width`
 * residues to a row, `colours` ([r, g, b] for each index) for image viewers,
 * and `meta` added to the required metadata.
 */
export async function pxilFile({ indices, alphabet, width, colours, meta = {}, sequence = null }) {
  const n = indices.length;
  const w = Math.max(1, Math.round(width));
  const h = Math.max(1, Math.ceil(n / w));
  const ihdr = new Uint8Array(13);
  const v = new DataView(ihdr.buffer);
  v.setUint32(0, w);
  v.setUint32(4, h);
  ihdr[8] = 8; // bit depth
  ihdr[9] = 3; // indexed colour
  const plte = new Uint8Array(3 * [...alphabet].length);
  for (let i = 0; i < plte.length / 3; i++) plte.set((colours[i] || [0, 0, 0]).map((x) => Math.round(x) & 0xff), 3 * i);
  const info = { pxil: PXIL_VERSION, ...meta, alphabet, length: n, width: w };
  if (sequence !== null && !info.sha256) info.sha256 = await sha256Hex(sequence);
  // iTXt: keyword, 0, compression flag 0, method 0, language tag "", 0, translated keyword "", 0, UTF-8 text
  const itxt = concat([ascii("pxil"), new Uint8Array([0, 0, 0, 0, 0]), new TextEncoder().encode(JSON.stringify(info))]);
  const title = [info.name, info.id].filter(Boolean).join(" ");
  const raw = new Uint8Array(h * (w + 1)); // each row: filter type 0, then its indices
  for (let r = 0; r < h; r++) raw.set(indices.subarray(r * w, Math.min(n, (r + 1) * w)), r * (w + 1) + 1);
  return concat([
    new Uint8Array([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk("IHDR", ihdr),
    chunk("PLTE", plte),
    ...(title && /^[\x20-\x7e]*$/.test(title) ? [chunk("tEXt", concat([ascii("Title"), new Uint8Array([0]), ascii(title)]))] : []),
    chunk("iTXt", itxt),
    chunk("IDAT", await zlib(raw)),
    chunk("IEND", new Uint8Array(0)),
  ]);
}

/** A zip archive of files ({name, bytes}), stored without compression (PNG data is compressed already). */
export function zip(files) {
  const now = new Date();
  const time = (now.getHours() << 11) | (now.getMinutes() << 5) | (now.getSeconds() >> 1);
  const date = ((now.getFullYear() - 1980) << 9) | ((now.getMonth() + 1) << 5) | now.getDate();
  const locals = [];
  const central = [];
  let offset = 0;
  for (const f of files) {
    const name = new TextEncoder().encode(f.name);
    const crc = crc32(f.bytes);
    const head = new Uint8Array(30 + name.length);
    const lv = new DataView(head.buffer);
    lv.setUint32(0, 0x04034b50, true);
    lv.setUint16(4, 20, true); // version needed
    lv.setUint16(6, 0x0800, true); // UTF-8 names
    lv.setUint16(8, 0, true); // stored
    lv.setUint16(10, time, true);
    lv.setUint16(12, date, true);
    lv.setUint32(14, crc, true);
    lv.setUint32(18, f.bytes.length, true);
    lv.setUint32(22, f.bytes.length, true);
    lv.setUint16(26, name.length, true);
    head.set(name, 30);
    const cd = new Uint8Array(46 + name.length);
    const cv = new DataView(cd.buffer);
    cv.setUint32(0, 0x02014b50, true);
    cv.setUint16(4, 20, true); // made by
    cv.setUint16(6, 20, true); // needed
    cv.setUint16(8, 0x0800, true);
    cv.setUint16(10, 0, true);
    cv.setUint16(12, time, true);
    cv.setUint16(14, date, true);
    cv.setUint32(16, crc, true);
    cv.setUint32(20, f.bytes.length, true);
    cv.setUint32(24, f.bytes.length, true);
    cv.setUint16(28, name.length, true);
    cv.setUint32(42, offset, true);
    cd.set(name, 46);
    locals.push(head, f.bytes);
    central.push(cd);
    offset += head.length + f.bytes.length;
  }
  const cdSize = central.reduce((n, c) => n + c.length, 0);
  const end = new Uint8Array(22);
  const ev = new DataView(end.buffer);
  ev.setUint32(0, 0x06054b50, true);
  ev.setUint16(8, files.length, true);
  ev.setUint16(10, files.length, true);
  ev.setUint32(12, cdSize, true);
  ev.setUint32(16, offset, true);
  return concat([...locals, ...central, end]);
}

/** A name safe for files: letters, digits, dot, dash and underscore. */
export const safeName = (s) => String(s).replace(/[^A-Za-z0-9._-]+/g, "_");
