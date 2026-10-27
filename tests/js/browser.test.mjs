// The browser's own modules, run in Node against the toy genome built by
// tests/build_fixture.py and served with range requests, as GitHub Pages does.
//
//   node --test tests/js/

import { test, before, after } from "node:test";
import assert from "node:assert/strict";
import http from "node:http";
import fs from "node:fs";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";

import { Genome } from "../../site/js/gpxil.js";
import { Overview, decodeBin } from "../../site/js/overview.js";
import { GeneDB } from "../../site/js/genes.js";
import { parsePlace } from "../../site/js/search.js";
import { cdsRange, gcFraction, lettersOf, markSelenocysteine, revcomp, spliced, translate, wrapFasta } from "../../site/js/seq.js";
import { PALETTES, buildLut, binColour, hasPalette, inkOn } from "../../site/js/palettes.js";
import { LIMITS, LiveTracks } from "../../site/js/live.js";
import { DNA, PROTEIN, crc32, indicesOf, pxilFile, safeName, zip } from "../../site/js/pxil.js";
import { LRU, busyWith, fetchJsonGz, fmtBp, lowerBound, track, upperBound } from "../../site/js/util.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..", "..");
const SITE = path.join(ROOT, "tests", ".fixture", "site");
const INPUTS = path.join(ROOT, "tests", ".fixture", "inputs");

if (!fs.existsSync(path.join(SITE, "data", "grch38", "manifest.json"))) {
  execFileSync(process.env.PYTHON || "python3", [path.join(ROOT, "tests", "build_fixture.py")], { stdio: "inherit" });
}

const expected = JSON.parse(fs.readFileSync(path.join(INPUTS, "expected.json"), "utf8"));
const fasta = (() => {
  const out = {};
  const names = { 1: "chr1", 2: "chr2", MT: "chrM" };
  let name = null;
  for (const line of fs.readFileSync(path.join(INPUTS, "toy.fa"), "utf8").split("\n")) {
    if (line.startsWith(">")) {
      const raw = line.slice(1).split(/\s/)[0];
      name = names[raw] || raw;
      out[name] = [];
    } else if (line) out[name].push(line);
  }
  for (const k of Object.keys(out)) out[k] = out[k].join("");
  return out;
})();

// a static file server that answers range requests, like GitHub Pages
let server;
let base;
let requests = 0;
before(async () => {
  server = http.createServer((req, res) => {
    requests++;
    const file = path.join(SITE, decodeURIComponent(new URL(req.url, "http://x").pathname));
    if (!file.startsWith(SITE) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) {
      res.writeHead(404);
      return res.end();
    }
    const size = fs.statSync(file).size;
    const m = /^bytes=(\d+)-(\d*)$/.exec(req.headers.range || "");
    if (m) {
      const a = +m[1];
      const b = Math.min(m[2] ? +m[2] : size - 1, size - 1);
      res.writeHead(206, { "Content-Range": `bytes ${a}-${b}/${size}`, "Content-Length": b - a + 1 });
      fs.createReadStream(file, { start: a, end: b }).pipe(res);
    } else {
      res.writeHead(200, { "Content-Length": size });
      fs.createReadStream(file).pipe(res);
    }
  });
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  base = new URL(`http://127.0.0.1:${server.address().port}/data/grch38/`);
});
after(() => server.close());

async function load() {
  const manifest = await (await fetch(new URL("manifest.json", base))).json();
  const cytobands = await fetchJsonGz(new URL(manifest.cytobands, base).href);
  return { manifest, cytobands, genome: new Genome(manifest, base) };
}

// a small deterministic random number generator
function rng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}

test("bases read over range requests match the FASTA", async () => {
  const { genome } = await load();
  const cases = [
    ["chr1", 0, 100],
    ["chr1", 8180, 8200], // across a row of the image
    ["chr1", 262_100, 262_200], // across a strip (32 rows of 8,192)
    ["chr1", 49_990, 50_010], // into the soft-masked repeat
    ["chr1", 199_990, 200_010], // into the run of Ns
    ["chr1", 0, 400_000], // the whole chromosome
    ["chr2", 119_000, 120_000],
    ["chrM", 0, 16_569],
    ["KI270728.1", 990, 1_510], // scaffolds share one file
    ["GL000009.2", 0, 3_000],
  ];
  for (const [c, a, b] of cases) assert.equal(lettersOf(await genome.bases(c, a, b), true), fasta[c].slice(a, b), `${c}:${a}-${b}`);
  const r = rng(7);
  for (let i = 0; i < 60; i++) {
    const a = Math.floor(r() * 399_000);
    const b = a + 1 + Math.floor(r() * 20_000);
    assert.equal(lettersOf(await genome.bases("chr1", a, b), true), fasta.chr1.slice(a, Math.min(b, 400_000)));
  }
  // intervals off either end are clamped
  assert.equal(lettersOf(await genome.bases("chr1", -5, 10), true), fasta.chr1.slice(0, 10));
  assert.equal((await genome.bases("chr1", 399_990, 400_100)).length, 10);
});

test("strips are cached: reading the same bases again makes no requests", async () => {
  const { genome } = await load();
  await genome.bases("chr2", 0, 5000);
  const before = requests;
  await genome.bases("chr2", 100, 4000);
  assert.equal(requests, before);
});

test("contig names and aliases resolve", async () => {
  const { genome } = await load();
  for (const [q, want] of [
    ["chr1", "chr1"],
    ["1", "chr1"],
    ["CHR2", "chr2"],
    ["MT", "chrM"],
    ["M", "chrM"],
    ["chrMT", "chrM"],
    ["KI270728.1", "KI270728.1"],
  ])
    assert.equal(genome.contig(q)?.name, want, q);
  assert.equal(genome.contig("chr3"), undefined);
});

// round half to even, as numpy's rint does in the pipeline
function rint(x) {
  const f = Math.floor(x);
  const d = x - f;
  return d > 0.5 ? f + 1 : d < 0.5 ? f : f % 2 === 0 ? f : f + 1;
}

test("overview bins match the bases they summarise", async () => {
  const { manifest } = await load();
  const ov = new Overview(manifest, base);
  for (const c of manifest.sequence.contigs) {
    const seq = fasta[c.name];
    for (let li = 0; li < ov.bins.length; li++) {
      const size = ov.bins[li];
      const bytes = await ov.level(c, li);
      assert.equal(bytes.length, 2 * Math.ceil(c.length / size));
      for (let i = 0; i < bytes.length / 2; i += li === 0 ? 37 : 1) {
        const part = seq.slice(i * size, (i + 1) * size);
        const up = part.toUpperCase();
        let acgt = 0;
        let gc = 0;
        let mask = 0;
        for (let k = 0; k < part.length; k++) {
          if ("ACGT".includes(up[k])) acgt++;
          if (up[k] === "G" || up[k] === "C") gc++;
          if (part[k] !== up[k]) mask++;
        }
        const [g, n, m] = decodeBin(bytes, i);
        assert.equal(Math.round(g * 255), acgt ? rint((255 * gc) / acgt) : 0, `${c.name} level ${li} bin ${i}`);
        assert.equal(Math.round(n * 15), rint((15 * (part.length - acgt)) / part.length));
        assert.equal(Math.round(m * 15), rint((15 * mask) / part.length));
      }
    }
  }
  assert.equal(ov.levelFor(1), 0);
  assert.equal(ov.levelFor(600), 1);
  assert.equal(ov.levelFor(1e9), ov.bins.length - 1);
});

test("genes: overlap queries, lookups and search", async () => {
  const { manifest } = await load();
  const db = new GeneDB(manifest, base);
  await db.loadContig("chr1");
  assert.deepEqual(
    db.overlapping("chr1", 100_500, 100_600).map((g) => g.name),
    ["PIXL1"]
  );
  assert.deepEqual(
    db.at("chr1", 100_010).map((g) => g.name),
    ["PIXL1-AS1", "PIXL1"] // shortest first
  );
  assert.deepEqual(db.overlapping("chr1", 0, 1000), []);
  assert.deepEqual(db.overlapping("chr2", 0, 120_000), []); // chr2 has no genes: loads as empty
  const p2 = await db.gene("ENSG00000000002", "chr1");
  assert.equal(p2.name, "PIXL2");
  assert.equal(p2.tx.strand, -1);
  assert.equal(p2.tx.canonical, true);
  assert.equal(p2.tx.mane, false);
  const p1 = await db.gene("ENSG00000000001"); // found through the index
  assert.equal(p1.tx.mane, true);
  assert.deepEqual(p1.xref, { hgnc: "90001", entrez: "999001", uniprot: "Q00001", omim: "600001" });
  assert.equal(p1.description, "pixel protein 1");

  assert.equal((await db.resolve("pixl2")).row[0], "PIXL2");
  assert.equal((await db.resolve("ENSG00000000003")).row[0], "MT-PIX");
  const viaAlias = await db.resolve("PIX1");
  assert.equal(viaAlias.row[0], "PIXL1");
  assert.equal(viaAlias.via, "PIX1");
  assert.equal((await db.resolve("OLDPIX")).row[0], "PIXL1");

  // every transcript comes from the site's own data, with one range request
  const before = requests;
  const t1 = await db.transcripts(p1);
  assert.equal(requests, before + 1);
  assert.deepEqual(
    t1.txs.map((t) => [t.id, t.canonical, t.mane, !!t.cds, t.protein]),
    [
      ["ENST00000000011", true, true, true, "ENSP00000000011"],
      ["ENST00000000012", false, false, false, null],
    ]
  );
  assert.equal(p1.tx.protein, "ENSP00000000011");
  assert.deepEqual(t1.txs[0].exons, p1.tx.exons);
  assert.deepEqual(t1.txs[0].cds, p1.tx.cds);
  assert.deepEqual(t1.txs[0].sec, p1.tx.sec);
  assert.deepEqual(t1.txs[1].exons, [
    [100_000, 100_100],
    [110_000, 110_300],
  ]);
  assert.equal(p1.txCount, 2);

  const sug = await db.suggest("pix", 10);
  assert.deepEqual(
    sug.map((s) => s.name),
    ["PIXL1", "PIXL2", "PIXL1-AS1"] // protein-coding first, then shorter names
  );
  assert.deepEqual(
    (await db.suggest("PXL", 10)).map((s) => s.name).sort(),
    ["PIXL1", "PIXL2"] // an alias shared by two genes
  );
});

test("spliced and translated proteins match the ones the toy genes encode", async () => {
  const { manifest, genome } = await load();
  const db = new GeneDB(manifest, base);
  for (const [id, e] of Object.entries(expected.genes)) {
    const g = await db.gene(id, e.contig);
    assert.ok(g, id);
    assert.deepEqual(g.tx.exons, e.exons);
    if (!e.protein) {
      assert.equal(g.tx.cds, null);
      continue;
    }
    const codes = await genome.bases(g.contig, g.tx.start, g.tx.end);
    const { seq, map } = spliced(codes, g.tx.start, g.tx);
    assert.equal(seq.length, e.exons.reduce((n, [a, b]) => n + b - a, 0));
    const range = cdsRange(g.tx, map);
    const start = range[0] + g.tx.frame;
    const protein = markSelenocysteine(translate(seq.slice(start, range[1]), g.contig === "chrM"), map, start, g.tx.sec);
    assert.equal(protein, e.protein, `${e.name} (${e.strand > 0 ? "+" : "-"} strand)`);
    assert.deepEqual(g.tx.sec, e.sec);
    // map[] gives the genomic position of every mRNA base, 5' to 3'
    assert.equal(map[range[0]], e.strand > 0 ? e.cds[0] : e.cds[1] - 1);
  }
});

test("places typed in the search box", async () => {
  const { genome, cytobands } = await load();
  const p = (q) => parsePlace(q, genome, cytobands);
  const pick = (r) => r && { kind: r.kind, contig: r.contig, start: r.start, end: r.end };
  assert.deepEqual(pick(p("chr1:1,000-2,000")), { kind: "locus", contig: "chr1", start: 999, end: 2000 });
  assert.deepEqual(pick(p("1:1k-2k")), { kind: "locus", contig: "chr1", start: 999, end: 2000 });
  assert.deepEqual(pick(p("chr1 1000 2000")), { kind: "locus", contig: "chr1", start: 999, end: 2000 });
  assert.deepEqual(pick(p("chr1:2000-1000")), { kind: "locus", contig: "chr1", start: 999, end: 2000 });
  assert.deepEqual(pick(p("MT:100..200")), { kind: "locus", contig: "chrM", start: 99, end: 200 });
  assert.deepEqual(pick(p("chr1:0.1M-0.2m")), { kind: "locus", contig: "chr1", start: 99_999, end: 200_000 });
  assert.deepEqual(pick(p("chr1:350000-900000")), { kind: "locus", contig: "chr1", start: 349_999, end: 400_000 });
  assert.deepEqual(pick(p("chr1:500000-600000")), { kind: "locus", contig: "chr1", start: 399_999, end: 400_000 });
  assert.deepEqual(pick(p("chr1:5000")), { kind: "point", contig: "chr1", start: 4999, end: 5000 });
  assert.deepEqual(pick(p("chrM")), { kind: "contig", contig: "chrM", start: 0, end: 16_569 });
  assert.deepEqual(pick(p("chromosome 2")), { kind: "contig", contig: "chr2", start: 0, end: 120_000 });
  assert.deepEqual(pick(p("1p11")), { kind: "band", contig: "chr1", start: 90_000, end: 190_000 });
  assert.deepEqual(pick(p("chr1q21")), { kind: "band", contig: "chr1", start: 215_000, end: 400_000 });
  assert.equal(p("chr3:1-100"), null);
  assert.equal(p("PIXL1"), null);
  assert.equal(p("1p99"), null);
});

test("sequence helpers", () => {
  assert.equal(revcomp("ACGTNacgtn"), "nacgtNACGT");
  assert.equal(revcomp("RYKM"), "KMRY");
  assert.equal(translate("ATGTGGTAA"), "MW*");
  assert.equal(translate("ATGTGAATAAGA"), "M*IR");
  assert.equal(translate("ATGTGAATAAGA", true), "MWM*"); // vertebrate mitochondrial code
  assert.equal(translate("ATGNNNTA"), "MX");
  const codes = Uint8Array.from([1, 2, 4, 8, 15, 0x11, 0x14]);
  assert.equal(lettersOf(codes), "ACGTNAG");
  assert.equal(lettersOf(codes, true), "ACGTNag");
  assert.equal(gcFraction(codes), 3 / 6);
  assert.equal(wrapFasta("x", "ACGTACGT", 3), ">x\nACG\nTAC\nGT\n");
  // a two-exon transcript on the minus strand
  const tx = { strand: -1, exons: [[10, 13], [20, 23]], cds: [11, 22] };
  const region = Uint8Array.from({ length: 30 }, (_, i) => [1, 2, 4, 8][i % 4]);
  const { seq, map } = spliced(region, 0, tx);
  assert.equal(seq, revcomp(lettersOf(region.subarray(10, 13)) + lettersOf(region.subarray(20, 23))));
  assert.deepEqual(map, [22, 21, 20, 12, 11, 10]);
  assert.deepEqual(cdsRange(tx, map), [1, 5]);
});

test("palettes", () => {
  for (const name of Object.keys(PALETTES)) {
    const lut = buildLut(name);
    assert.equal(lut.rgb.length, 768);
    assert.equal(lut.packed.length, 256);
    // N stays the palette's N colour whatever the stretch
    const out = binColour(lut, 0.9, 1, 0, [0, 0, 0], 2.6);
    assert.deepEqual(out.map(Math.round), PALETTES[name].rgb.N, name);
    // soft-masked bases fade, in the palettes that fade them (GPxI greys doesn't)
    const a = [...lut.rgb.subarray(3, 6)];
    const masked = [...lut.rgb.subarray(3 * 0x11, 3 * 0x11 + 3)];
    if (PALETTES[name].fade) assert.notDeepEqual(a, masked, name);
    else assert.deepEqual(a, masked, name);
  }
  // names every object has are not palettes (a link with pal=toString must not break the site)
  assert.equal(hasPalette("toString"), false);
  assert.equal(hasPalette("__proto__"), false);
  assert.equal(hasPalette("neon"), true);
  assert.equal(buildLut("toString").palette, PALETTES.classic);
  assert.equal(inkOn([255, 255, 255]), "#000");
  assert.equal(inkOn([0, 0, 0]), "#fff");
});

test("utilities", async () => {
  assert.equal(fmtBp(950), "950 bp");
  assert.match(fmtBp(1_500), /^1\.5 kb$/);
  assert.match(fmtBp(83_300_000), /^83\.3 Mb$/);
  assert.equal(upperBound([1, 2, 2, 5], 2), 3);
  assert.equal(lowerBound([1, 2, 2, 5], 2), 1);
  const lru = new LRU(10, (v) => v.length);
  lru.set("a", "12345");
  lru.set("b", "12345");
  lru.get("a");
  lru.set("c", "123");
  assert.equal(lru.get("b"), undefined); // least recently used goes first
  assert.equal(lru.get("a"), "12345");
  let release;
  const p = track(new Promise((r) => (release = r)), "DNA");
  assert.deepEqual(busyWith(), ["DNA"]);
  release();
  await p;
  assert.deepEqual(busyWith(), []);
});

// ---------------------------------------------------------------- live Ensembl tracks, with fetch mocked

function mockFetch(handler) {
  const real = globalThis.fetch;
  const calls = [];
  globalThis.fetch = (url, opts = {}) => {
    const call = { url: String(url), opts };
    calls.push(call);
    return new Promise((resolve, reject) => {
      call.resolve = (rows) => resolve(new Response(JSON.stringify(rows), { status: 200 }));
      opts.signal?.addEventListener("abort", () => reject(Object.assign(new Error("aborted"), { name: "AbortError" })));
      handler?.(call);
    });
  };
  return { calls, restore: () => (globalThis.fetch = real) };
}

const region = (url) => {
  const m = /region\/human\/(\w+):(\d+)-(\d+)/.exec(url);
  return { chrom: m[1], start: +m[2] - 1, end: +m[3] };
};

test("live tracks: the request covers the whole window", async () => {
  const f = mockFetch();
  try {
    const live = new LiveTracks();
    live.enable("regulatory", true);
    for (const len of [1000, 2_900_000, LIMITS.regulatory]) {
      const d = { contig: "chr1", start: 10_000_000, end: 10_000_000 + len, length: 248_956_422 };
      live.fetchKind("regulatory", d);
      const r = region(f.calls.at(-1).url);
      assert.ok(r.start <= d.start && r.end >= d.end, `window of ${len}`);
      assert.ok(r.end - r.start <= LIMITS.regulatory);
    }
    // chrM is MT at Ensembl
    live.fetchKind("regulatory", { contig: "chrM", start: 0, end: 16_569, length: 16_569 });
    assert.equal(region(f.calls.at(-1).url).chrom, "MT");
  } finally {
    f.restore();
  }
});

test("live tracks: answers for an old window, or a track switched off, are dropped", async () => {
  const f = mockFetch();
  try {
    const live = new LiveTracks();
    live.enable("variation", true);
    const a = { contig: "chr7", start: 100_000, end: 110_000, length: 159_345_973 };
    const b = { contig: "chr7", start: 500_000, end: 510_000, length: 159_345_973 };
    const pa = live.fetchKind("variation", a);
    const pb = live.fetchKind("variation", b); // moved on: the request for a is aborted
    f.calls[1].resolve([{ id: "rsB", start: 505_000, end: 505_000, alleles: ["A", "G"], clinical_significance: ["pathogenic"] }]);
    f.calls[0].resolve([{ id: "rsA", start: 105_000, end: 105_000, alleles: ["C", "T"] }]);
    await Promise.all([pa, pb]);
    assert.deepEqual(live.items("variation", "chr7", 0, 1e9).map((v) => v.id), ["rsB"]);
    assert.equal(live.items("variation", "chr7", 0, 1e9)[0].severe, true);
    // switched off while a request is on its way: nothing is kept or drawn
    const c = { contig: "chr7", start: 900_000, end: 910_000, length: 159_345_973 };
    const pc = live.fetchKind("variation", c);
    live.enable("variation", false);
    f.calls[2].resolve([{ id: "rsC", start: 905_000, end: 905_000, alleles: ["G", "T"] }]);
    await pc;
    assert.equal(live.data.variation, null);
    assert.deepEqual(live.items("variation", "chr7", 0, 1e9), []);
    // regulatory elements are typed by their description
    live.enable("regulatory", true);
    const pr = live.fetchKind("regulatory", { contig: "chr7", start: 0, end: 50_000, length: 159_345_973 });
    f.calls.at(-1).resolve([{ id: "ENSR1", start: 101, end: 400, feature_type: "regulatory", description: "enhancer" }]);
    await pr;
    assert.deepEqual(live.items("regulatory", "chr7", 0, 1e9), [{ id: "ENSR1", start: 100, end: 400, type: "enhancer" }]);
  } finally {
    f.restore();
  }
});

// ---------------------------------------------------------------- pxil files and zips, checked with Python (Pillow, zipfile)

test("pxil files are PNGs whose pixels read back as the sequence", async () => {
  const out = path.join(ROOT, "tests", ".fixture", "out");
  fs.mkdirSync(out, { recursive: true });
  // a gene with soft-masked bases, ambiguity codes and Ns, and a protein with selenocysteine and a stop
  const dna = "ACGTacgtNNRYacgtTTGA".repeat(13) + "AC";
  const protein = "MEEPQSDPSVUEPPLSQ*";
  const files = [
    { name: "x.gene.pxil", seq: dna, alphabet: DNA, width: 30, kind: "gene" },
    { name: "x.protein.pxil", seq: protein, alphabet: PROTEIN, width: 7, kind: "protein" },
  ];
  for (const f of files) {
    const colours = [...f.alphabet].map((_, i) => [i * 8, 255 - i * 8, 100]);
    f.bytes = await pxilFile({ indices: indicesOf(f.seq, f.alphabet), alphabet: f.alphabet, width: f.width, colours, meta: { kind: f.kind, id: "ENSX1", name: "X" }, sequence: f.seq });
    fs.writeFileSync(path.join(out, f.name), f.bytes);
  }
  fs.writeFileSync(path.join(out, "x.zip"), zip(files.map((f) => ({ name: f.name, bytes: f.bytes }))));
  const py = `
import json, sys, zipfile, hashlib
from PIL import Image
out = {}
for name in sys.argv[2:]:
    im = Image.open(sys.argv[1] + "/" + name)
    meta = json.loads(im.info["pxil"])
    idx = im.tobytes()[: meta["length"]]
    seq = "".join(meta["alphabet"][i] for i in idx)
    out[name] = {"mode": im.mode, "size": list(im.size), "meta": meta, "seq": seq, "pad": set(im.tobytes()[meta["length"]:]) <= {0},
                 "sha_ok": meta["sha256"] == hashlib.sha256(seq.encode()).hexdigest(), "bits": im.info.get("bits", 8)}
z = zipfile.ZipFile(sys.argv[1] + "/x.zip")
out["zip"] = {"names": z.namelist(), "bad": z.testzip(), "same": all(z.read(n) == open(sys.argv[1] + "/" + n, "rb").read() for n in z.namelist())}
print(json.dumps(out))
`;
  const r = JSON.parse(execFileSync(process.env.PYTHON || "python3", ["-c", py, out, ...files.map((f) => f.name)], { encoding: "utf8" }));
  for (const f of files) {
    const g = r[f.name];
    assert.equal(g.mode, "P");
    assert.deepEqual(g.size, [f.width, Math.ceil(f.seq.length / f.width)]);
    assert.equal(g.seq, f.seq);
    assert.equal(g.pad, true);
    assert.equal(g.sha_ok, true);
    assert.equal(g.meta.pxil, "0.1");
    assert.equal(g.meta.kind, f.kind);
    assert.equal(g.meta.length, f.seq.length);
  }
  assert.deepEqual(r.zip, { names: ["x.gene.pxil", "x.protein.pxil"], bad: null, same: true });
  assert.equal(crc32(new TextEncoder().encode("123456789")), 0xcbf43926); // the standard check value
  assert.equal(safeName("TP53-201 / α"), "TP53-201_");
  assert.throws(() => indicesOf("ACGX", DNA));
});
