// Genes: a static index built from Ensembl's gene set (every transcript, HGNC
// names and aliases), plus live lookups in the Ensembl REST API for transcript
// IDs and symbols the index doesn't know, and for variant IDs.

import { Emitter, fetchJson, fetchJsonGz, fetchRange, gunzip, upperBound } from "./util.js";

export const REST = "https://rest.ensembl.org";
const REST_HEADERS = { headers: { Accept: "application/json" } };

// transcript flags in the transcripts files (pipeline/genes.py)
const CANONICAL = 1;
const MANE_SELECT = 2;
const MANE_CLINICAL = 4;
const BASIC = 8;

function normaliseGene(row, contig) {
  const [id, name, start, end, strand, biotype, description, tx, xref, where] = row;
  const g = { id, name, contig, start, end, strand, biotype, description: description || "", xref: null, tx: null, txIndex: where || null, txCount: where ? where[2] : tx ? 1 : 0 };
  if (xref) g.xref = { hgnc: xref[0], entrez: xref[1], uniprot: xref[2], omim: xref[3] };
  if (tx) {
    const [tid, tname, tbio, mane, ts, te, flat, cds, incomplete, sec] = tx;
    const exons = [];
    for (let i = 0; i < flat.length; i += 2) exons.push([ts + flat[i], ts + flat[i + 1]]);
    g.tx = {
      id: tid,
      name: tname,
      biotype: tbio,
      mane: !!mane,
      canonical: true,
      start: ts,
      end: te,
      strand,
      exons,
      cds: cds ? [cds[0], cds[1]] : null,
      frame: cds ? cds[2] : 0,
      incomplete: incomplete || "",
      sec: sec || [], // selenocysteine codons (genomic position of their first base)
    };
  }
  return g;
}

/** One transcript from a transcripts file (coordinates relative to the gene's start). */
function normaliseTranscript(row, gene) {
  const [id, name, biotype, flags, s, e, flat, cds, incomplete, sec] = row;
  const g0 = gene.start;
  const exons = [];
  for (let i = 0; i < flat.length; i += 2) exons.push([g0 + flat[i], g0 + flat[i + 1]]);
  return {
    id,
    name,
    biotype,
    canonical: !!(flags & CANONICAL),
    mane: !!(flags & MANE_SELECT),
    maneClinical: !!(flags & MANE_CLINICAL),
    basic: !!(flags & BASIC),
    start: g0 + s,
    end: g0 + e,
    strand: gene.strand,
    exons,
    cds: cds ? [g0 + cds[0], g0 + cds[1]] : null,
    frame: cds ? cds[2] : 0,
    incomplete: incomplete || "",
    sec: (sec || []).map((p) => g0 + p),
  };
}

export class GeneDB extends Emitter {
  constructor(manifest, base) {
    super();
    this.base = base;
    this.manifest = manifest;
    this.contigFiles = new Map(manifest.sequence.contigs.filter((c) => c.genes).map((c) => [c.name, c.genes]));
    this.txFiles = new Map(manifest.sequence.contigs.filter((c) => c.transcripts).map((c) => [c.name, c.transcripts]));
    this.byContig = new Map(); // contig -> {genes, starts, maxEnd}
    this.loading = new Map();
    this.indexPromise = null;
    this.aliasPromise = null;
    this.txCache = new Map();
  }

  // ---------------------------------------------------------------- per contig
  loadContig(contig) {
    if (this.byContig.has(contig)) return Promise.resolve(this.byContig.get(contig));
    if (!this.contigFiles.has(contig)) {
      const empty = { genes: [], starts: [], maxEnd: [] };
      this.byContig.set(contig, empty);
      return Promise.resolve(empty);
    }
    if (!this.loading.has(contig)) {
      const p = fetchJsonGz(new URL(this.contigFiles.get(contig), this.base).href, "genes").then((d) => {
        const genes = d.genes.map((r) => normaliseGene(r, contig));
        genes.sort((a, b) => a.start - b.start);
        const starts = genes.map((g) => g.start);
        const maxEnd = [];
        let m = -1;
        for (const g of genes) maxEnd.push((m = Math.max(m, g.end)));
        const entry = { genes, starts, maxEnd, byId: new Map(genes.map((g) => [g.id, g])) };
        this.byContig.set(contig, entry);
        this.loading.delete(contig);
        this.emit("loaded", contig);
        return entry;
      });
      p.catch(() => this.loading.delete(contig));
      this.loading.set(contig, p);
    }
    return this.loading.get(contig);
  }

  /** Genes overlapping [start, end) if the contig is loaded (and starts loading it if not). */
  overlapping(contig, start, end) {
    const e = this.byContig.get(contig);
    if (!e) {
      this.loadContig(contig).catch(() => {});
      return [];
    }
    const out = [];
    let i = upperBound(e.starts, end - 1) - 1;
    for (; i >= 0 && e.maxEnd[i] > start; i--) if (e.genes[i].end > start) out.push(e.genes[i]);
    return out.reverse();
  }

  /** Genes at one position, shortest first. */
  at(contig, pos) {
    return this.overlapping(contig, pos, pos + 1).sort((a, b) => a.end - a.start - (b.end - b.start));
  }

  async gene(id, contig) {
    if (contig) {
      const e = await this.loadContig(contig);
      if (e.byId && e.byId.has(id)) return e.byId.get(id);
    }
    const idx = await this.index();
    const hit = idx.byId.get(id);
    if (!hit) return null;
    const e = await this.loadContig(hit[2]);
    return e.byId.get(id) || null;
  }

  // ---------------------------------------------------------------- search index
  index() {
    if (!this.indexPromise) {
      this.indexPromise = fetchJsonGz(new URL(this.manifest.genes.index, this.base).href, "gene index").then((rows) => {
        const byName = new Map();
        const byId = new Map();
        for (const r of rows) {
          const k = r[0].toLowerCase();
          if (!byName.has(k)) byName.set(k, []);
          byName.get(k).push(r);
          byId.set(r[1], r);
        }
        const names = [...byName.keys()].sort();
        return { rows, byName, byId, names };
      });
      this.indexPromise.catch(() => (this.indexPromise = null));
    }
    return this.indexPromise;
  }

  aliases() {
    if (!this.aliasPromise) {
      this.aliasPromise = fetchJsonGz(new URL(this.manifest.genes.aliases, this.base).href, "aliases").then((rows) => {
        const m = new Map();
        for (const [alias, syms] of rows) m.set(alias.toLowerCase(), { alias, syms });
        return m;
      });
      this.aliasPromise.catch(() => (this.aliasPromise = null));
    }
    return this.aliasPromise;
  }

  /** Suggestions for a typed query: exact symbols, then prefixes, then aliases. */
  async suggest(q, limit = 12) {
    q = q.trim();
    if (!q) return [];
    const idx = await this.index();
    const k = q.toLowerCase();
    const out = [];
    const seen = new Set();
    const add = (r, via) => {
      if (seen.has(r[1]) || out.length >= limit) return;
      seen.add(r[1]);
      out.push({ name: r[0], id: r[1], contig: r[2], start: r[3], end: r[4], strand: r[5], biotype: r[6], via });
    };
    if (/^ensg\d+/i.test(q)) {
      const r = idx.byId.get(q.toUpperCase().split(".")[0]);
      if (r) add(r);
    }
    for (const r of idx.byName.get(k) || []) add(r);
    const al = this.aliasPromise ? await this.aliases() : null;
    if (al && al.has(k)) for (const sym of al.get(k).syms) for (const r of idx.byName.get(sym.toLowerCase()) || []) add(r, al.get(k).alias);
    // prefix matches, protein-coding first
    let i = upperBound(idx.names, k) - 1;
    if (i < 0 || idx.names[i] !== k) i++;
    const pref = [];
    for (; i < idx.names.length && idx.names[i].startsWith(k) && pref.length < 400; i++) for (const r of idx.byName.get(idx.names[i])) pref.push(r);
    pref.sort((a, b) => (a[6] === "protein_coding" ? 0 : 1) - (b[6] === "protein_coding" ? 0 : 1) || a[0].length - b[0].length || a[0].localeCompare(b[0]));
    for (const r of pref) add(r);
    return out;
  }

  /** Resolve a free-text gene query to one gene (index, then aliases, then Ensembl REST). */
  async resolve(q) {
    const idx = await this.index();
    const k = q.trim().toLowerCase();
    if (/^ensg\d+/i.test(q)) {
      const r = idx.byId.get(q.trim().toUpperCase().split(".")[0]);
      if (r) return { row: r };
    }
    const exact = idx.byName.get(k);
    if (exact) return { row: pickBest(exact) };
    const al = await this.aliases();
    if (al.has(k)) {
      const rows = al.get(k).syms.flatMap((s) => idx.byName.get(s.toLowerCase()) || []);
      if (rows.length) return { row: pickBest(rows), via: al.get(k).alias };
    }
    if (/^enst\d+/i.test(q)) {
      const t = await fetchJson(`${REST}/lookup/id/${encodeURIComponent(q.trim().split(".")[0])}?content-type=application/json`, "Ensembl", REST_HEADERS);
      const r = idx.byId.get(t.Parent);
      if (r) return { row: r, transcript: t.id };
    }
    // last resort: Ensembl's cross-references (other synonyms, RefSeq IDs ...)
    try {
      const xs = await fetchJson(`${REST}/xrefs/symbol/homo_sapiens/${encodeURIComponent(q.trim())}?content-type=application/json`, "Ensembl", REST_HEADERS);
      for (const x of xs) if (x.type === "gene" && idx.byId.has(x.id)) return { row: idx.byId.get(x.id), via: q.trim() };
    } catch {
      /* not found */
    }
    return null;
  }

  // ---------------------------------------------------------------- transcripts
  /** All transcripts of a gene, canonical first: from the site's own data, or failing that from Ensembl. */
  transcripts(gene) {
    if (this.txCache.has(gene.id)) return this.txCache.get(gene.id);
    const p = gene.txIndex && this.txFiles.has(gene.contig) ? this.storedTranscripts(gene) : this.ensemblTranscripts(gene);
    this.txCache.set(gene.id, p);
    p.catch(() => this.txCache.delete(gene.id));
    return p;
  }

  /** One range request into the contig's transcripts file: the gene's own gzip member. */
  async storedTranscripts(gene) {
    const [offset, length] = gene.txIndex;
    const url = new URL(this.txFiles.get(gene.contig), this.base).href;
    const bytes = await fetchRange(url, offset, offset + length, "transcripts");
    const rows = JSON.parse(new TextDecoder().decode(await gunzip(bytes)));
    return { description: "", txs: rows.map((r) => normaliseTranscript(r, gene)) };
  }

  /** All transcripts of a gene from the Ensembl REST API (for data built without transcripts files). */
  ensemblTranscripts(gene) {
    return fetchJson(`${REST}/lookup/id/${gene.id}?expand=1;utr=1;content-type=application/json`, "Ensembl transcripts", REST_HEADERS).then((d) => {
      const txs = (d.Transcript || []).map((t) => {
        const exons = (t.Exon || []).map((e) => [e.start - 1, e.end]).sort((a, b) => a[0] - b[0]);
        const tr = t.Translation;
        return {
          id: t.id,
          name: t.display_name || t.id,
          biotype: t.biotype,
          canonical: !!t.is_canonical,
          mane: gene.tx && gene.tx.id === t.id ? gene.tx.mane : false,
          start: t.start - 1,
          end: t.end,
          strand: t.strand,
          exons,
          cds: tr ? [tr.start - 1, tr.end] : null, // Ensembl translations include the stop codon
          frame: 0,
          length: t.length,
          protein: tr ? tr.length : 0,
        };
      });
      txs.sort((a, b) => b.canonical - a.canonical || (b.cds ? 1 : 0) - (a.cds ? 1 : 0) || a.name.localeCompare(b.name));
      return { description: d.description || "", txs };
    });
  }

  async variant(rsid) {
    const d = await fetchJson(`${REST}/variation/human/${encodeURIComponent(rsid)}?content-type=application/json`, "Ensembl variation", REST_HEADERS);
    const m = (d.mappings || []).find((x) => x.assembly_name === "GRCh38" && /^([0-9]+|X|Y|MT)$/.test(x.seq_region_name)) || (d.mappings || [])[0];
    if (!m) throw new Error(`${rsid} has no GRCh38 position`);
    return { id: d.name, contig: m.seq_region_name, start: m.start - 1, end: m.end, alleles: m.allele_string, consequence: d.most_severe_consequence, maf: d.MAF, minor: d.minor_allele };
  }

  /** NCBI's summary of a gene (RefSeq), via MyGene.info. */
  async summary(gene) {
    if (!gene.xref || !gene.xref.entrez) return null;
    const d = await fetchJson(`https://mygene.info/v3/gene/${gene.xref.entrez}?fields=summary,name`, "NCBI summary");
    return d && d.summary ? d.summary : null;
  }
}

function pickBest(rows) {
  const order = (r) => (/^(1?[0-9]|2[0-2]|X|Y|M)$/.test(r[2].replace(/^chr/, "")) ? 0 : 1) * 2 + (r[6] === "protein_coding" ? 0 : 1);
  return [...rows].sort((a, b) => order(a) - order(b))[0];
}
