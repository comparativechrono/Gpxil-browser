// Gpxil browser: wires the pieces together, keeps the URL in step with the
// view, and handles the panel, search, keyboard and the landing screen.

import { Genome } from "./gpxil.js";
import { Overview } from "./overview.js";
import { GeneDB } from "./genes.js";
import { PixelView, SCALES } from "./view.js";
import { Karyotype, bandAt } from "./karyotype.js";
import { Ideogram } from "./ideogram.js";
import { Search, parsePlace } from "./search.js";
import { GenePanel } from "./genepanel.js";
import { Landing } from "./landing.js";
import { DnaLoader } from "./loader.js";
import { LIMITS, LiveTracks } from "./live.js";
import { TranscriptLanes } from "./lanes.js";
import { safeName } from "./pxil.js";
import { IUPAC, PALETTES, buildLut, hasPalette, inkOn } from "./palettes.js";
import { gcFraction, lettersOf, wrapFasta } from "./seq.js";
import { clamp, debounce, escapeHtml, fetchJson, fetchJsonGz, fmt, fmtBp, storageGet, storageSet } from "./util.js";

const DATA = new URL("data/grch38/", location.href);
const $ = (s) => document.querySelector(s);
const NAV_KEYS = new Set(["ArrowDown", "ArrowUp", "ArrowLeft", "ArrowRight", "PageDown", "PageUp", " ", "Home", "End", "+", "=", "-", "_"]);

// display options are remembered in this browser; the palette is kept on its own
function savedOptions() {
  const o = storageGet("gpxil.opts", {});
  if (!o || typeof o !== "object") return {};
  delete o.palette;
  return o;
}

function saveOptions(opts) {
  const { palette, ...rest } = opts;
  void palette;
  storageSet("gpxil.opts", rest);
}

class App {
  constructor({ manifest, karyotype, cytobands }) {
    this.manifest = manifest;
    this.cytobands = cytobands;
    this.genome = new Genome(manifest, DATA);
    this.overview = new Overview(manifest, DATA);
    this.genes = new GeneDB(manifest, DATA);
    this.palette = storageGet("gpxil.palette", "classic");
    if (!hasPalette(this.palette)) this.palette = "classic";
    this.view = new PixelView({ canvas: $("#view"), gutter: $("#gutter"), genome: this.genome, overview: this.overview, genes: this.genes });
    this.karyo = new Karyotype({ canvas: $("#karyo"), genome: this.genome, karyotype, cytobands, genes: this.genes, overview: this.overview });
    this.ideo = new Ideogram({ canvas: $("#ideogram"), cytobands, karyotype });
    this.search = new Search({ input: $("#search"), list: $("#suggest"), genome: this.genome, genes: this.genes, cytobands });
    this.panel = new GenePanel({ root: $("#gene-panel"), genome: this.genome, genes: this.genes, cytobands });
    const release = /release (\d+)/.exec(manifest.sources.sequence || "");
    this.panel.source = release ? `GRCh38, Ensembl release ${release[1]}` : "GRCh38";
    this.live = new LiveTracks();
    this.view.live = this.live;
    this.lanes = new TranscriptLanes({ canvas: $("#lanes"), genome: this.genome, lut: () => this.view.lut, gutterWidth: () => $("#gutter").getBoundingClientRect().width || 118 });
    this.highlight = null;
    this.screen = "genome";
    this.applying = false;
    this.routeToken = 0; // tokens: only the latest of each kind of request is applied
    this.geneToken = 0;
    this.selToken = 0;
    this.pointer = null; // where the mouse last was over the view, for refreshing the tooltip
    this.bindView();
    this.bindPanels();
    this.bindKeys();
    this.fillContigSelect();
    this.setPalette(this.palette, false);
    this.view.setOptions(savedOptions());
    this.syncOptionInputs();
    this.panel.clear();
    this.renderBookmarks();
    $("#gene-count").textContent = `${fmt(manifest.genes.count)} genes`;
    $("#data-sources").innerHTML = Object.entries(manifest.sources)
      .map(([k, v]) => `<b>${escapeHtml(k)}</b>: ${escapeHtml(v)}`)
      .join("<br>");
    this.genes.on("loaded", () => this.view.schedule());
  }

  // ------------------------------------------------------------------ routing
  start() {
    addEventListener("popstate", () => this.route(location.hash));
    this.route(location.hash, true);
  }

  async route(hash, first = false) {
    let h = (hash || "").replace(/^#/, "");
    try {
      h = decodeURIComponent(h);
    } catch {
      /* a stray % in a hand-made link: use it as it is */
    }
    const [loc, ...rest] = h.split("&");
    const p = new URLSearchParams(rest.join("&"));
    const token = ++this.routeToken;
    this.applying = true;
    try {
      if (hasPalette(p.get("pal"))) this.setPalette(p.get("pal"), false);
      const wrap = p.get("wrap");
      if (wrap && /^\d+$/.test(wrap)) this.view.setWrap("fixed", +wrap);
      else if (wrap === "auto" || (!wrap && this.view.state.wrap !== "auto")) this.view.setWrap("auto");
      this.syncOptionInputs();
      if (!loc || loc === "genome") return this.showGenome(false);
      const place = parsePlace(loc, this.genome, this.cytobands);
      if (!place) {
        if (first) this.toast(`Couldn't read the place “${loc}”`);
        return this.showGenome(false);
      }
      const c = this.genome.contig(place.contig);
      if (place.kind === "point") this.showContig(place.contig, place.start - 60, place.end + 60, false);
      else {
        // a window this browser wrote (or any typed interval) comes back at the same zoom;
        // a whole chromosome or a band gets a margin
        const whole = place.start === 0 && place.end === c.length;
        this.showContig(place.contig, place.start, place.end, false, { exact: place.kind === "locus" && !whole });
      }
      if (place.kind === "point") {
        this.view.setSelection([place.start, place.end]);
        this.selectionInfo([place.start, place.end]);
      }
      const gene = p.get("gene");
      if (gene) await this.selectGene(gene, place.contig, { mode: p.get("mode") || "gene", transcript: p.get("tx"), zoom: false, push: false });
      else if (this.highlight) this.clearGene(false);
    } finally {
      if (token === this.routeToken) this.applying = false;
    }
  }

  hash() {
    if (this.screen === "genome") return "#genome";
    const d = this.view.describe();
    if (!d) return "#genome";
    let h = `#${d.contig}:${d.start + 1}-${d.end}`;
    if (this.highlight && this.highlight.gene.contig === d.contig) {
      h += `&gene=${this.highlight.gene.id}`;
      if (this.highlight.mode !== "gene") h += `&mode=${this.highlight.mode}`;
      if (this.highlight.tx && !this.highlight.tx.canonical) h += `&tx=${this.highlight.tx.id}`;
    }
    if (this.view.state.wrap === "fixed") h += `&wrap=${this.view.state.width}`;
    if (this.palette !== "classic") h += `&pal=${this.palette}`;
    return h;
  }

  pushState() {
    if (this.applying) return;
    const h = this.hash();
    if (h !== location.hash) history.pushState(null, "", h);
  }

  replaceState = debounce(() => {
    if (this.applying) return;
    const h = this.hash();
    if (h !== location.hash) history.replaceState(null, "", h);
  }, 300);

  // ------------------------------------------------------------------ screens
  showGenome(push = true) {
    this.screen = "genome";
    this.syncLanes();
    $("#genome-view").hidden = false;
    $("#chrom-view").hidden = true;
    $("#contig-select").value = "";
    this.karyo.setPalette(this.palette);
    this.karyo.show();
    document.title = "Gpxil Browser · GRCh38";
    if (push) this.pushState();
  }

  showContig(name, start = null, end = null, push = true, opts = {}) {
    const c = this.genome.contig(name);
    if (!c) return this.toast(`No sequence called ${name}`);
    this.screen = "contig";
    $("#genome-view").hidden = true;
    $("#chrom-view").hidden = false;
    this.view.resize();
    this.view.show(c.name, start, end, opts);
    this.syncLanes();
    $("#contig-select").value = c.name;
    document.title = `${c.name} · Gpxil Browser`;
    if (push) this.pushState();
  }

  fillContigSelect() {
    const sel = $("#contig-select");
    const chroms = this.manifest.sequence.contigs.filter((c) => c.group === "chromosome");
    const other = this.manifest.sequence.contigs.filter((c) => c.group !== "chromosome");
    sel.innerHTML =
      `<option value="">chromosome…</option>` +
      chroms.map((c) => `<option value="${c.name}">${c.name} · ${fmtBp(c.length)}</option>`).join("") +
      (other.length ? `<optgroup label="Unplaced and unlocalised">${other.map((c) => `<option value="${c.name}">${c.name} · ${fmtBp(c.length)}</option>`).join("")}</optgroup>` : "");
    sel.addEventListener("change", () => sel.value && this.showContig(sel.value));
  }

  // ------------------------------------------------------------------ view events
  bindView() {
    const v = this.view;
    v.on("navigate", (d) => {
      if (!d) return;
      this.ideo.set(this.genome.contig(d.contig), [d.start, d.end], this.markers());
      $("#st-win").textContent = `${d.contig}:${fmt(d.start + 1)}-${fmt(d.end)} · ${fmtBp(d.end - d.start)}`;
      $("#st-zoom").textContent = v.statusText();
      this.live.update(d);
      this.replaceState();
    });
    this.live.on("update", () => v.schedule());
    this.live.on("error", (m) => this.toast(m));
    v.on("hover", (pos, x, y) => {
      this.pointer = pos === null ? null : [x, y];
      this.hoverInfo(pos, x, y);
    });
    v.on("data", () => {
      if (v.hover !== null) this.hoverInfo(v.hover, ...(this.pointer || []));
    });
    v.on("click", (pos) => {
      if (pos === null) return;
      const variant = this.live.items("variation", v.state.contig, pos, pos + 1)[0];
      if (variant) {
        v.setSelection([variant.start, variant.end]);
        return this.selectionInfo([variant.start, variant.end], { id: variant.id, alleles: variant.alleles, consequence: variant.consequence, clinical: variant.clinical });
      }
      const genes = this.genes.at(v.state.contig, pos);
      if (!genes.length) return;
      const pick = genes.find((g) => g.biotype === "protein_coding") || genes[0];
      if (!this.highlight || this.highlight.gene.id !== pick.id) this.selectGene(pick.id, pick.contig, { zoom: false });
    });
    v.on("select", (sel) => this.selectionInfo(sel));
    // Expanded mode: the lanes take the viewer's place
    const lanes = this.lanes;
    lanes.on("navigate", (d) => {
      this.ideo.set(this.genome.contig(d.contig), [d.start, d.end], this.markers());
      $("#st-win").textContent = `${d.contig}:${fmt(d.start + 1)}-${fmt(d.end)} · ${fmtBp(d.end - d.start)}`;
      $("#st-zoom").textContent = `${fmt(d.transcripts)} transcripts · ${d.bpp >= 1 ? `${fmt(d.bpp)} bases per pixel` : `${(1 / d.bpp).toFixed(1)} px per base`}`;
    });
    lanes.on("hover", (text, x, y, pos) => {
      this.tooltip(text, x, y);
      $("#st-pos").textContent = text && pos !== null && pos !== undefined ? `${lanes.gene.contig}:${fmt(pos + 1)}` : "";
    });
    lanes.on("pick", (id) => {
      this.panel.pickTx(id);
      this.panel.setMode("transcript");
      this.tooltip(null);
      this.toast(`${this.panel.tx.name} · press 3 to see every transcript again`, 3500);
    });
    lanes.on("error", (e) => this.toast(`Couldn't load: ${e.message}`));
    v.on("error", (e) => this.toast(`Couldn't load: ${e.message}`));
    this.karyo.on("open", (contig, pos) => {
      if (pos === null) this.showContig(contig);
      else this.showContig(contig, Math.max(0, pos - 4e6), pos + 4e6);
    });
    this.karyo.on("hover", (text, x, y) => this.tooltip(text, x, y));
    this.ideo.on("seek", (pos) => {
      if (this.lanes.active) return this.lanes.centreOn(pos);
      const d = v.describe();
      const half = Math.floor((d.end - d.start) / 2);
      v.show(d.contig, pos - half, pos + half, { zoom: v.state.zi });
    });
    this.ideo.on("hover", (text, x, y) => this.tooltip(text, x, y));
    this.search.on("place", (p) => {
      if (p.kind === "point") {
        this.showContig(p.contig, p.start - 60, p.end + 60);
        v.setSelection([p.start, p.end]);
        this.selectionInfo([p.start, p.end]);
      } else this.showContig(p.contig, p.start, p.end);
    });
    this.search.on("gene", (id, contig, tx) => this.selectGene(id, contig, { mode: tx ? "transcript" : "gene", transcript: tx }));
    this.search.on("variant", (id) => this.gotoVariant(id));
    this.search.on("notfound", (q) => this.toast(`Nothing called “${q}” in GRCh38 genes`));
    window.addEventListener("resize", () => {
      if (this.screen === "contig") v.resize();
    });
  }

  markers() {
    const out = [];
    if (this.highlight && this.highlight.gene.contig === this.view.state.contig) out.push({ start: this.highlight.gene.start, end: this.highlight.gene.end, colour: "rgba(255,62,165,0.9)" });
    if (this.view.selection) out.push({ start: this.view.selection[0], end: this.view.selection[1], colour: "rgba(255,225,77,0.9)" });
    return out;
  }

  // ------------------------------------------------------------------ genes
  async selectGene(id, contig, { mode = "gene", transcript = null, zoom = true, push = true } = {}) {
    const token = ++this.geneToken;
    let gene;
    try {
      gene = await this.genes.gene(id, contig);
    } catch (e) {
      return this.toast(`Couldn't load genes: ${e.message}`);
    }
    if (token !== this.geneToken) return; // another gene was picked meanwhile
    if (!gene) return this.toast(`Gene ${id} not found`);
    if (zoom || this.view.state.contig !== gene.contig || this.screen !== "contig") {
      const pad = Math.max(200, Math.round((gene.end - gene.start) * 0.08));
      this.showContig(gene.contig, gene.start - pad, gene.end + pad, false);
    }
    this.showTab("gene");
    this.openPanel(true);
    await this.panel.open(gene, { mode, transcript });
    if (push && token === this.geneToken) this.pushState();
  }

  clearGene(push = true) {
    this.highlight = null;
    this.view.setHighlight(null);
    this.panel.clear();
    this.syncLanes();
    this.refreshIdeo();
    if (push) this.pushState();
  }

  refreshIdeo() {
    const d = this.lanes.active ? this.lanes.describe() : this.view.describe();
    if (d) this.ideo.set(this.genome.contig(d.contig), [d.start, d.end], this.markers());
  }

  /** Expanded mode shows every transcript as lanes over the viewer; other modes, the genome. */
  syncLanes() {
    const h = this.highlight;
    const on = !!h && h.mode === "expanded" && this.screen === "contig" && h.gene.contig === this.view.state.contig;
    if (on) {
      if (!this.lanes.active) this.tooltip(null);
      this.view.paused = true;
      this.lanes.show(h.gene, h.txs || (h.gene.tx ? [h.gene.tx] : []), h.tx ? h.tx.id : null);
    } else if (this.lanes.active) {
      this.tooltip(null);
      this.lanes.hide();
      this.view.paused = false;
      this.view.dirtyFrame = true;
      this.view.schedule();
      const d = this.view.describe();
      if (d) this.view.emit("navigate", d);
    }
  }

  async gotoVariant(id) {
    try {
      const v = await this.genes.variant(id);
      const contig = this.genome.resolve(v.contig);
      if (!contig) throw new Error(`${id} maps to ${v.contig}, which is not in this assembly`);
      this.showContig(contig, v.start - 40, v.end + 40);
      this.view.setSelection([v.start, Math.max(v.end, v.start + 1)]);
      this.selectionInfo([v.start, Math.max(v.end, v.start + 1)], v);
      this.toast(`${v.id} · ${v.alleles}${v.consequence ? ` · ${v.consequence.replace(/_/g, " ")}` : ""}${v.maf ? ` · MAF ${v.maf}` : ""}`, 5000);
    } catch (e) {
      this.toast(`Couldn't find ${id}: ${e.message}`);
    }
  }

  async randomGene() {
    const idx = await this.genes.index();
    const pool = idx.rows.filter((r) => r[6] === "protein_coding" && !r[0].startsWith("ENSG") && /^chr/.test(r[2]));
    const r = pool[Math.floor(Math.random() * pool.length)];
    this.toast(`🎲 ${r[0]}`, 1600);
    this.selectGene(r[1], r[2]);
  }

  stepGene(dir) {
    const d = this.view.describe();
    if (!d) return;
    const e = this.genes.byContig.get(d.contig);
    if (!e || !e.genes.length) return;
    const ref = this.highlight && this.highlight.gene.contig === d.contig ? this.highlight.gene.start : Math.floor((d.start + d.end) / 2);
    const coding = e.genes.filter((g) => g.biotype === "protein_coding");
    const list = coding.length ? coding : e.genes;
    const g = dir > 0 ? list.find((x) => x.start > ref) : [...list].reverse().find((x) => x.start < ref);
    if (g) this.selectGene(g.id, g.contig);
    else this.toast(dir > 0 ? "No more genes this way" : "That was the first gene");
  }

  // ------------------------------------------------------------------ info panels
  hoverInfo(pos, x, y) {
    const v = this.view;
    if (pos === null || pos === undefined) {
      $("#st-pos").textContent = "";
      this.tooltip(null);
      return;
    }
    const c = v.state.contig;
    const codes = v.codesFor(Math.max(0, pos - 12), Math.min(v.contig.length, pos + 13));
    const at = Math.min(12, pos);
    const code = codes ? codes[at] : null;
    const base = code === null ? null : IUPAC[code & 15];
    const band = bandAt(this.cytobands[c], pos);
    const genes = this.genes.at(c, pos);
    const L = v.layout();
    const span = L.bpp > 1 ? ` (+${fmt(L.bpp - 1)} bases)` : "";
    $("#st-pos").textContent = `${c}:${fmt(pos + 1)}${span}${base ? ` ${base}` : ""}`;
    const geneText = genes
      .slice(0, 3)
      .map((g) => `${g.name}${this.exonText(g, pos)}`)
      .join(", ");
    const vars = this.live.items("variation", c, pos, pos + 1);
    const regs = this.live.items("regulatory", c, pos, pos + 1);
    const liveText = [...vars.map((x) => `${x.id} ${x.alleles}${x.clinical ? ` (${x.clinical})` : ""}`), ...regs.map((r) => r.type.replace(/_/g, " "))].join(" · ");
    if (x === undefined) this.tooltip(null);
    else this.tooltip(`${c}:${fmt(pos + 1)}${base ? ` · ${base}${code & 0x10 ? " (repeat)" : ""}` : ""}${band ? ` · ${c.replace("chr", "")}${band[2]}` : ""}${geneText ? ` · ${geneText}` : ""}${liveText ? ` · ${liveText}` : ""}`, x, y, true);
    const lut = buildLut(this.palette);
    const chip = code === null ? "" : `<span class="base-chip" style="background:rgb(${lut.rgb[3 * code]},${lut.rgb[3 * code + 1]},${lut.rgb[3 * code + 2]});color:${inkOn([lut.rgb[3 * code], lut.rgb[3 * code + 1], lut.rgb[3 * code + 2]])}">${code & 0x10 ? base.toLowerCase() : base}</span>`;
    const peek = codes ? `${lettersOf(codes.subarray(0, at), true)}<b style="color:var(--yellow)">${lettersOf(codes.subarray(at, at + 1), true)}</b>${lettersOf(codes.subarray(at + 1), true)}` : "";
    $("#here").innerHTML = `
      <p>${chip}<b>${c}:${fmt(pos + 1)}</b>${span}</p>
      <dl class="kv">
        ${band ? `<dt>band</dt><dd>${c.replace("chr", "")}${band[2]} <small>(${band[3]})</small></dd>` : ""}
        ${code !== null ? `<dt>base</dt><dd>${base}${code & 0x10 ? " · in a repeat (soft-masked)" : ""}${(code & 15) === 15 ? " · unknown (N)" : ""}</dd>` : ""}
        <dt>genes</dt><dd class="here-genes">${genes.length ? genes.map((g) => `<a href="#" data-gene="${g.id}">${escapeHtml(g.name)}</a><small>${escapeHtml(g.biotype.replace(/_/g, " "))}${this.exonText(g, pos)}</small>`).join("<br>") : "<span class='empty'>none (intergenic)</span>"}</dd>
        ${vars.length ? `<dt>variants</dt><dd>${vars.map((x) => `<a target="_blank" rel="noopener" href="https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=${x.id}">${x.id}</a> ${escapeHtml(x.alleles)} <small>${escapeHtml(x.consequence.replace(/_/g, " "))}${x.clinical ? ` · ${escapeHtml(x.clinical)}` : ""}</small>`).join("<br>")}</dd>` : ""}
        ${regs.length ? `<dt>regulation</dt><dd>${regs.map((r) => `${escapeHtml(r.type.replace(/_/g, " "))} <small>${escapeHtml(r.id)}</small>`).join("<br>")}</dd>` : ""}
        <dt>row</dt><dd>${fmt(Math.floor(pos / L.w) + 1)} of ${fmt(v.totalRows(L))} · ${fmtBp(L.w)} per row</dd>
      </dl>
      ${peek ? `<p class="seq-peek">${peek}</p>` : ""}`;
  }

  exonText(g, pos) {
    const tx = this.highlight && this.highlight.gene.id === g.id && this.highlight.tx ? this.highlight.tx : g.tx;
    if (!tx) return "";
    const i = tx.exons.findIndex(([a, b]) => pos >= a && pos < b);
    if (i < 0) return pos >= tx.start && pos < tx.end ? " · intron" : "";
    const n = tx.strand > 0 ? i + 1 : tx.exons.length - i;
    const coding = tx.cds && pos >= tx.cds[0] && pos < tx.cds[1];
    return ` · exon ${n}/${tx.exons.length}${tx.cds ? (coding ? " · coding" : " · UTR") : ""}`;
  }

  async selectionInfo(sel, variant = null) {
    const el = $("#selection");
    if (!sel) {
      el.innerHTML = "";
      return;
    }
    const c = this.view.state.contig;
    const [a, b] = sel;
    const token = ++this.selToken;
    this.refreshIdeo();
    const len = b - a;
    let codes = this.view.codesFor(a, b);
    if (!codes && len <= 2e6) {
      try {
        codes = await this.genome.bases(c, a, b);
      } catch {
        codes = null;
      }
      if (token !== this.selToken) return; // the selection moved on while the bases came in
    }
    const ens = c === "chrM" ? "MT" : c.replace(/^chr/, "");
    const stats = codes && len >= 10 ? ` · GC ${(100 * gcFraction(codes)).toFixed(1)}%` : "";
    const where = len === 1 ? `${c}:${fmt(a + 1)}` : `${c}:${fmt(a + 1)}-${fmt(b)}`;
    el.innerHTML = `
      <h3>${variant ? escapeHtml(variant.id) : "Selection"}</h3>
      <p><b>${where}</b> · ${fmtBp(len)}${stats}</p>
      ${variant ? `<p class="meta">${escapeHtml(variant.alleles || "")} · ${escapeHtml((variant.consequence || "").replace(/_/g, " "))}${variant.maf ? ` · minor allele ${escapeHtml(variant.minor || "")} at ${variant.maf}` : ""}${variant.clinical ? ` · ${escapeHtml(variant.clinical)}` : ""} · <a target="_blank" rel="noopener" href="https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=${encodeURIComponent(variant.id)}">Ensembl</a> · <a target="_blank" rel="noopener" href="https://www.ncbi.nlm.nih.gov/snp/${encodeURIComponent(variant.id)}">dbSNP</a></p>` : ""}
      ${codes && len <= 600 ? `<p class="seq-peek">${lettersOf(codes, true)}</p>` : ""}
      <div class="actions">
        <button data-sel="zoom">Zoom to it</button>
        ${len <= 2e6 ? `<button data-sel="copy">Copy sequence</button><button data-sel="fasta">FASTA</button>` : ""}
        <button data-sel="clear">Clear</button>
      </div>
      <p class="links"><a target="_blank" rel="noopener" href="https://www.ensembl.org/Homo_sapiens/Location/View?r=${ens}:${a + 1}-${b}">Ensembl</a> · <a target="_blank" rel="noopener" href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${c}:${a + 1}-${b}">UCSC</a> · <a target="_blank" rel="noopener" href="https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000001405.40&chr=${ens}&from=${a + 1}&to=${b}">NCBI</a></p>`;
    this.selection = { contig: c, start: a, end: b, variant: variant && variant.id ? variant.id : null };
    this.showTab("here");
  }

  async selectionAction(act) {
    if (act === "clear") {
      this.selToken++;
      this.view.setSelection(null);
      this.selection = null;
      $("#selection").innerHTML = "";
      return this.refreshIdeo();
    }
    const s = this.selection;
    if (!s) return;
    if (act === "zoom") return this.showContig(s.contig, s.start, s.end);
    const codes = await this.genome.bases(s.contig, s.start, s.end);
    const seq = lettersOf(codes, true);
    if (act === "copy") {
      try {
        await navigator.clipboard.writeText(seq);
        this.toast(`Copied ${fmtBp(seq.length)} of sequence`);
      } catch {
        this.toast("The browser didn't allow copying; use FASTA instead");
      }
    } else if (act === "fasta") {
      // the file is named after what was selected: a variant's ID, or the highlighted gene it lies in
      const g = s.variant ? null : this.geneIn(s.contig, s.start, s.end);
      const who = s.variant ? safeName(s.variant) : g ? `${safeName(g.name)}_${g.id}` : "";
      const head = `${s.contig}:${s.start + 1}-${s.end} GRCh38${s.variant ? ` ${s.variant}` : g ? ` ${g.id} ${g.name}` : ""}`;
      this.download(`${who ? `${who}_` : ""}${s.contig}_${s.start + 1}-${s.end}.fa`, wrapFasta(head, seq));
    }
  }

  // ------------------------------------------------------------------ panels & settings
  bindPanels() {
    document.querySelectorAll(".tabs button").forEach((b) => b.addEventListener("click", () => this.showTab(b.dataset.tab)));
    $("#panel-toggle").addEventListener("click", () => this.openPanel());
    $("#to-genome").addEventListener("click", () => this.showGenome());
    $(".logo").addEventListener("click", (e) => {
      e.preventDefault();
      this.showGenome();
    });
    $("#zoom-in").addEventListener("click", () => this.zoom(1));
    $("#zoom-out").addEventListener("click", () => this.zoom(-1));
    $("#random").addEventListener("click", () => this.randomGene());
    $("#snap").addEventListener("click", () => this.snapshot());
    $("#help-btn").addEventListener("click", () => $("#help").showModal());
    $("#downloads-btn").addEventListener("click", () => this.showDownloads());
    $("#downloads").addEventListener("click", async (e) => {
      const b = e.target.closest("[data-copy]");
      if (!b) return;
      const text = $(`#${b.dataset.copy}`).textContent.trim();
      try {
        await navigator.clipboard.writeText(text);
        this.toast("Copied the command");
      } catch {
        getSelection().selectAllChildren($(`#${b.dataset.copy}`));
        this.toast("Press Ctrl+C (⌘C) to copy");
      }
    });
    $("#downloads-link").addEventListener("click", (e) => {
      e.preventDefault();
      this.showDownloads();
    });
    $("#about-link").addEventListener("click", (e) => {
      e.preventDefault();
      $("#help").showModal();
    });
    $("#here").addEventListener("click", (e) => {
      const a = e.target.closest("[data-gene]");
      if (a) {
        e.preventDefault();
        this.selectGene(a.dataset.gene, this.view.state.contig, { zoom: false });
      }
    });
    $("#selection").addEventListener("click", (e) => {
      const b = e.target.closest("[data-sel]");
      if (b) this.selectionAction(b.dataset.sel);
    });
    // palettes
    $("#palettes").innerHTML = Object.entries(PALETTES)
      .map(([k, p]) => `<button data-pal="${k}" title="${escapeHtml(p.about)}"><span class="swatches">${["A", "C", "G", "T"].map((b) => `<i style="background:rgb(${p.rgb[b].join(",")})"></i>`).join("")}</span>${escapeHtml(p.label)}</button>`)
      .join("");
    $("#palettes").addEventListener("click", (e) => {
      const b = e.target.closest("[data-pal]");
      if (b) this.setPalette(b.dataset.pal);
    });
    document.querySelectorAll("[data-opt]").forEach((cb) =>
      cb.addEventListener("change", () => {
        this.view.setOptions({ [cb.dataset.opt]: cb.checked });
        saveOptions(this.view.opts);
      })
    );
    $("#repeats").addEventListener("click", (e) => {
      const b = e.target.closest("button");
      if (!b) return;
      this.view.setOptions({ repeats: b.dataset.v });
      saveOptions(this.view.opts);
      this.lanes.schedule();
      this.syncOptionInputs();
    });
    $("#wrap").addEventListener("click", (e) => {
      const b = e.target.closest("button");
      if (!b) return;
      this.view.setWrap(b.dataset.v, +$("#wrap-width").value || 171);
      this.syncOptionInputs();
      this.pushState();
    });
    $("#wrap-width").addEventListener("change", () => {
      this.view.setWrap("fixed", +$("#wrap-width").value || 171);
      this.syncOptionInputs();
      this.pushState();
    });
    document.querySelectorAll(".presets button").forEach((b) =>
      b.addEventListener("click", () => {
        $("#wrap-width").value = b.dataset.w;
        this.view.setWrap("fixed", +b.dataset.w);
        this.syncOptionInputs();
        this.pushState();
      })
    );
    document.querySelectorAll("[data-live]").forEach((cb) =>
      cb.addEventListener("change", () => {
        this.live.enable(cb.dataset.live, cb.checked);
        const d = this.view.describe();
        if (cb.checked && d) {
          const lim = LIMITS[cb.dataset.live];
          if (d.end - d.start > lim) this.toast(`Zoom in to ${fmtBp(lim)} or less to see ${cb.dataset.live === "variation" ? "variants" : "regulatory elements"}`);
          this.live.update(d);
        }
      })
    );
    $("#find-period").addEventListener("click", () => this.findPeriod());
    $("#bookmark").addEventListener("click", () => this.addBookmark());
    $("#bookmarks").addEventListener("click", (e) => {
      const b = e.target.closest("[data-bm]");
      if (!b) return;
      e.preventDefault(); // the link is for opening in a new tab; here we route ourselves
      const list = storageGet("gpxil.bookmarks", []);
      const i = +b.dataset.bm;
      if (b.dataset.act === "del") {
        list.splice(i, 1);
        storageSet("gpxil.bookmarks", list);
        this.renderBookmarks();
      } else {
        history.pushState(null, "", list[i].hash);
        this.route(list[i].hash);
      }
    });
    // gene panel
    this.panel.on("change", (h) => {
      this.highlight = h.gene ? h : null;
      this.view.setHighlight(this.highlight);
      this.refreshIdeo();
      this.syncLanes();
      if (!this.applying) this.replaceState();
    });
    this.panel.on("close", () => this.clearGene());
    this.panel.on("zoom", (g) => (this.lanes.active ? this.lanes.reset() : this.showContig(g.contig, g.start, g.end)));
    this.panel.on("file", (name, blob) => {
      this.downloadBlob(name, blob);
      this.toast(`⇩ ${name}`, 2500);
    });
    this.panel.on("warn", (m) => this.toast(m));
  }

  showTab(name) {
    document.querySelectorAll(".tabs button").forEach((b) => {
      b.classList.toggle("on", b.dataset.tab === name);
      b.setAttribute("aria-selected", String(b.dataset.tab === name));
    });
    document.querySelectorAll(".tab").forEach((t) => (t.hidden = t.dataset.tab !== name));
    if (name === "gene" && this.panel.gene) this.panel.renderPixels();
  }

  openPanel(force) {
    const app = $("#app");
    const small = matchMedia("(max-width: 900px)").matches;
    if (small) app.classList.toggle("panel-open", force ?? !app.classList.contains("panel-open"));
    else {
      app.classList.toggle("panel-closed", force === undefined ? !app.classList.contains("panel-closed") : !force);
      requestAnimationFrame(() => {
        this.view.resize();
        if (this.screen === "genome") this.karyo.draw();
      });
    }
  }

  setPalette(name, push = true) {
    if (!hasPalette(name)) name = "classic";
    this.palette = name;
    storageSet("gpxil.palette", name);
    this.view.setOptions({ palette: name });
    this.karyo.palette = name;
    if (this.screen === "genome") this.karyo.draw();
    this.ideo.palette = name;
    this.ideo.draw();
    this.panel.palette = name;
    if (this.panel.gene) this.panel.renderPixels();
    this.lanes.schedule();
    document.querySelectorAll("[data-pal]").forEach((b) => b.classList.toggle("on", b.dataset.pal === name));
    if (push) this.replaceState();
  }

  syncOptionInputs() {
    const o = this.view.opts;
    document.querySelectorAll("[data-opt]").forEach((cb) => (cb.checked = !!o[cb.dataset.opt]));
    document.querySelectorAll("#repeats button").forEach((b) => b.classList.toggle("on", b.dataset.v === o.repeats));
    document.querySelectorAll("#wrap button").forEach((b) => b.classList.toggle("on", b.dataset.v === this.view.state.wrap));
    $("#wrap-width").value = this.view.state.width;
  }

  /** Guess the period of a tandem repeat in view: the shift at which bases most often match. */
  async findPeriod() {
    if (this.screen !== "contig") return this.toast("Open a chromosome first");
    const d = this.view.describe();
    const mid = Math.floor((d.start + d.end) / 2);
    const n = Math.min(24000, d.end - d.start);
    const codes = await this.genome.bases(d.contig, mid - (n >> 1), mid + (n >> 1));
    const maxP = Math.min(4000, Math.floor(codes.length / 3));
    let best = { p: 0, score: 0 };
    const scores = new Float32Array(maxP + 1);
    for (let p = 2; p <= maxP; p++) {
      let same = 0;
      let all = 0;
      for (let i = 0; i + p < codes.length; i += 1) {
        const a = codes[i] & 15;
        const b = codes[i + p] & 15;
        if (a === 15 || b === 15) continue;
        all++;
        if (a === b) same++;
      }
      scores[p] = all ? same / all : 0;
      if (scores[p] > best.score) best = { p, score: scores[p] };
    }
    if (best.score < 0.4) return this.toast(`No strong repeat here (best: ${Math.round(100 * best.score)}% of bases match at a shift of ${best.p})`, 4500);
    // prefer the smallest shift that scores nearly as well (multiples of the period score alike)
    let p = best.p;
    for (let q = 2; q < best.p; q++) if (scores[q] >= best.score * 0.97) {
      p = q;
      break;
    }
    this.view.setWrap("fixed", p);
    this.syncOptionInputs();
    this.pushState();
    this.toast(`Repeat period ≈ ${p} bp: ${Math.round(100 * scores[p])}% of bases match one row down`, 4500);
  }

  zoom(steps) {
    if (this.lanes.active) this.lanes.zoom(steps);
    else if (this.screen === "contig") this.view.zoomAt(steps);
  }

  addBookmark() {
    if (this.screen !== "contig") return this.toast("Open a chromosome first");
    const d = this.view.describe();
    const label = this.highlight ? `${this.highlight.gene.name} · ${d.contig}` : `${d.contig}:${fmt(d.start + 1)}-${fmt(d.end)}`;
    const list = storageGet("gpxil.bookmarks", []);
    list.unshift({ label, hash: this.hash(), when: new Date().toISOString().slice(0, 10) });
    storageSet("gpxil.bookmarks", list.slice(0, 50));
    this.renderBookmarks();
    this.toast(`★ Saved ${label}`);
  }

  renderBookmarks() {
    const list = storageGet("gpxil.bookmarks", []);
    $("#bookmarks").innerHTML = list.length
      ? list.map((b, i) => `<li><a href="${escapeHtml(b.hash)}" data-bm="${i}">${escapeHtml(b.label)}</a><span><small>${b.when}</small> <button class="icon" data-bm="${i}" data-act="del" title="Forget">✕</button></span></li>`).join("")
      : `<li class="empty">Nothing saved yet.</li>`;
  }

  async snapshot() {
    if (this.screen !== "contig") return this.toast("Open a chromosome to take a picture");
    // with a highlighted gene in the picture, its name carries the gene's Ensembl ID
    const d = this.lanes.active ? this.lanes.describe() : this.view.describe();
    const g = this.geneIn(d.contig, d.start, d.end);
    const who = g ? `${safeName(g.name)}_${g.id}_` : "";
    if (this.lanes.active) return this.downloadBlob(`gpxil_${who}transcripts_${d.contig}_${d.start + 1}-${d.end}.png`, await this.lanes.snapshot());
    this.downloadBlob(`gpxil_${who}${d.contig}_${d.start + 1}-${d.end}.png`, await this.view.snapshot());
  }

  /** The highlighted gene, if it overlaps [start, end) of a contig: downloads of that stretch carry its Ensembl ID. */
  geneIn(contig, start, end) {
    const g = this.highlight && this.highlight.gene;
    return g && g.contig === contig && g.start < end && start < g.end ? g : null;
  }

  /** The downloads dialog: the genome as gpxil files, and links to the formats. */
  showDownloads() {
    const list = $("#dl-files");
    if (!list.childElementCount) {
      const mb = (b) => (!b ? "" : b < 1e6 ? `${Math.ceil(b / 1e3)} kB` : `${(b / 1e6).toFixed(b < 1e7 ? 1 : 0)} MB`);
      const files = new Map();
      for (const c of this.manifest.sequence.contigs) {
        const f = files.get(c.file) || { file: c.file, names: [], bytes: c.bytes };
        f.names.push(c.name);
        files.set(c.file, f);
      }
      list.innerHTML = [...files.values()]
        .map((f) => {
          const label = f.names.length > 1 ? `${f.names.length} unplaced scaffolds` : f.names[0];
          return `<li><a href="${new URL(f.file, DATA).href}" download>${escapeHtml(label)}</a> <small>${mb(f.bytes)}</small></li>`;
        })
        .join("");
    }
    $("#downloads").showModal();
  }

  download(name, text) {
    this.downloadBlob(name, new Blob([text], { type: "text/plain" }));
  }

  downloadBlob(name, blob) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => {
      URL.revokeObjectURL(a.href);
      a.remove();
    }, 1000);
  }

  toast(msg, ms = 3200) {
    const t = $("#toast");
    t.textContent = msg;
    t.classList.add("on");
    clearTimeout(this.toastTimer);
    this.toastTimer = setTimeout(() => t.classList.remove("on"), ms);
  }

  tooltip(text, x, y, nearCursor = false) {
    const tip = $("#tooltip");
    if (!text) {
      tip.hidden = true;
      return;
    }
    tip.textContent = text;
    tip.hidden = false;
    if (x === undefined) return;
    const box = nearCursor ? $("#view").getBoundingClientRect() : { left: 0, top: 0 };
    const px = box.left + x;
    const py = box.top + y;
    const w = tip.offsetWidth;
    tip.style.left = `${Math.min(innerWidth - w - 8, px + 16)}px`;
    tip.style.top = `${Math.min(innerHeight - 40, py + 18)}px`;
  }

  // ------------------------------------------------------------------ keyboard
  bindKeys() {
    const konami = ["ArrowUp", "ArrowUp", "ArrowDown", "ArrowDown", "ArrowLeft", "ArrowRight", "ArrowLeft", "ArrowRight", "b", "a"];
    let k = 0;
    addEventListener("keydown", (e) => {
      k = e.key === konami[k] ? k + 1 : e.key === konami[0] ? 1 : 0;
      if (k === konami.length) {
        k = 0;
        this.party();
      }
      if (e.target.matches("input, select, textarea") || e.metaKey || e.ctrlKey || e.altKey) return;
      if ($("#help").open) return;
      // Space and Enter press the focused button or link, as usual
      if ((e.key === " " || e.key === "Enter") && e.target.closest("button, a, summary, [role=tab], [tabindex]:not(canvas)")) return;
      const v = this.view;
      const inView = this.screen === "contig";
      const page = inView ? v.visibleRows() * 0.85 : 0;
      if (this.lanes.active) {
        const ln = this.lanes;
        const laneKeys = {
          ArrowLeft: () => ln.pan(e.shiftKey ? -300 : -60),
          ArrowRight: () => ln.pan(e.shiftKey ? 300 : 60),
          ArrowUp: () => ln.scroll(-40),
          ArrowDown: () => ln.scroll(40),
          PageUp: () => ln.scroll(-ln.css.h * 0.8),
          PageDown: () => ln.scroll(ln.css.h * 0.8),
          " ": () => ln.scroll((e.shiftKey ? -1 : 1) * ln.css.h * 0.8),
          Home: () => ln.scroll(-1e9),
          End: () => ln.scroll(1e9),
        }[e.key];
        if (laneKeys) {
          e.preventDefault();
          return laneKeys();
        }
      }
      const act = {
        "/": () => $("#search").focus(),
        g: () => $("#search").focus(),
        k: () => this.showGenome(),
        r: () => this.randomGene(),
        p: () => this.snapshot(),
        "?": () => $("#help").showModal(),
        "+": () => this.zoom(1),
        "=": () => this.zoom(1),
        "-": () => this.zoom(-1),
        _: () => this.zoom(-1),
        "[": () => inView && this.stepGene(-1),
        "]": () => inView && this.stepGene(1),
        1: () => this.panel.gene && this.panel.setMode("gene"),
        2: () => this.panel.gene && this.panel.setMode("transcript"),
        3: () => this.panel.gene && this.panel.setMode("expanded"),
        ArrowDown: () => inView && v.scrollRows(e.shiftKey ? 10 : Math.max(1, 30 / v.layout().rowH)),
        ArrowUp: () => inView && v.scrollRows(-(e.shiftKey ? 10 : Math.max(1, 30 / v.layout().rowH))),
        ArrowRight: () => inView && v.scrollPx(e.shiftKey ? 200 : 40, 0),
        ArrowLeft: () => inView && v.scrollPx(e.shiftKey ? -200 : -40, 0),
        PageDown: () => inView && v.scrollRows(page),
        PageUp: () => inView && v.scrollRows(-page),
        " ": () => inView && v.scrollRows(e.shiftKey ? -page : page),
        Home: () => inView && v.scrollRows(-1e12),
        End: () => inView && v.scrollRows(1e12),
        Escape: () => {
          if (v.selection) this.selectionAction("clear");
          else if (this.highlight) this.clearGene();
        },
      }[e.key];
      if (act) {
        e.preventDefault();
        act();
        if (inView && NAV_KEYS.has(e.key)) this.keyboardCursor();
      }
    });
  }

  /** Without a pointer, the base in the middle of the view is "here" (keyboard navigation). */
  keyboardCursor() {
    requestAnimationFrame(() => {
      const v = this.view;
      if (this.screen !== "contig" || !v.contig) return;
      const pos = v.centreBase();
      const r = v.rects(pos, pos + 1)[0];
      v.hover = pos;
      v.schedule();
      this.pointer = r ? [r[0] + r[2] / 2, r[1] + r[3] / 2] : null;
      this.hoverInfo(pos, ...(this.pointer || []));
    });
  }

  party() {
    document.body.classList.add("party");
    const names = Object.keys(PALETTES);
    const before = this.palette;
    let i = 0;
    this.toast("★ PIXEL PARTY ★ (↑↑↓↓←→←→BA)", 4000);
    const t = setInterval(() => this.setPalette(names[i++ % names.length], false), 350);
    setTimeout(() => {
      clearInterval(t);
      document.body.classList.remove("party");
      this.setPalette(before, false);
    }, 4200);
  }
}

// ------------------------------------------------------------------ boot
async function boot() {
  const landing = new Landing($("#landing"));
  const loader = new DnaLoader($("#loader"));
  void loader;
  const deepLink = location.hash && location.hash !== "#genome";
  const skip = storageGet("gpxil.skipIntro", false);
  const data = (async () => {
    const manifest = await fetchJson(new URL("manifest.json", DATA), "manifest");
    const [karyotype, cytobands] = await Promise.all([fetchJson(new URL(manifest.karyotype, DATA), "karyotype"), fetchJsonGz(new URL(manifest.cytobands, DATA).href, "cytobands")]);
    return { manifest, karyotype, cytobands };
  })();
  let entered;
  if (skip) entered = Promise.resolve();
  else {
    landing.start();
    entered = new Promise((resolve) => {
      const go = () => {
        removeEventListener("keydown", onKey);
        resolve();
      };
      const onKey = (e) => {
        if (e.key === "Enter" || e.key === " ") go();
      };
      $("#enter").addEventListener("click", go, { once: true });
      addEventListener("keydown", onKey);
      $("#enter").focus({ preventScroll: true });
    });
  }
  await entered;
  let d;
  try {
    d = await data;
  } catch (e) {
    landing.stop();
    document.body.innerHTML = `<p class="noscript">Couldn't load the genome data (${escapeHtml(e.message)}). The data files are built by the site's GitHub Actions workflow; see the README.</p>`;
    return;
  }
  await landing.exit();
  $("#app").hidden = false;
  $("#app").classList.add("entering");
  const app = new App(d);
  window.gpxil = app; // for the curious, in the console
  app.start();
  if (!deepLink && !skip) app.toast("Click a chromosome, or search for a gene", 4000);
}

boot();
void SCALES;
void clamp;
