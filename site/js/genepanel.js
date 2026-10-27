// The gene panel: what the highlighted gene is, its transcripts, the spliced
// mRNA and the protein as pixels (one pixel per base or amino acid), and
// downloads of the gene, transcripts and proteins as pxil files (pictures of
// sequences: https://github.com/comparativechrono/pxil).

import { AMINO, buildLut } from "./palettes.js";
import { bandAt } from "./karyotype.js";
import { cdsRange, markSelenocysteine, spliced, translate } from "./seq.js";
import { DNA, PROTEIN, indicesOf, pxilFile, safeName, zip } from "./pxil.js";
import { Emitter, escapeHtml, fmt, fmtBp } from "./util.js";

const COMP4 = new Uint8Array(32); // complement of a gpxil code, keeping the soft-mask bit
for (let c = 0; c < 32; c++) COMP4[c] = (c & 16) | ((c & 1) << 3) | ((c & 2) << 1) | ((c & 4) >> 1) | ((c & 8) >> 3);

const MODES = [
  ["gene", "Gene", "The whole gene lit up in the genome."],
  ["transcript", "Transcript", "One transcript lit up in the genome: exons bright, UTRs half-lit, introns dimmed. Pick another below."],
  ["expanded", "Expanded", "Every transcript on its own lane in the viewer. Click a lane to open that transcript."],
];

export class GenePanel extends Emitter {
  constructor({ root, genome, genes, cytobands }) {
    super();
    this.root = root;
    this.genome = genome;
    this.genes = genes;
    this.cytobands = cytobands;
    this.gene = null;
    this.tx = null;
    this.mode = "gene";
    this.txs = null;
    this.codes = null;
    this.palette = "classic";
    this.token = 0; // answers that arrive after another gene was opened are dropped
    this.picks = { gene: false, transcript: new Set(), protein: new Set() }; // for the zip of pxil files
    this.manyOpen = false;
    this.source = "GRCh38"; // set by the app from the data's manifest
    root.addEventListener("click", (e) => this.click(e));
    root.addEventListener("change", (e) => {
      if (e.target.matches("select[data-role=tx]")) this.pickTx(e.target.value);
      else if (e.target.matches("[data-pick], [data-pick-all]")) this.pick(e.target);
    });
    root.addEventListener("toggle", (e) => e.target.matches("details.pxil-many") && (this.manyOpen = e.target.open), true);
  }

  clear() {
    this.gene = null;
    this.tx = null;
    this.root.innerHTML = `<p class="empty">Search for a gene, or click one in the viewer, and it lights up here.</p>`;
  }

  async open(gene, { mode = this.mode, transcript = null } = {}) {
    const token = ++this.token;
    const same = this.gene && this.gene.id === gene.id;
    this.gene = gene;
    this.mode = mode;
    if (!same) {
      this.txs = null;
      this.codes = null;
      this.summary = undefined;
      this.picks = { gene: false, transcript: new Set(), protein: new Set() };
    }
    this.tx = gene.tx;
    if (transcript && (!this.tx || this.tx.id !== transcript)) {
      await this.loadTranscripts();
      if (token !== this.token) return;
      this.tx = (this.txs || []).find((t) => t.id === transcript) || this.tx;
    }
    this.render();
    if (!this.txs)
      this.loadTranscripts().then(() => {
        if (token !== this.token) return;
        this.render();
        this.changed();
      });
    this.loadCodes()
      .then(() => token === this.token && this.renderPixels())
      .catch((e) => this.emit("warn", `Couldn't read ${gene.name}'s bases: ${e.message}`));
    if (this.summary === undefined) {
      this.summary = null;
      this.genes
        .summary(gene)
        .then((s) => {
          if (this.gene !== gene) return;
          this.summary = s;
          const el = this.root.querySelector(".summary");
          if (el && s) el.innerHTML = `${escapeHtml(s)} <small>— NCBI Gene, via MyGene.info</small>`;
        })
        .catch(() => {});
    }
    this.changed();
  }

  changed() {
    this.emit("change", { gene: this.gene, tx: this.tx, mode: this.mode, txs: this.txs });
  }

  async loadTranscripts() {
    if (this.txs) return this.txs;
    const g = this.gene;
    let txs;
    try {
      const d = await this.genes.transcripts(g);
      if (!g.description && d.description) g.description = d.description.replace(/\s*\[Source:.*\]$/, "");
      txs = d.txs;
    } catch (e) {
      txs = g.tx ? [g.tx] : [];
      if (this.gene === g) this.emit("warn", "Ensembl is not answering, so only the canonical transcript is shown.");
    }
    if (this.gene !== g) return null; // another gene was opened meanwhile
    this.txs = txs;
    return txs;
  }

  async loadCodes() {
    if (this.codes && this.codes.id === this.gene.id) return;
    const g = this.gene;
    const codes = await this.genome.bases(g.contig, g.start, g.end);
    if (this.gene === g) this.codes = { id: g.id, start: g.start, codes };
  }

  pickTx(id) {
    const t = (this.txs || []).find((x) => x.id === id) || (this.gene.tx && this.gene.tx.id === id ? this.gene.tx : null);
    if (!t) return;
    this.tx = t;
    if (this.mode === "gene") this.mode = "transcript";
    this.render();
    this.renderPixels();
    this.changed();
  }

  setMode(m) {
    this.mode = m;
    const g = this.gene;
    if (!this.txs)
      this.loadTranscripts().then(() => {
        if (this.gene !== g) return;
        this.render();
        this.changed();
      });
    this.render();
    this.renderPixels();
    this.changed();
  }

  render() {
    const g = this.gene;
    if (!g) return this.clear();
    const band = bandAt(this.cytobands[g.contig], g.start);
    const tx = this.tx;
    const txs = this.txs || (g.tx ? [g.tx] : []);
    const tag = (t) => (t.mane ? "MANE Select" : t.canonical ? "canonical" : "");
    const exonLen = tx ? tx.exons.reduce((s, [a, b]) => s + b - a, 0) : 0;
    const cdsLen = tx && tx.cds ? tx.exons.reduce((s, [a, b]) => s + Math.max(0, Math.min(b, tx.cds[1]) - Math.max(a, tx.cds[0])), 0) : 0;
    const x = g.xref || {};
    const links = [
      ["Ensembl", `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${g.id}`],
      x.entrez && ["NCBI Gene", `https://www.ncbi.nlm.nih.gov/gene/${x.entrez}`],
      x.hgnc && ["HGNC", `https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:${x.hgnc}`],
      x.uniprot && ["UniProt", `https://www.uniprot.org/uniprotkb/${x.uniprot}`],
      x.omim && ["OMIM", `https://omim.org/entry/${x.omim}`],
      ["GeneCards", `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(g.name)}`],
      ["UCSC", `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${g.contig}:${g.start + 1}-${g.end}`],
    ].filter(Boolean);
    this.root.innerHTML = `
      <div class="gene-head">
        <h2>${escapeHtml(g.name)} <span class="strand" title="${g.strand > 0 ? "forward" : "reverse"} strand">${g.strand > 0 ? "▶" : "◀"}</span></h2>
        <button class="icon close" data-act="close" title="Clear the highlight (Esc)">✕</button>
      </div>
      <p class="desc">${escapeHtml(g.description || "")}</p>
      <p class="meta">${escapeHtml(g.biotype.replace(/_/g, " "))} · <a href="#" data-act="zoom">${g.contig}:${fmt(g.start + 1)}-${fmt(g.end)}</a> · ${fmtBp(g.end - g.start)}${band ? ` · ${g.contig.replace("chr", "")}${band[2]}` : ""}</p>
      <div class="modes" role="tablist">${MODES.map(([k, l, about]) => `<button role="tab" data-mode="${k}" aria-selected="${this.mode === k}" title="${about}" class="${this.mode === k ? "on" : ""}">${l}</button>`).join("")}</div>
      <p class="hint mode-hint">${MODES.find(([k]) => k === this.mode)[2]}</p>
      ${
        txs.length
          ? `<label class="txpick">Transcript${this.txs ? ` <small>${fmt(this.txs.length)} in all</small>` : ""} <select data-role="tx">${txs
              .map((t) => `<option value="${t.id}" ${tx && t.id === tx.id ? "selected" : ""}>${escapeHtml(t.name)} · ${t.id}${tag(t) ? ` · ${tag(t)}` : ""}</option>`)
              .join("")}</select></label>`
          : ""
      }
      ${!this.txs && g.txCount > 1 ? `<p class="hint">loading all ${fmt(g.txCount)} transcripts…</p>` : ""}
      ${
        tx
          ? `<p class="meta">${escapeHtml(tx.biotype.replace(/_/g, " "))} · ${tx.exons.length} exons · ${fmt(exonLen)} nt${cdsLen ? ` · CDS ${fmt(cdsLen)} nt · ${fmt(cdsLen / 3 - 1)} aa` : ""}${tx.incomplete ? ` · ${escapeHtml(tx.incomplete)}` : ""}</p>`
          : ""
      }
      <div class="pixels">
        ${this.mode === "expanded" ? `<p class="hint">In the lanes, coding exons fill the lane, UTRs and non-coding exons half of it, and introns are a thin line with arrows for the strand; the top lane is the DNA. Hover for exon numbers. Zoom with + and −, or ctrl + scroll.</p>` : ""}
        <h3>mRNA <small>${tx ? `${escapeHtml(tx.name)}, ` : ""}exons spliced, 5′→3′</small></h3><canvas class="mrna"></canvas>
        ${tx && tx.cds ? `<h3>Protein <small>one pixel per amino acid</small></h3><canvas class="protein"></canvas>` : ""}
        <p class="pixel-info" aria-live="polite"></p>
      </div>
      ${this.downloadsHtml(txs)}
      <p class="summary">${this.summary ? `${escapeHtml(this.summary)} <small>— NCBI Gene, via MyGene.info</small>` : ""}</p>
      <div class="actions"><button data-act="zoom">Zoom to gene</button></div>
      <p class="links">${links.map(([l, u]) => `<a href="${u}" target="_blank" rel="noopener">${l}</a>`).join(" · ")}</p>`;
    this.renderPixels();
  }

  /** The pxil downloads: this transcript, its protein and the gene, or several at once (zipped). */
  downloadsHtml(txs) {
    const g = this.gene;
    const tx = this.tx;
    const p = this.picks;
    const names = this.fileNames(tx);
    const rows = txs
      .map(
        (t) => `<tr><td>${escapeHtml(t.name)}${t.mane ? " ★" : ""}<small>${escapeHtml(t.biotype.replace(/_/g, " "))}</small></td>
          <td><input type="checkbox" data-pick="transcript" value="${t.id}" aria-label="${escapeHtml(t.name)} transcript" ${p.transcript.has(t.id) ? "checked" : ""}></td>
          <td>${t.cds ? `<input type="checkbox" data-pick="protein" value="${t.id}" aria-label="${escapeHtml(t.name)} protein" ${p.protein.has(t.id) ? "checked" : ""}>` : "–"}</td></tr>`
      )
      .join("");
    return `
      <h3>Download as pixels <small><a href="https://github.com/comparativechrono/pxil" target="_blank" rel="noopener">pxil files</a></small></h3>
      <div class="actions pxil">
        ${tx ? `<button data-act="pxil-transcript" title="${escapeHtml(names.transcript)}">transcript.pxil</button>` : ""}
        ${tx && tx.cds ? `<button data-act="pxil-protein" title="${escapeHtml(names.protein)}">protein.pxil</button>` : ""}
        <button data-act="pxil-gene" title="${escapeHtml(names.gene)}">gene.pxil</button>
      </div>
      <details class="pxil-many" ${this.manyOpen ? "open" : ""}>
        <summary>Several transcripts and proteins at once</summary>
        <div class="pxil-table"><table>
          <thead><tr><th>${fmt(txs.length)} transcript${txs.length === 1 ? "" : "s"}</th>
            <th><label><input type="checkbox" data-pick-all="transcript"> transcript</label></th>
            <th><label><input type="checkbox" data-pick-all="protein"> protein</label></th></tr></thead>
          <tbody>${rows}</tbody>
        </table></div>
        <label class="pick-gene"><input type="checkbox" data-pick="gene" ${p.gene ? "checked" : ""}> and ${escapeHtml(g.name)}'s genomic sequence (gene.pxil)</label>
        <button data-act="pxil-zip" class="wide">${this.zipLabel()}</button>
      </details>
      <p class="hint">The pictures above, one pixel per base or amino acid, as PNG images whose pixel values are the sequence: they open in image software and read back exactly. File names carry the Ensembl IDs.</p>`;
  }

  zipLabel() {
    const n = (this.picks.gene ? 1 : 0) + this.picks.transcript.size + this.picks.protein.size;
    return n ? `Download ${n} file${n === 1 ? "" : "s"} (.zip)` : "Tick the files you want";
  }

  pick(input) {
    const p = this.picks;
    const txs = this.txs || (this.gene.tx ? [this.gene.tx] : []);
    if (input.dataset.pickAll) {
      const kind = input.dataset.pickAll;
      for (const t of txs) if (kind === "transcript" || t.cds) input.checked ? p[kind].add(t.id) : p[kind].delete(t.id);
      this.root.querySelectorAll(`[data-pick="${kind}"]`).forEach((cb) => (cb.checked = input.checked));
    } else if (input.dataset.pick === "gene") p.gene = input.checked;
    else input.checked ? p[input.dataset.pick].add(input.value) : p[input.dataset.pick].delete(input.value);
    const b = this.root.querySelector("[data-act=pxil-zip]");
    if (b) b.textContent = this.zipLabel();
  }

  /** File names, with the Ensembl IDs: gene ENSG, transcript ENST, protein ENSP. */
  fileNames(t) {
    const g = this.gene;
    return {
      gene: `${safeName(g.name)}_${g.id}.gene.pxil`,
      transcript: t ? `${safeName(t.name)}_${t.id}.transcript.pxil` : "",
      protein: t ? `${safeName(t.name)}_${t.protein || `${t.id}_protein`}.protein.pxil` : "",
    };
  }

  panelWidth() {
    const w = this.root.getBoundingClientRect().width;
    return w ? Math.max(160, Math.floor(w - 8)) : 340;
  }

  common(palette = this.palette) {
    const g = this.gene;
    return { source: this.source, assembly: "GRCh38", contig: g.contig, created: new Date().toISOString().replace(/\.\d+Z$/, "Z"), generator: "Gpxil browser", palette };
  }

  dnaColours() {
    const lut = buildLut(this.palette);
    return Array.from({ length: 32 }, (_, c) => [lut.rgb[3 * c], lut.rgb[3 * c + 1], lut.rgb[3 * c + 2]]);
  }

  /** The gene's genomic sequence, 5′ to 3′ of the gene, soft-masked repeats kept. */
  async genePxil() {
    await this.loadCodes();
    const g = this.gene;
    const n = this.codes.codes.length;
    const idx = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      const c = this.codes.codes[g.strand > 0 ? i : n - 1 - i] & 31;
      idx[i] = g.strand > 0 ? c : COMP4[c];
    }
    const seq = Array.from(idx, (i) => DNA[i]).join("");
    const meta = {
      kind: "gene",
      id: g.id,
      name: g.name,
      description: g.description || undefined,
      biotype: g.biotype,
      ...this.common(),
      start: g.start,
      end: g.end,
      strand: g.strand,
      orientation: g.strand > 0 ? "5' to 3' of the gene: the reference's forward strand" : "5' to 3' of the gene: the reverse complement of the reference",
    };
    const bytes = await pxilFile({ indices: idx, alphabet: DNA, width: pickWrap(n, this.panelWidth(), 4).cols, colours: this.dnaColours(), meta, sequence: seq });
    return { name: this.fileNames(null).gene, bytes };
  }

  /** A transcript's spliced sequence and, for coding ones, its protein. */
  async spliceOf(t) {
    await this.loadCodes();
    const { seq, map } = spliced(this.codes.codes, this.codes.start, t);
    const range = cdsRange(t, map);
    let protein = null;
    if (range) {
      const start = range[0] + (t.frame || 0);
      protein = markSelenocysteine(translate(seq.slice(start, range[1]), this.gene.contig === "chrM"), map, start, t.sec);
    }
    return { seq, range, protein };
  }

  async transcriptPxil(t) {
    const g = this.gene;
    const { seq, range } = await this.spliceOf(t);
    const meta = {
      kind: "transcript",
      id: t.id,
      name: t.name,
      biotype: t.biotype,
      gene: g.id,
      gene_name: g.name,
      protein: t.protein || undefined,
      tags: [t.mane && "MANE Select", t.maneClinical && "MANE Plus Clinical", t.canonical && "Ensembl canonical"].filter(Boolean),
      ...this.common(),
      start: t.start,
      end: t.end,
      strand: t.strand,
      orientation: "5' to 3', exons spliced",
      exons: t.exons,
      cds: range ? [range[0], range[1]] : undefined,
    };
    const bytes = await pxilFile({ indices: indicesOf(seq, DNA), alphabet: DNA, width: pickWrap(seq.length, this.panelWidth(), 4).cols, colours: this.dnaColours(), meta, sequence: seq });
    return { name: this.fileNames(t).transcript, bytes };
  }

  async proteinPxil(t) {
    const g = this.gene;
    const { protein } = await this.spliceOf(t);
    if (!protein) return null;
    const meta = {
      kind: "protein",
      id: t.protein || undefined,
      name: t.name,
      transcript: t.id,
      gene: g.id,
      gene_name: g.name,
      genetic_code: g.contig === "chrM" ? "vertebrate mitochondrial" : "standard",
      description: "translated in the Gpxil browser from the transcript's coding sequence in GRCh38; * is the stop codon",
      ...this.common("Taylor (1997) amino acid colours"),
    };
    const colours = [...PROTEIN].map((ch, i) => (i === 0 ? [7, 7, 15] : AMINO[ch] || AMINO.X));
    const bytes = await pxilFile({ indices: indicesOf(protein, PROTEIN), alphabet: PROTEIN, width: pickWrap(protein.length, this.panelWidth(), 6).cols, colours, meta, sequence: protein });
    return { name: this.fileNames(t).protein, bytes };
  }

  async downloadPxil(act) {
    const g = this.gene;
    try {
      if (act === "pxil-gene") return this.emit("file", ...(await this.fileOf(this.genePxil())));
      if (act === "pxil-transcript" && this.tx) return this.emit("file", ...(await this.fileOf(this.transcriptPxil(this.tx))));
      if (act === "pxil-protein" && this.tx) return this.emit("file", ...(await this.fileOf(this.proteinPxil(this.tx))));
      if (act === "pxil-zip") {
        const p = this.picks;
        const txs = this.txs || (g.tx ? [g.tx] : []);
        const files = [];
        if (p.gene) files.push(await this.genePxil());
        for (const t of txs) {
          if (p.transcript.has(t.id)) files.push(await this.transcriptPxil(t));
          if (p.protein.has(t.id)) files.push(await this.proteinPxil(t));
        }
        const ok = files.filter(Boolean);
        if (!ok.length) return this.emit("warn", "Tick the transcripts and proteins you want first");
        if (ok.length === 1) return this.emit("file", ok[0].name, new Blob([ok[0].bytes], { type: "image/png" }));
        this.emit("file", `${safeName(g.name)}_${g.id}.pxil.zip`, new Blob([zip(ok)], { type: "application/zip" }));
      }
    } catch (e) {
      this.emit("warn", `Couldn't make the file: ${e.message}`);
    }
  }

  async fileOf(promise) {
    const f = await promise;
    return [f.name, new Blob([f.bytes], { type: "image/png" })];
  }

  renderPixels() {
    if (!this.gene || !this.codes || this.codes.id !== this.gene.id) return;
    const lut = buildLut(this.palette);
    const width = this.panelWidth();
    const info = this.root.querySelector(".pixel-info");
    const tx = this.tx;
    const mrna = this.root.querySelector("canvas.mrna");
    if (!tx || !mrna) return;
    const { seq, map } = spliced(this.codes.codes, this.codes.start, tx);
    const range = cdsRange(tx, map);
    this.mrna = { seq, map, range };
    const wrap = pickWrap(seq.length, width, 4);
    drawLetters(mrna, seq, wrap, lut, (i) => (range && (i < range[0] || i >= range[1]) ? 0.45 : 1));
    mrna.onpointermove = (e) => {
      const i = cellAt(mrna, e, wrap);
      if (i === null || i >= seq.length) return;
      const frame = tx.frame || 0; // bases of an incomplete first codon
      const where = range ? (i < range[0] ? "5′ UTR" : i >= range[1] ? "3′ UTR" : i < range[0] + frame ? "incomplete first codon" : `codon ${Math.floor((i - range[0] - frame) / 3) + 1}`) : "non-coding";
      info.textContent = `mRNA ${fmt(i + 1)} · ${seq[i]} · ${this.gene.contig}:${fmt(map[i] + 1)} · ${where}`;
    };
    const prot = this.root.querySelector("canvas.protein");
    if (prot && range) {
      const mito = this.gene.contig === "chrM";
      const start = range[0] + (tx.frame || 0);
      const aa = markSelenocysteine(translate(seq.slice(start, range[1]), mito), map, start, tx.sec);
      this.protein = aa;
      const pw = pickWrap(aa.length, width, 6);
      drawAmino(prot, aa, pw);
      prot.onpointermove = (e) => {
        const i = cellAt(prot, e, pw);
        if (i === null || i >= aa.length) return;
        const codon = seq.slice(start + 3 * i, start + 3 * i + 3);
        const name = aa[i] === "*" ? "stop" : aa[i] === "U" ? "U, selenocysteine" : aa[i];
        info.textContent = `residue ${fmt(i + 1)} · ${name} (${codon})${mito ? " · mitochondrial code" : ""}`;
      };
    }
  }

  async click(e) {
    const b = e.target.closest("[data-act],[data-mode]");
    if (!b) return;
    if (b.dataset.mode) return this.setMode(b.dataset.mode);
    e.preventDefault();
    const g = this.gene;
    const act = b.dataset.act;
    if (act === "close") return this.emit("close");
    if (act === "zoom") return this.emit("zoom", g);
    if (act.startsWith("pxil")) return this.downloadPxil(act);
  }
}

function pickWrap(n, width, maxPx) {
  // the largest pixel size (<= maxPx) that keeps the raster under ~360 px tall
  for (let px = maxPx; px >= 1; px--) {
    const cols = Math.floor(width / px);
    if (Math.ceil(n / cols) * px <= 360 || px === 1) return { px, cols };
  }
  return { px: 1, cols: width };
}

function drawLetters(cv, seq, { px, cols }, lut, alphaAt) {
  const rows = Math.max(1, Math.ceil(seq.length / cols));
  cv.width = cols;
  cv.height = rows;
  cv.style.width = `${cols * px}px`;
  cv.style.height = `${rows * px}px`;
  const ctx = cv.getContext("2d");
  const img = ctx.createImageData(cols, rows);
  const codeOf = { A: 1, C: 2, G: 4, T: 8 };
  for (let i = 0; i < seq.length; i++) {
    const c = (codeOf[seq[i]] || 15) * 3;
    const k = alphaAt(i);
    const o = 4 * i;
    img.data[o] = lut.rgb[c] * k + 7 * (1 - k);
    img.data[o + 1] = lut.rgb[c + 1] * k + 7 * (1 - k);
    img.data[o + 2] = lut.rgb[c + 2] * k + 15 * (1 - k);
    img.data[o + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
}

function drawAmino(cv, aa, { px, cols }) {
  const rows = Math.max(1, Math.ceil(aa.length / cols));
  cv.width = cols;
  cv.height = rows;
  cv.style.width = `${cols * px}px`;
  cv.style.height = `${rows * px}px`;
  const ctx = cv.getContext("2d");
  const img = ctx.createImageData(cols, rows);
  for (let i = 0; i < aa.length; i++) {
    const c = AMINO[aa[i]] || AMINO.X;
    img.data.set([c[0], c[1], c[2], 255], 4 * i);
  }
  ctx.putImageData(img, 0, 0);
}

function cellAt(cv, e, { px, cols }) {
  const r = cv.getBoundingClientRect();
  const x = Math.floor((e.clientX - r.left - cv.clientLeft) / px); // inside the border
  const y = Math.floor((e.clientY - r.top - cv.clientTop) / px);
  if (x < 0 || x >= cols || y < 0) return null;
  return y * cols + x;
}
