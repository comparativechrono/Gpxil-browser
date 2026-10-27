// The gene panel: what the highlighted gene is, its transcripts, and the
// same gene as pixels three more ways: every transcript as a pixel track
// (expanded), the spliced mRNA, and the protein (one pixel per amino acid).

import { AMINO, buildLut, inkOn } from "./palettes.js";
import { bandAt } from "./karyotype.js";
import { cdsRange, gcFraction, markSelenocysteine, spliced, translate, wrapFasta } from "./seq.js";
import { Emitter, escapeHtml, fmt, fmtBp } from "./util.js";

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
    root.addEventListener("click", (e) => this.click(e));
    root.addEventListener("change", (e) => {
      if (e.target.matches("select[data-role=tx]")) this.pickTx(e.target.value);
    });
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
      <p class="summary">${this.summary ? `${escapeHtml(this.summary)} <small>— NCBI Gene, via MyGene.info</small>` : ""}</p>
      <div class="actions">
        <button data-act="zoom">Zoom to gene</button>
        <button data-act="fasta-genomic">Genomic FASTA</button>
        ${tx ? `<button data-act="fasta-mrna">mRNA FASTA</button>` : ""}
        ${tx && tx.cds ? `<button data-act="fasta-protein">Protein FASTA</button>` : ""}
      </div>
      <p class="links">${links.map(([l, u]) => `<a href="${u}" target="_blank" rel="noopener">${l}</a>`).join(" · ")}</p>`;
    this.renderPixels();
  }

  renderPixels() {
    if (!this.gene || !this.codes || this.codes.id !== this.gene.id) return;
    const lut = buildLut(this.palette);
    const width = Math.max(160, Math.floor(this.root.getBoundingClientRect().width - 8));
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
    if (act.startsWith("fasta")) {
      await this.loadCodes();
      let text = "";
      if (act === "fasta-genomic") {
        const { lettersOf } = await import("./seq.js");
        text = wrapFasta(`${g.name} ${g.contig}:${g.start + 1}-${g.end} GRCh38 (forward strand, soft-masked) GC=${(100 * gcFraction(this.codes.codes)).toFixed(1)}%`, lettersOf(this.codes.codes, true));
      } else if (act === "fasta-mrna" && this.mrna) text = wrapFasta(`${this.tx.id} ${this.tx.name} mRNA`, this.mrna.seq);
      else if (act === "fasta-protein" && this.protein) text = wrapFasta(`${this.tx.id} ${this.tx.name} protein (translated here)`, this.protein.replace(/\*$/, ""));
      this.emit("download", `${g.name}-${act.replace("fasta-", "")}.fa`, text);
    }
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
  void inkOn;
}

function cellAt(cv, e, { px, cols }) {
  const r = cv.getBoundingClientRect();
  const x = Math.floor((e.clientX - r.left - cv.clientLeft) / px); // inside the border
  const y = Math.floor((e.clientY - r.top - cv.clientTop) / px);
  if (x < 0 || x >= cols || y < 0) return null;
  return y * cols + x;
}
