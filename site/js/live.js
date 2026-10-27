// Live tracks from the Ensembl REST API, fetched for the part of the genome
// on screen: known variants (dbSNP, with consequences and clinical
// significance) and the Ensembl regulatory build (promoters, enhancers, CTCF
// sites, open chromatin). Drawn by the viewer over the pixels.

import { REST } from "./genes.js";
import { Emitter, debounce, fetchJson } from "./util.js";

// the widest window each track is fetched for (dbSNP has a variant every few bases)
export const LIMITS = { variation: 30000, regulatory: 3000000 };
const LABELS = { variation: "Ensembl variants", regulatory: "Ensembl regulation" };
const SEVERE = new Set(["pathogenic", "likely pathogenic", "risk factor", "drug response"]);

export const REG_COL = {
  promoter: "rgba(255, 62, 165, 0.95)",
  enhancer: "rgba(255, 225, 77, 0.95)",
  CTCF_binding_site: "rgba(32, 227, 255, 0.95)",
  open_chromatin_region: "rgba(160, 255, 160, 0.9)",
  TF_binding_site: "rgba(200, 150, 255, 0.95)",
  emar: "rgba(255, 160, 80, 0.95)",
};

export class LiveTracks extends Emitter {
  constructor() {
    super();
    this.on_ = { variation: false, regulatory: false };
    this.data = { variation: null, regulatory: null }; // {contig, start, end, items}
    this.pending = { variation: null, regulatory: null }; // {contig, start, end, ctrl}
    this.request = debounce((d) => this.fetchAll(d), 450);
  }

  enable(kind, on) {
    this.on_[kind] = on;
    if (!on) {
      this.data[kind] = null;
      this.cancel(kind);
    }
    this.emit("update");
  }

  enabled(kind) {
    return this.on_[kind];
  }

  cancel(kind) {
    const p = this.pending[kind];
    if (p) p.ctrl.abort();
    this.pending[kind] = null;
  }

  /** Called whenever the view moves. */
  update(d) {
    if (d) this.request(d);
  }

  fetchAll(d) {
    for (const kind of Object.keys(LIMITS)) this.fetchKind(kind, d); // in parallel
  }

  async fetchKind(kind, d) {
    if (!this.on_[kind]) return;
    const len = d.end - d.start;
    if (len > LIMITS[kind]) {
      this.cancel(kind);
      if (this.data[kind]) {
        this.data[kind] = null;
        this.emit("update");
      }
      return;
    }
    const covers = (x) => x && x.contig === d.contig && x.start <= d.start && x.end >= d.end;
    if (covers(this.data[kind]) || covers(this.pending[kind])) return;
    // the window, with a margin on each side so small moves need no new request
    const pad = Math.floor(Math.min(len / 2, (LIMITS[kind] - len) / 2));
    const start = Math.max(0, d.start - pad);
    const end = Math.min(d.length, d.end + pad);
    this.cancel(kind); // a request for somewhere else is no longer wanted
    const mine = { contig: d.contig, start, end, ctrl: new AbortController() };
    this.pending[kind] = mine;
    const chrom = d.contig.replace(/^chr/, "").replace(/^M$/, "MT");
    const url = `${REST}/overlap/region/human/${chrom}:${start + 1}-${end}?feature=${kind}&content-type=application/json`;
    try {
      const rows = await fetchJson(url, LABELS[kind], { headers: { Accept: "application/json" }, signal: mine.ctrl.signal });
      if (this.pending[kind] !== mine || !this.on_[kind]) return;
      this.data[kind] = {
        contig: d.contig,
        start,
        end,
        items: rows.map((r) =>
          kind === "variation"
            ? {
                id: r.id,
                start: r.start - 1,
                end: Math.max(r.end, r.start),
                alleles: (r.alleles || []).join("/"),
                consequence: r.consequence_type || "",
                clinical: (r.clinical_significance || []).join(", "),
                severe: (r.clinical_significance || []).some((c) => SEVERE.has(c)),
              }
            : // the element type (promoter, enhancer, ...) is in "description"; "feature_type" is just "regulatory"
              { id: r.id, start: r.start - 1, end: r.end, type: r.description || r.feature_type || "regulatory" }
        ),
      };
      this.emit("update");
    } catch (e) {
      if (e.name === "AbortError" || this.pending[kind] !== mine) return;
      this.emit("error", `Ensembl ${kind} lookup failed: ${e.message}`);
    } finally {
      if (this.pending[kind] === mine) this.pending[kind] = null;
    }
  }

  /** Items of a track overlapping [start, end) on a contig. */
  items(kind, contig, start, end) {
    const d = this.data[kind];
    if (!this.on_[kind] || !d || d.contig !== contig) return [];
    return d.items.filter((x) => x.end > start && x.start < end);
  }
}
