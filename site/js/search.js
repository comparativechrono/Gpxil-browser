// Search: coordinates, cytobands, chromosomes, genes (symbols, aliases,
// Ensembl IDs) and variants (rsIDs, via Ensembl).

import { Emitter, debounce, escapeHtml, fmtBp } from "./util.js";

const NUM = String.raw`(\d[\d,_]*(?:\.\d+)?)\s*(k|kb|m|mb|g|gb)?`;
const LOCUS = new RegExp(String.raw`^([\w.]+?)\s*[:\s]\s*${NUM}\s*(?:-|–|\.\.|\s)\s*${NUM}$`, "i");
const POINT = new RegExp(String.raw`^([\w.]+?)\s*[:\s]\s*${NUM}$`, "i");
const BAND = /^(?:chr)?(\d{1,2}|x|y)([pq])(\d+(?:\.\d+)?)?$/i;

function num(text, unit) {
  let v = parseFloat(text.replace(/[,_]/g, ""));
  const u = (unit || "").toLowerCase();
  if (u.startsWith("k")) v *= 1e3;
  else if (u.startsWith("m")) v *= 1e6;
  else if (u.startsWith("g")) v *= 1e9;
  return Math.round(v);
}

/**
 * Parse a query that names a place directly. Returns null for anything that
 * needs the gene index or Ensembl (symbols, IDs, rsIDs).
 */
export function parsePlace(q, genome, cytobands) {
  q = q.trim();
  let m = q.match(LOCUS);
  if (m) {
    const c = genome.contig(m[1]);
    if (c) {
      let a = num(m[2], m[3]);
      let b = num(m[4], m[5]);
      if (a > b) [a, b] = [b, a];
      a = Math.min(c.length, Math.max(1, a));
      b = Math.min(c.length, Math.max(1, b));
      return { kind: "locus", contig: c.name, start: a - 1, end: Math.max(a, b), label: `${c.name}:${a.toLocaleString("en-US")}-${b.toLocaleString("en-US")}` };
    }
  }
  m = q.match(POINT);
  if (m) {
    const c = genome.contig(m[1]);
    if (c) {
      const p = Math.min(c.length, Math.max(1, num(m[2], m[3])));
      return { kind: "point", contig: c.name, start: p - 1, end: p, label: `${c.name}:${p.toLocaleString("en-US")}` };
    }
  }
  const c = genome.contig(q.replace(/^chromosome\s*/i, ""));
  if (c) return { kind: "contig", contig: c.name, start: 0, end: c.length, label: `${c.name} (${fmtBp(c.length)})` };
  m = q.match(BAND);
  if (m) {
    const contig = genome.resolve(m[1]);
    const bands = contig && cytobands[contig];
    if (bands) {
      const prefix = `${m[2].toLowerCase()}${m[3] || ""}`;
      const use = bands.filter((b) => b[2].startsWith(prefix)); // 7q3 -> q31.1 ... q36.3; 7q31 -> q31.1, q31.2, q31.3
      if (use.length) {
        const start = Math.min(...use.map((b) => b[0]));
        const end = Math.max(...use.map((b) => b[1]));
        return { kind: "band", contig, start, end, label: `${contig.replace("chr", "")}${prefix} · ${fmtBp(end - start)}` };
      }
    }
  }
  return null;
}

export class Search extends Emitter {
  constructor({ input, list, genome, genes, cytobands }) {
    super();
    this.input = input;
    this.list = list;
    this.genome = genome;
    this.genes = genes;
    this.cytobands = cytobands;
    this.items = [];
    this.active = -1;
    this.seq = 0;
    const update = debounce(() => this.update(), 110);
    input.addEventListener("input", update);
    input.addEventListener("focus", () => {
      this.genes.index().catch(() => {}); // warm the index
      this.genes.aliases().catch(() => {});
      if (input.value.trim()) update();
    });
    input.addEventListener("keydown", (e) => this.key(e));
    input.addEventListener("blur", () => setTimeout(() => this.close(), 150));
    list.addEventListener("pointerdown", (e) => {
      const li = e.target.closest("li[data-i]");
      if (li) {
        e.preventDefault();
        this.choose(this.items[+li.dataset.i]);
      }
    });
  }

  async update() {
    const q = this.input.value.trim();
    const seq = ++this.seq;
    if (!q) return this.close();
    const items = [];
    const place = parsePlace(q, this.genome, this.cytobands);
    if (place) items.push({ ...place, title: place.label, sub: place.kind === "band" ? "cytoband" : place.kind === "contig" ? "whole chromosome" : "go to position" });
    if (/^rs\d+$/i.test(q)) items.push({ kind: "variant", id: q.toLowerCase(), title: q.toLowerCase(), sub: "variant · look up in Ensembl" });
    if (!place || place.kind === "contig") {
      try {
        if (!this.indexReady) this.render([...items, { kind: "wait", title: "loading the gene index…", sub: "" }], q);
        const sug = await this.genes.suggest(q, 10);
        this.indexReady = true;
        if (seq !== this.seq) return;
        for (const s of sug)
          items.push({
            kind: "gene",
            gene: s,
            title: s.name,
            sub: `${s.via ? `alias ${escapeHtml(s.via)} · ` : ""}${s.biotype.replace(/_/g, " ")} · ${s.contig}:${(s.start + 1).toLocaleString("en-US")}`,
          });
      } catch (e) {
        console.warn(e);
      }
      if (!items.some((i) => i.kind === "gene") && !/^rs\d+$/i.test(q) && q.length >= 2)
        items.push({ kind: "lookup", query: q, title: `Search Ensembl for “${escapeHtml(q)}”`, sub: "synonyms, transcript IDs, RefSeq IDs" });
    }
    if (seq !== this.seq) return;
    this.render(items, q);
  }

  render(items, query) {
    this.items = items;
    this.itemsFor = query; // the text these suggestions are for
    this.active = items.length ? 0 : -1;
    if (!items.length) return this.close();
    this.list.innerHTML = items
      .map(
        (it, i) =>
          `<li data-i="${i}" class="${i === this.active ? "active" : ""} ${it.kind}" role="option"><b>${escapeHtml(it.title)}</b><span>${it.sub || ""}</span></li>`
      )
      .join("");
    this.list.hidden = false;
    this.input.setAttribute("aria-expanded", "true");
  }

  close() {
    this.list.hidden = true;
    this.input.setAttribute("aria-expanded", "false");
  }

  key(e) {
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      if (!this.items.length) return;
      e.preventDefault();
      this.active = (this.active + (e.key === "ArrowDown" ? 1 : -1) + this.items.length) % this.items.length;
      [...this.list.children].forEach((li, i) => li.classList.toggle("active", i === this.active));
    } else if (e.key === "Enter") {
      e.preventDefault();
      this.submit();
    } else if (e.key === "Escape") {
      this.close();
      this.input.blur();
    }
  }

  async submit() {
    const q = this.input.value.trim();
    if (!q) return;
    // suggestions count only if they are for what is in the box now (not a moment ago)
    const fresh = this.itemsFor === q && !this.list.hidden && this.active >= 0 && this.items.length;
    if (fresh && this.items[this.active].kind !== "wait") return this.choose(this.items[this.active]);
    const place = parsePlace(q, this.genome, this.cytobands);
    if (place) return this.choose(place);
    if (/^rs\d+$/i.test(q)) return this.choose({ kind: "variant", id: q.toLowerCase() });
    return this.choose({ kind: "lookup", query: q });
  }

  async choose(it) {
    this.close();
    this.input.blur();
    if (!it || it.kind === "wait") return;
    if (it.kind === "gene") return this.emit("gene", it.gene.id, it.gene.contig);
    if (it.kind === "variant") return this.emit("variant", it.id);
    if (it.kind === "lookup") {
      try {
        const r = await this.genes.resolve(it.query);
        if (!r) return this.emit("notfound", it.query);
        return this.emit("gene", r.row[1], r.row[2], r.transcript);
      } catch (e) {
        return this.emit("notfound", it.query, e);
      }
    }
    this.emit("place", it);
  }
}
