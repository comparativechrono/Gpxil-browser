# Gpxil browser: enter the pixelverse

A genome browser for the human reference genome (GRCh38) in which the view is
nothing but pixels: one pixel per base, read straight from
[gpxil](https://github.com/comparativechrono/Gpxil) files, the image format for
genomes that follows Albrecht-Buehler's genome pixel images
([PLoS ONE, 2009](https://doi.org/10.1371/journal.pone.0004701)).

It is a static site for GitHub Pages: **https://comparativechrono.github.io/Gpxil-browser/**
(live once Pages is switched on, see [Deploying](#deploying)).

![The landing screen: "Enter the pixelverse" assembling from flying bases above a turning pixel double helix](docs/img/landing.png)

| | |
|---|---|
| ![The karyotype: every chromosome as a column of pixels, with cytobands and gene density](docs/img/karyotype.png) | ![TP53 highlighted, with its Ensembl transcripts drawn as pixel tracks](docs/img/gene.png) |
| The karyotype, one pixel per 131,072 bases | A gene, its transcripts and its exons |
| ![The chromosome 1 centromere wrapped at 340 bases per row, where the alpha-satellite repeat lines up into columns](docs/img/repeat.png) | |
| Alpha satellite on chr1, wrapped at the repeat period the browser found (340 bp) | |

## What you can do

**Move around the genome**

- **Karyotype**: all 24 chromosomes as columns of pixels, coloured by base
  content, with cytobands on the left and protein-coding gene density on the
  right. chrM, far smaller than one karyotype pixel, is drawn as a ring
  (it is circular). Click anywhere to dive in.
- **Chromosome view**: the bases wrap into rows like lines of text. Zoom runs
  from 262,144 bases per pixel to 64 pixels per base; from 12 pixels per base
  the letters appear. Scroll, drag, pinch, double-click, or use the keyboard.
- **Ideogram** above the view: cytobands, the chromosome's own pixels, and a
  bracket for the part on screen; click or drag it to jump.
- **Gutter tracks** beside the row numbers: GC content of each row, genes and
  exons, and live data when it is switched on.
- **Wrap**: fit the screen, or a fixed number of bases per row, which lines up
  tandem repeats as stripes. *Find the repeat period here* picks the width for
  you by autocorrelation.

**Search** (press `/`)

- Genes by HGNC symbol, alias or previous symbol (HER2 finds ERBB2), or
  Ensembl gene ID; anything else (transcript IDs, other synonyms) is looked
  up in Ensembl.
- Places: `chr7:117,480,025-117,668,665`, `X:73.8M`, `17:7.66M-7.69M`, a
  whole chromosome, or a cytoband (`17q21.31`, `7q3`).
- Variants by rsID (`rs334`), placed through Ensembl.

**Genes**

- The gene's pixels stay lit while the rest of the view dims.
- Three modes (`1` `2` `3`): **gene**, **transcript** (introns dimmed, UTRs
  lightly dimmed, exons boxed) and **expanded** (every transcript as a pixel
  track). The canonical transcript (MANE Select where there is one) comes with
  the site; all the others are fetched from Ensembl and can be picked from a list.
- The spliced mRNA and the protein drawn as pixels (one pixel per amino acid,
  translated here, with the mitochondrial code on chrM and selenocysteine where
  Ensembl marks it); hover for codon and position.
- Genomic, mRNA and protein FASTA downloads; the NCBI Gene summary; links to
  Ensembl, NCBI Gene, HGNC, UniProt and OMIM; previous and next gene (`[` `]`).

**Like a traditional browser, but in pixels**

- Hover any pixel for its base, position, band, gene, exon/intron/UTR and, with
  live data on, variants and regulatory elements. Moving with the keyboard, the
  base in the middle of the view is the one described.
- Shift-drag to select bases: length, GC, the sequence, copy, FASTA, and the
  same place in Ensembl, UCSC or NCBI.
- Live tracks from Ensembl: known variants (dbSNP, with consequences;
  clinically significant ones in red) and the regulatory build (promoters,
  enhancers, CTCF sites, open chromatin).
- CpG sites, base letters, a pixel grid, six palettes (Classic, GPxI greys,
  Okabe-Ito, Process inks, Neon, GC heat), and repeats faded, plain or hidden.
- Every view has its own link (the address bar), saved places, a picture of the
  view (`P`), and a random gene (`R`). There may be other surprises.
- Whenever something takes a while, a pixel double helix spins and says what it
  is waiting for.

## How it works

The site has no build step: plain HTML, CSS and JavaScript modules in
[`site/`](site/). The data is built by [`pipeline/`](pipeline/) into
`site/data/grch38/`:

| Files | What | Size |
|---|---|---|
| `gpxil/chr*.gpxil.tif`, `gpxil/unplaced.gpxil.tif` | The reference as gpxil: one reference-only file per chromosome, 8,192 bases per row, Deflate, 256 KiB strips; the 169 unplaced and unlocalised scaffolds share one file | 849 MB |
| `overview/*.ov.gz` | Summary bins for zoomed-out views: for bins of 128 to 131,072 bases, GC, N and soft-masked fractions | 39 MB |
| `genes/*.json.gz` | 78,941 genes (Ensembl 116) with HGNC names, cross-references and canonical transcripts; a search index; 58,149 aliases | 8 MB |
| `karyotype.json`, `cytobands.json.gz`, `manifest.json` | The karyotype bins, cytobands, and the index of everything | 0.1 MB |

The browser never downloads a whole chromosome. Like a cloud-optimised
GeoTIFF, a gpxil file can be read in pieces: the manifest gives each contig's
TIFF page, and the browser reads the page's strip table and then only the
strips under the view, with HTTP range requests, and inflates them with the
browser's own `DecompressionStream`. Zoomed out, it reads the summary bins
instead.

Some things are not stored on the site but fetched live, only when asked for:

| From | What |
|---|---|
| [Ensembl REST API](https://rest.ensembl.org) | All transcripts of a gene, transcript IDs and other gene names, rsIDs, variants and regulatory elements in view |
| [MyGene.info](https://mygene.info) | Gene summaries from NCBI Gene |

## Building the data

Python 3.10 or later:

```sh
pip install -r pipeline/requirements.txt
python -I pipeline/build.py --out site/data/grch38
```

This downloads the GRCh38 soft-masked primary assembly and the gene set of
Ensembl release 116, the HGNC complete set and the cytobands (about 1 GB, kept
in `.cache/`), and takes about ten minutes. Any input can be given as a local
file instead (`--fasta`, `--gtf`, `--hgnc`, `--bands`); `--only chr21,chrM`
makes a quick partial build.

Then serve the site with range requests (Python's own `http.server` ignores
them, which makes everything slow):

```sh
python dev/serve.py        # http://localhost:8000
```

## Deploying

[`.github/workflows/pages.yml`](.github/workflows/pages.yml) builds the data,
adds the site and deploys it to GitHub Pages on every push to `main`. The data
is cached between runs and rebuilt only when `pipeline/` changes, when the
cache has expired, or when the workflow is run by hand with *rebuild* ticked.

One-time setup: in the repository's **Settings → Pages → Build and
deployment**, set **Source** to **GitHub Actions**.

The whole site is about 900 MB, inside GitHub Pages' 1 GB limit for a site.

## Tests

```sh
pip install -r pipeline/requirements.txt -r tests/requirements.txt
python -m pytest tests/test_pipeline.py      # builds a toy genome's data and checks it
node --test tests/js/browser.test.mjs        # the browser's modules on that data
```

The toy genome ([`tests/fixture.py`](tests/fixture.py)) has soft-masked
repeats, a run of Ns, scaffolds, and genes on both strands and on chrM whose
proteins are known in advance, so the tests check reading bases over range
requests, the summary bins, gene search, and splicing and translation end to
end. [`.github/workflows/tests.yml`](.github/workflows/tests.yml) runs them on
every push.

With the full data built and `dev/serve.py` running, three Playwright scripts
drive a real browser: `python tests/smoke.py` and `python tests/scenarios.py`
click through the site and save screenshots to `tests/output/`, and
`python tests/regressions.py` checks behaviour that once went wrong (links that
keep their zoom, late Ensembl answers, pinches, keyboard use and so on).

## Data sources

- Sequence and genes: [Ensembl](https://www.ensembl.org) release 116, GRCh38.p14.
  Ensembl imposes no restrictions on access to, or use of, its data
  ([disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html)).
- Gene names, aliases and cross-references: [HGNC](https://www.genenames.org),
  released under CC0 ([licence](https://www.genenames.org/about/license/)).
  Please cite *Genenames.org: the HGNC resources in 2023*, Nucleic Acids Research.
- Cytobands, transcripts, variants and the regulatory build: the
  [Ensembl REST API](https://rest.ensembl.org).
- Gene summaries: NCBI Gene, through [MyGene.info](https://mygene.info).

## Related

- [Gpxil](https://github.com/comparativechrono/Gpxil): the gpxil specification.
- [Gpxil-tools](https://github.com/comparativechrono/Gpxil-tools): `gpxiltools`,
  the Python reference implementation the pipeline uses.

## Licence

No licence has been chosen for this repository yet.
