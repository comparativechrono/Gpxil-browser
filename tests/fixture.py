"""A small synthetic genome for the tests: FASTA, GTF, HGNC table and cytobands.

    python tests/fixture.py --out tests/.fixture/inputs

It looks like a tiny GRCh38 to the pipeline: chromosomes 1, 2 and MT and two
unplaced scaffolds, with soft-masked repeats (a 171-base tandem repeat, like
alpha satellite), a run of Ns, and three genes whose proteins are known in
advance, so the tests can check splicing and translation end to end:

- PIXL1 (+ strand, chr1): three exons, a MANE Select canonical transcript and
  a second, non-canonical one; one TGA codon is read as selenocysteine (U).
- PIXL2 (- strand, chr1): three exons.
- MT-PIX (+ strand, chrM): uses the vertebrate mitochondrial code (TGA = W,
  ATA = M, AGA = stop).

Also writes ``expected.json`` with the proteins and the gene coordinates.
"""

from __future__ import annotations

import argparse
import json
import os
import random

BASES = "TCAG"
STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
VERT_MITO = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
AMINO = "ACDEFGHIKLMNPQRSTVWY"


def codons_for(table: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for i, aa in enumerate(table):
        codon = BASES[i // 16] + BASES[(i // 4) % 4] + BASES[i % 4]
        out.setdefault(aa, []).append(codon)
    return out


def translate(cds: str, table: str) -> str:
    return "".join(table[BASES.index(cds[i]) * 16 + BASES.index(cds[i + 1]) * 4 + BASES.index(cds[i + 2])] for i in range(0, len(cds) - 2, 3))


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def random_dna(rng: random.Random, n: int, gc: float = 0.42) -> list[str]:
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return rng.choices("ACGT", weights=w, k=n)


class Gene:
    def __init__(self, gid, name, contig, strand, exons, cds, biotype="protein_coding", mito=False, tags=("Ensembl_canonical",), sec=()):
        self.gid, self.name, self.contig, self.strand = gid, name, contig, strand
        self.exons = sorted(exons)  # genomic, 0-based half-open
        self.cds = cds  # genomic [start, end) including the stop codon, or None
        self.biotype, self.mito, self.tags = biotype, mito, tags
        self.sec = tuple(sec)  # codon numbers (0-based) written as TGA and read as selenocysteine
        self.sec_genomic = []  # (start0, end) of each such codon, for the GTF
        self.protein = None

    def tx_order(self):
        return self.exons if self.strand > 0 else self.exons[::-1]

    def mrna_positions(self):
        pos = []
        for a, b in self.tx_order():
            pos += list(range(a, b)) if self.strand > 0 else list(range(b - 1, a - 1, -1))
        return pos

    def place(self, rng: random.Random, chrom: list[str]):
        """Write an mRNA with a known protein into the chromosome."""
        pos = self.mrna_positions()
        mrna = random_dna(rng, len(pos), 0.5)
        if self.cds:
            idx = [i for i, p in enumerate(pos) if self.cds[0] <= p < self.cds[1]]
            first, last = idx[0], idx[-1] + 1
            n = last - first
            assert n % 3 == 0, f"{self.name}: CDS of {n} bases"
            table = VERT_MITO if self.mito else STANDARD
            syn = codons_for(table)
            aas = ["M"] + rng.choices(AMINO, k=n // 3 - 2)
            if self.mito:
                aas[3] = aas[9] = "W"  # will be written as TGA, a stop in the standard code
                aas[5] = "M"
            for i in self.sec:
                aas[i] = "U"
            codons = []
            for i, aa in enumerate(aas):
                if aa == "U":
                    codons.append("TGA")
                    gpos = pos[first + 3 * i : first + 3 * i + 3]
                    assert max(gpos) - min(gpos) == 2, "keep selenocysteine codons inside one exon"
                    self.sec_genomic.append((min(gpos), max(gpos) + 1))
                elif self.mito and aa == "W" and i in (3, 9):
                    codons.append("TGA")
                elif self.mito and i == 5:
                    codons.append("ATA")  # Ile in the standard code, Met here
                else:
                    codons.append(rng.choice(syn[aa]))
            codons.append("AGA" if self.mito else rng.choice(syn["*"]))
            cds = "".join(codons)
            assert len(cds) == n
            mrna[first:last] = list(cds)
            protein = list(translate(cds, table))
            for i in self.sec:
                assert protein[i] == "*"
                protein[i] = "U"
            self.protein = "".join(protein)
            assert self.protein == "".join(aas) + "*"
        seq = "".join(mrna)
        for i, p in enumerate(pos):
            chrom[p] = seq[i] if self.strand > 0 else seq[i].translate(COMP)

    def gtf(self, contig_name: str, txid: str, txname: str, extra_tx=None) -> list[str]:
        st = "+" if self.strand > 0 else "-"
        a = f'gene_id "{self.gid}"; gene_version "1"; gene_name "{self.name}"; gene_source "ensembl_havana"; gene_biotype "{self.biotype}";'
        lo, hi = self.exons[0][0], self.exons[-1][1]
        out = [f"{contig_name}\tensembl_havana\tgene\t{lo + 1}\t{hi}\t.\t{st}\t.\t{a}"]
        out += self._tx_lines(contig_name, st, a, txid, txname, self.exons, self.cds, self.tags)
        if extra_tx:
            out += self._tx_lines(contig_name, st, a, *extra_tx)
        return out

    def _tx_lines(self, contig_name, st, a, txid, txname, exons, cds, tags):
        tagtext = " ".join(f'tag "{t}";' for t in tags)
        t = f'{a} transcript_id "{txid}"; transcript_version "1"; transcript_name "{txname}"; transcript_source "ensembl_havana"; transcript_biotype "{self.biotype}"; {tagtext}'
        lo, hi = exons[0][0], exons[-1][1]
        lines = [f"{contig_name}\tensembl_havana\ttranscript\t{lo + 1}\t{hi}\t.\t{st}\t.\t{t}"]
        order = exons if self.strand > 0 else exons[::-1]
        done = 0  # CDS bases so far, for the frame column
        for k, (s, e) in enumerate(order, 1):
            lines.append(f'{contig_name}\tensembl_havana\texon\t{s + 1}\t{e}\t.\t{st}\t.\t{t} exon_number "{k}"; exon_id "{txid}E{k}";')
            if cds:
                # CDS lines exclude the stop codon, which gets its own line (as in Ensembl's GTF)
                stop = (cds[1] - 3, cds[1]) if self.strand > 0 else (cds[0], cds[0] + 3)
                cs, ce = (cds[0], cds[1] - 3) if self.strand > 0 else (cds[0] + 3, cds[1])
                ps, pe = max(s, cs), min(e, ce)
                if ps < pe:
                    frame = (3 - done % 3) % 3
                    protein = txid.replace("ENST", "ENSP")
                    lines.append(f'{contig_name}\tensembl_havana\tCDS\t{ps + 1}\t{pe}\t.\t{st}\t{frame}\t{t} exon_number "{k}"; protein_id "{protein}"; protein_version "1";')
                    done += pe - ps
                qs, qe = max(s, stop[0]), min(e, stop[1])
                if qs < qe:
                    lines.append(f'{contig_name}\tensembl_havana\tstop_codon\t{qs + 1}\t{qe}\t.\t{st}\t0\t{t} exon_number "{k}";')
                if "Ensembl_canonical" in tags:
                    for a, b in self.sec_genomic:
                        if s <= a and b <= e:
                            lines.append(f'{contig_name}\tensembl_havana\tSelenocysteine\t{a + 1}\t{b}\t.\t{st}\t.\t{t} exon_number "{k}";')
        return lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=2009)
    a = ap.parse_args(argv)
    rng = random.Random(a.seed)
    os.makedirs(a.out, exist_ok=True)

    lengths = {"1": 400_000, "2": 120_000, "MT": 16_569, "KI270728.1": 5_000, "GL000009.2": 3_000}
    seqs = {name: random_dna(rng, n) for name, n in lengths.items()}

    # a soft-masked tandem repeat with a 171-base unit (a few mutations), and a gap of Ns
    unit = random_dna(rng, 171, 0.38)
    for i in range(50_000, 60_000):
        b = unit[(i - 50_000) % 171]
        seqs["1"][i] = (rng.choice("ACGT") if rng.random() < 0.02 else b).lower()
    for i in range(200_000, 210_000):
        seqs["1"][i] = "N"
    # a soft-masked stretch on the scaffold too
    for i in range(1_000, 1_500):
        seqs["KI270728.1"][i] = seqs["KI270728.1"][i].lower()

    genes = [
        Gene("ENSG00000000001", "PIXL1", "1", 1, [(100_000, 100_100), (105_000, 105_200), (110_000, 110_300)], (100_050, 110_101), tags=("basic", "Ensembl_canonical", "MANE_Select", "seleno"), sec=(40,)),
        Gene("ENSG00000000002", "PIXL2", "1", -1, [(300_000, 300_150), (302_000, 302_090), (304_000, 304_240)], (300_100, 304_199)),
        Gene("ENSG00000000003", "MT-PIX", "MT", 1, [(3_000, 3_600)], (3_010, 3_571), mito=True),
        Gene("ENSG00000000004", "PIXL1-AS1", "1", -1, [(99_000, 99_400), (99_800, 100_020)], None, biotype="lncRNA"),
    ]
    for g in sorted(genes, key=lambda g: g.cds is not None):  # non-coding first: coding genes win overlaps
        g.place(rng, seqs[g.contig])

    with open(os.path.join(a.out, "toy.fa"), "w") as fh:
        for name, s in seqs.items():
            kind = "chromosome" if name in ("1", "2", "MT") else "scaffold"
            fh.write(f">{name} dna_sm:{kind} {kind}:GRCh38:{name}:1:{len(s)}:1 REF\n")
            text = "".join(s)
            for i in range(0, len(text), 60):
                fh.write(text[i : i + 60] + "\n")

    with open(os.path.join(a.out, "toy.gtf"), "w") as fh:
        fh.write("#!genome-build GRCh38.p14\n")
        lines = []
        lines += genes[0].gtf("1", "ENST00000000011", "PIXL1-201", extra_tx=("ENST00000000012", "PIXL1-202", [(100_000, 100_100), (110_000, 110_300)], None, ("basic",)))
        lines += genes[1].gtf("1", "ENST00000000021", "PIXL2-201")
        lines += genes[2].gtf("MT", "ENST00000000031", "MT-PIX-201")
        lines += genes[3].gtf("1", "ENST00000000041", "PIXL1-AS1-201")
        # a gene on a contig the FASTA doesn't have, which the pipeline must skip
        lines.append('KI270999.1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000009"; gene_name "GHOST"; gene_biotype "protein_coding";')
        fh.write("\n".join(lines) + "\n")

    cols = ["hgnc_id", "symbol", "name", "locus_group", "status", "location", "alias_symbol", "prev_symbol", "entrez_id", "ensembl_gene_id", "uniprot_ids", "omim_id"]
    rows = [
        ["HGNC:90001", "PIXL1", "pixel protein 1", "protein-coding gene", "Approved", "1p11", "PIX1|PXL", "OLDPIX", "999001", "ENSG00000000001", "Q00001|Q00009", "600001"],
        ["HGNC:90002", "PIXL2", "pixel protein 2", "protein-coding gene", "Approved", "1q21.2", "PXL", "", "999002", "ENSG00000000002", "Q00002", ""],
        ["HGNC:90003", "MT-PIX", "mitochondrially encoded pixel", "protein-coding gene", "Approved", "mitochondria", "", "", "999003", "ENSG00000000003", "", ""],
        ["HGNC:90004", "PIXL1-AS1", "PIXL1 antisense RNA 1", "non-coding RNA", "Approved", "1p11", "", "", "", "ENSG00000000004", "", ""],
        ["HGNC:90005", "WITHDRAWN1", "an entry that was withdrawn", "other", "Entry Withdrawn", "", "PIX1", "", "", "", "", ""],
    ]
    with open(os.path.join(a.out, "hgnc.txt"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    bands = {
        "assembly_name": "GRCh38.toy",
        "karyotype": ["1", "2", "MT"],
        "top_level_region": [
            {
                "name": "1",
                "coord_system": "chromosome",
                "length": lengths["1"],
                "bands": [
                    {"id": "p12", "start": 1, "end": 90_000, "stain": "gneg"},
                    {"id": "p11", "start": 90_001, "end": 190_000, "stain": "gpos50"},
                    {"id": "q11", "start": 190_001, "end": 215_000, "stain": "acen"},
                    {"id": "q21.1", "start": 215_001, "end": 290_000, "stain": "gpos100"},
                    {"id": "q21.2", "start": 290_001, "end": 400_000, "stain": "gneg"},
                ],
            },
            {"name": "2", "coord_system": "chromosome", "length": lengths["2"], "bands": [{"id": "p11", "start": 1, "end": 60_000, "stain": "gneg"}, {"id": "q11", "start": 60_001, "end": 120_000, "stain": "gpos25"}]},
            {"name": "MT", "coord_system": "chromosome", "length": lengths["MT"]},
            {"name": "KI270728.1", "coord_system": "scaffold", "length": lengths["KI270728.1"]},
        ],
    }
    with open(os.path.join(a.out, "bands.json"), "w") as fh:
        json.dump(bands, fh)

    expected = {
        "genes": {
            g.gid: {
                "name": g.name,
                "contig": "chrM" if g.contig == "MT" else f"chr{g.contig}",
                "start": g.exons[0][0],
                "end": g.exons[-1][1],
                "strand": g.strand,
                "exons": g.exons,
                "cds": g.cds,
                "sec": [a if g.strand > 0 else b - 1 for a, b in g.sec_genomic],
                "protein": g.protein,
            }
            for g in genes
        },
        "lengths": {("chrM" if k == "MT" else f"chr{k}" if k in ("1", "2") else k): v for k, v in lengths.items()},
    }
    with open(os.path.join(a.out, "expected.json"), "w") as fh:
        json.dump(expected, fh, indent=1)
    print(f"wrote the toy genome to {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
