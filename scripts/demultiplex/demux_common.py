"""Shared helpers for the demultiplex branch (scripts/demultiplex/, route 00b).

Imported by the numbered demultiplex steps; not run directly.
- chromosomes(species): nuclear chromosomes used for genotype calling. Autosomes
  plus X; Y and chrM are left out (Y is absent in female donors, and the
  mitochondrial genome is shared-ish within a pool and heteroplasmic).
- fold_of(barcode): the stable A/B split used to hold out singlets. md5 of the
  barcode, so the same cell always lands in the same fold across reruns.
- read_barcodes / read_sites: the small file formats passed between steps.
"""
import hashlib

_AUTOSOMES = {"human": 22, "mouse": 19}


def chromosomes(species):
    if species not in _AUTOSOMES:
        raise SystemExit(f"species must be one of {sorted(_AUTOSOMES)} (config/pipeline.config), got: {species}")
    return ["chr%s" % c for c in list(range(1, _AUTOSOMES[species] + 1)) + ["X"]]


def fold_of(barcode):
    return "A" if int(hashlib.md5(barcode.encode()).hexdigest(), 16) % 2 == 0 else "B"


def read_barcodes(path):
    return [line.strip() for line in open(path)]


def read_sites(path):
    """(chrom, pos, ref, alt) for each record of a sites-only VCF."""
    sites = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t", 5)
        sites.append((f[0], int(f[1]), f[3], f[4]))
    return sites
