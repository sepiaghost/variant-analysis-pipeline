"""Synthetic data with known truth. No real genomes or patients.

Germline: one sample. Somatic: a NORMAL and TUMOR column in one VCF, as many
callers emit. The generator records exactly which sites are true somatic
mutations (clonal and subclonal) so the caller can be scored.
"""
import random

from .vcf import Variant

CONTIGS = {"chr1": 1_000_000, "chr2": 800_000, "chr3": 600_000}
BASES = "ACGT"
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def _binom(rng, n, p):
    return sum(rng.random() < p for _ in range(n))


def _depth(rng, mean):
    return max(0, int(rng.gauss(mean, mean ** 0.5 * 1.3)))


def _alt_for(rng, ref, ti_frac=2 / 3):
    if rng.random() < ti_frac:
        return TRANSITION[ref]
    return rng.choice([b for b in BASES if b not in (ref, TRANSITION[ref])])


def _positions(rng, n):
    seen, out = set(), []
    while len(out) < n:
        c = rng.choice(list(CONTIGS))
        p = rng.randint(1000, CONTIGS[c] - 1000)
        if (c, p) not in seen:
            seen.add((c, p))
            out.append((c, p))
    out.sort(key=lambda x: (x[0], x[1]))
    return out


def make_genes(seed=0, n=40):
    """Fake gene intervals (chrom, start, end, name) for the annotation step."""
    rng = random.Random(seed)
    genes = []
    for i in range(n):
        c = rng.choice(list(CONTIGS))
        s = rng.randint(1000, CONTIGS[c] - 20000)
        genes.append((c, s, s + rng.randint(2000, 15000), f"GENE{i + 1:02d}"))
    return sorted(genes)


def make_germline(seed=1, n_sites=3000, depth=30):
    rng = random.Random(seed)
    variants = []
    for c, p in _positions(rng, n_sites):
        ref = rng.choice(BASES)
        r = rng.random()
        if r < 0.10:                                   # insertion
            alt = ref + "".join(rng.choice(BASES) for _ in range(rng.randint(1, 3)))
        elif r < 0.20:                                 # deletion
            alt, ref = ref, ref + "".join(rng.choice(BASES) for _ in range(rng.randint(1, 3)))
        else:
            alt = _alt_for(rng, ref)
        het = rng.random() < 0.6
        junk = rng.random() < 0.06                     # low-confidence site to exercise filters
        dp = _depth(rng, 5 if junk else depth)
        k = _binom(rng, dp, 0.5 if het else 0.98)
        gt = "0/1" if het else "1/1"
        qual = rng.uniform(3, 18) if junk else rng.uniform(40, 99)
        gq = rng.randint(5, 18) if junk else rng.randint(35, 99)
        variants.append(Variant(c, p, ref, alt, qual, "PASS",
                                {"SAMPLE": {"GT": gt, "DP": dp, "AD": (dp - k, k), "GQ": gq}}))
    return variants


def make_tumor_normal(seed=2, n_germline=1500, n_clonal=150, n_subclonal=40,
                      n_noise=250, purity=0.7, tumor_depth=60, normal_depth=35,
                      normal_contam=0.01, noise_alt_max=3):
    """Return (variants, truth) where truth maps (chrom,pos,ref,alt) -> 'clonal'|'subclonal'.

    * germline het/hom sites appear in both samples (VAF ~0.5 / ~1.0),
    * clonal somatic VAF ~ purity/2, subclonal ~ 0.17 * purity (diploid, copy-neutral),
    * `normal_contam` is the chance a somatic read leaks into the normal (tumor-in-normal),
    * noise sites: 1..noise_alt_max alt reads in the tumor only (sequencing error; raise
      noise_alt_max to model stronger recurrent artifacts).
    """
    rng = random.Random(seed)
    sites = _positions(rng, n_germline + n_clonal + n_subclonal + n_noise)
    rng.shuffle(sites)
    kinds = (["germline"] * n_germline + ["clonal"] * n_clonal +
             ["subclonal"] * n_subclonal + ["noise"] * n_noise)
    variants, truth = [], {}
    for (c, p), kind in zip(sites, kinds):
        ref = rng.choice(BASES)
        alt = _alt_for(rng, ref, ti_frac=0.55 if kind in ("clonal", "subclonal") else 2 / 3)
        td, nd = _depth(rng, tumor_depth), _depth(rng, normal_depth)
        if kind == "germline":
            p_alt = 0.5 if rng.random() < 0.6 else 0.98
            ta, na = _binom(rng, td, p_alt), _binom(rng, nd, p_alt)
        elif kind in ("clonal", "subclonal"):
            vaf = purity / 2 if kind == "clonal" else 0.17 * purity
            ta = _binom(rng, td, vaf)
            na = _binom(rng, nd, normal_contam * vaf)
            truth[(c, p, ref, alt)] = kind
        else:
            ta, na = rng.randint(1, noise_alt_max), 0
        gt_t = "0/1" if ta / max(td, 1) > 0.2 else "0/0"
        gt_n = "0/1" if na / max(nd, 1) > 0.2 else "0/0"
        variants.append(Variant(c, p, ref, alt, rng.uniform(30, 90), "PASS", {
            "NORMAL": {"GT": gt_n, "DP": nd, "AD": (nd - na, na), "GQ": 60},
            "TUMOR": {"GT": gt_t, "DP": td, "AD": (td - ta, ta), "GQ": 60}}))
    return variants, truth
