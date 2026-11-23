"""Germline QC and summary for a single-sample VCF."""
import bisect
from collections import Counter

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def qc_filter(variants, sample, min_qual=30.0, min_dp=10, min_gq=20):
    """Apply hard filters. Returns (kept, reasons Counter). A record can fail several
    filters; it is counted once, under its first failing reason."""
    kept, reasons = [], Counter()
    for v in variants:
        s = v.samples[sample]
        if v.filter not in ("PASS", "."):
            reasons["caller_filter"] += 1
        elif v.qual is not None and v.qual < min_qual:
            reasons["low_qual"] += 1
        elif s["DP"] is None or s["DP"] < min_dp:
            reasons["low_depth"] += 1
        elif s["GQ"] is not None and s["GQ"] < min_gq:
            reasons["low_gq"] += 1
        else:
            kept.append(v)
    return kept, reasons


def zygosity(gt):
    alleles = gt.replace("|", "/").split("/")
    if "." in alleles:
        return "missing"
    n_alt = sum(a != "0" for a in alleles)
    return "hom_alt" if n_alt == len(alleles) else ("het" if n_alt else "hom_ref")


def summarize(variants, sample):
    kinds = Counter(v.kind for v in variants)
    snvs = [v for v in variants if v.kind == "SNV"]
    ti = sum((v.ref, v.alt) in TRANSITIONS for v in snvs)
    tv = len(snvs) - ti
    zyg = Counter(zygosity(v.samples[sample]["GT"]) for v in variants)
    per_chrom = Counter(v.chrom for v in variants)
    return {
        "total": len(variants),
        "by_type": dict(kinds),
        "ti_tv": ti / tv if tv else None,
        "het": zyg["het"], "hom_alt": zyg["hom_alt"],
        "het_hom_ratio": zyg["het"] / zyg["hom_alt"] if zyg["hom_alt"] else None,
        "per_chrom": dict(sorted(per_chrom.items())),
    }


def annotate(variants, genes):
    """Label each variant with the gene interval it falls in, else 'intergenic'.
    genes: list of (chrom, start, end, name), sorted."""
    by_chrom = {}
    for c, s, e, n in genes:
        by_chrom.setdefault(c, []).append((s, e, n))
    starts = {c: [g[0] for g in gs] for c, gs in by_chrom.items()}
    out = {}
    for v in variants:
        name = "intergenic"
        gs = by_chrom.get(v.chrom, [])
        i = bisect.bisect_right(starts.get(v.chrom, []), v.pos) - 1
        # genes may overlap, so scan back a few entries
        for j in range(i, max(i - 5, -1), -1):
            if gs[j][0] <= v.pos <= gs[j][1]:
                name = gs[j][2]
                break
        out[(v.chrom, v.pos, v.ref, v.alt)] = name
    return out
