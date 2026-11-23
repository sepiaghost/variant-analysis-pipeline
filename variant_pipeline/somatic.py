"""Tumor-normal somatic calling from a two-sample VCF.

This filters candidate sites that a caller already emitted; it is not a
variant caller. A candidate is SOMATIC when the tumor shows enough alt support and
the matched normal shows essentially none.
"""
from collections import Counter

from .germline import TRANSITIONS

DEFAULTS = dict(min_tumor_dp=10, min_normal_dp=8, min_tumor_alt=4,
                min_tumor_vaf=0.05, max_normal_vaf=0.03, max_normal_alt=1)

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
SPECTRUM_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def vaf(ad):
    n = ad[0] + ad[1]
    return ad[1] / n if n else 0.0


def classify(v, tumor, normal, **kw):
    p = {**DEFAULTS, **kw}
    t, n = v.samples[tumor], v.samples[normal]
    if t["AD"] is None or n["AD"] is None:
        return "no_allele_depth"
    td, nd = sum(t["AD"]), sum(n["AD"])
    tv, nv = vaf(t["AD"]), vaf(n["AD"])
    if td < p["min_tumor_dp"]:
        return "low_tumor_depth"
    if nd < p["min_normal_dp"]:
        return "low_normal_depth"
    if nv >= 0.3:
        return "germline"            # clearly present in the normal
    if n["AD"][1] > p["max_normal_alt"] or nv > p["max_normal_vaf"]:
        return "normal_support"      # too much normal evidence: contamination / mosaic
    if t["AD"][1] < p["min_tumor_alt"]:
        return "low_tumor_alt"
    if tv < p["min_tumor_vaf"]:
        return "low_vaf"
    return "SOMATIC"


def call_somatic(variants, tumor="TUMOR", normal="NORMAL", **kw):
    calls, reasons = [], Counter()
    for v in variants:
        c = classify(v, tumor, normal, **kw)
        reasons[c] += 1
        if c == "SOMATIC":
            calls.append(v)
    return calls, reasons


def spectrum(calls):
    """Six-class SNV substitution spectrum, expressed on the pyrimidine strand."""
    cnt = Counter()
    for v in calls:
        if v.kind != "SNV":
            continue
        r, a = v.ref, v.alt
        if r in "AG":
            r, a = COMPLEMENT[r], COMPLEMENT[a]
        cnt[f"{r}>{a}"] += 1
    return {k: cnt.get(k, 0) for k in SPECTRUM_ORDER}


def vaf_histogram(calls, tumor="TUMOR", bins=20):
    h = [0] * bins
    for v in calls:
        i = min(int(vaf(v.samples[tumor]["AD"]) * bins), bins - 1)
        h[i] += 1
    return h


def clonal_summary(calls, tumor="TUMOR", bins=20):
    """Mode of the VAF histogram, and the purity it implies *if* the tumor is diploid and
    copy-neutral at those sites (an assumption, stated in the report)."""
    h = vaf_histogram(calls, tumor, bins)
    if not any(h):
        return {"mode_vaf": None, "implied_purity": None}
    mode = (h.index(max(h)) + 0.5) / bins
    return {"mode_vaf": mode, "implied_purity": min(1.0, 2 * mode)}


def score_against_truth(calls, truth):
    called = {(v.chrom, v.pos, v.ref, v.alt) for v in calls}
    tp = called & truth.keys()
    out = {"called": len(called), "true_positive": len(tp),
           "false_positive": len(called - truth.keys()),
           "precision": len(tp) / len(called) if called else None}
    for cls in ("clonal", "subclonal"):
        want = {k for k, v in truth.items() if v == cls}
        out[f"recall_{cls}"] = len(want & called) / len(want) if want else None
    out["recall_all"] = len(tp) / len(truth) if truth else None
    return out
