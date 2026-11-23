"""Command line: python -m variant_pipeline <command> ...

  demo                              build synthetic data, run both stages, write reports to out/
  germline  VCF [--sample S]        QC + summary + HTML report
  somatic   VCF --tumor T --normal N   tumor-normal filtering + HTML report
"""
import argparse
import os
import sys

from . import germline as G
from . import somatic as S
from . import report as R
from .synth import make_germline, make_tumor_normal, make_genes, CONTIGS
from .vcf import read_vcf, write_vcf


def run_germline(vcf, sample, out, genes=None, **qc):
    names, variants = read_vcf(vcf)
    sample = sample or names[0]
    if sample not in names:
        sys.exit(f"sample {sample!r} not in VCF (found: {', '.join(names)})")
    kept, reasons = G.qc_filter(variants, sample, **qc)
    summary = G.summarize(kept, sample)
    hit = 0
    if genes:
        lab = G.annotate(kept, genes)
        hit = sum(v != "intergenic" for v in lab.values())
    with open(out, "w") as f:
        f.write(R.germline_report(summary, reasons, hit, len(variants)))
    return summary, reasons


def run_somatic(vcf, tumor, normal, out, truth=None, **kw):
    names, variants = read_vcf(vcf)
    for n in (tumor, normal):
        if n not in names:
            sys.exit(f"sample {n!r} not in VCF (found: {', '.join(names)})")
    calls, reasons = S.call_somatic(variants, tumor, normal, **kw)
    score = S.score_against_truth(calls, truth) if truth else None
    with open(out, "w") as f:
        f.write(R.somatic_report(reasons, S.spectrum(calls), S.vaf_histogram(calls, tumor),
                                 S.clonal_summary(calls, tumor), score, len(variants),
                                 {**S.DEFAULTS, **kw}))
    return calls, reasons, score


def main(argv=None):
    ap = argparse.ArgumentParser(prog="variant_pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    g = sub.add_parser("germline")
    g.add_argument("vcf"); g.add_argument("--sample"); g.add_argument("--out", default="germline_report.html")
    s = sub.add_parser("somatic")
    s.add_argument("vcf"); s.add_argument("--tumor", required=True); s.add_argument("--normal", required=True)
    s.add_argument("--out", default="somatic_report.html")
    a = ap.parse_args(argv)

    if a.cmd == "demo":
        os.makedirs("out", exist_ok=True)
        write_vcf("out/germline.vcf", make_germline(), ["SAMPLE"], CONTIGS)
        tn, truth = make_tumor_normal()
        write_vcf("out/tumor_normal.vcf", tn, ["NORMAL", "TUMOR"], CONTIGS)
        summ, reasons = run_germline("out/germline.vcf", None, "out/germline_report.html", make_genes())
        calls, sreasons, score = run_somatic("out/tumor_normal.vcf", "TUMOR", "NORMAL",
                                             "out/somatic_report.html", truth)
        print(f"germline: kept {summ['total']}, Ti/Tv {summ['ti_tv']:.2f}, "
              f"het/hom {summ['het_hom_ratio']:.2f}, removed {dict(reasons)}")
        print(f"somatic : {score['called']} called, precision {score['precision']:.2f}, "
              f"recall clonal {score['recall_clonal']:.2f}, subclonal {score['recall_subclonal']:.2f}")
        print("wrote out/germline_report.html, out/somatic_report.html")
    elif a.cmd == "germline":
        summ, reasons = run_germline(a.vcf, a.sample, a.out)
        print(f"kept {summ['total']} variants; removed {dict(reasons)}; wrote {a.out}")
    else:
        calls, reasons, _ = run_somatic(a.vcf, a.tumor, a.normal, a.out)
        print(f"{len(calls)} somatic calls; outcomes {dict(reasons)}; wrote {a.out}")


if __name__ == "__main__":
    main()
