"""Score the somatic filter on synthetic scenarios, averaged over seeds.

Run: python -m variant_pipeline.evaluate     (writes out/evaluation.json)
"""
import json
import os
import statistics as st

from .somatic import call_somatic, score_against_truth
from .synth import make_tumor_normal

SCENARIOS = {
    "baseline (purity 0.7, depth 60/35)": {},
    "low purity 0.3": {"purity": 0.3},
    "low purity 0.15": {"purity": 0.15},
    "high purity 0.9": {"purity": 0.9},
    "low depth 25/15": {"tumor_depth": 25, "normal_depth": 15},
    "high depth 120/60": {"tumor_depth": 120, "normal_depth": 60},
    "tumor-in-normal leak 0.10": {"normal_contam": 0.10},
    "tumor-in-normal leak 0.30": {"normal_contam": 0.30},
    "stronger artifacts (up to 6 alt reads)": {"noise_alt_max": 6},
}


def run(n_seeds=20):
    out = {}
    for name, kw in SCENARIOS.items():
        rows = []
        for seed in range(n_seeds):
            variants, truth = make_tumor_normal(seed=100 + seed, **kw)
            calls, _ = call_somatic(variants)
            rows.append(score_against_truth(calls, truth))
        def mean(k):
            vals = [r[k] for r in rows if r[k] is not None]
            return st.mean(vals) if vals else None
        out[name] = {"precision": mean("precision"), "recall_clonal": mean("recall_clonal"),
                     "recall_subclonal": mean("recall_subclonal"),
                     "false_positives_per_run": mean("false_positive")}
    return {"n_seeds": n_seeds, "scenarios": out}


if __name__ == "__main__":
    res = run()
    os.makedirs("out", exist_ok=True)
    json.dump(res, open("out/evaluation.json", "w"), indent=1)
    print(f"{'scenario':42s} prec  rec_clonal rec_subcl FP/run")
    for k, v in res["scenarios"].items():
        f = lambda x: "  n/a" if x is None else f"{x:5.2f}"
        print(f"{k:42s} {f(v['precision'])} {f(v['recall_clonal'])}     {f(v['recall_subclonal'])}    {v['false_positives_per_run']:.1f}")
