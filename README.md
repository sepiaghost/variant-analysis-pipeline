# Variant Analysis Pipeline

Pure-Python (standard library only) VCF analysis in two stages:

1. **Germline**: hard QC filters (QUAL, depth, GQ), variant typing, Ti/Tv, het/hom ratio,
   per-chromosome counts, gene-interval annotation.
2. **Somatic (tumor-normal)**: classifies every candidate site in a two-sample VCF as somatic or
   rejects it with a named reason (low depth, germline, normal support, low alt reads, low VAF);
   substitution spectrum, tumor VAF histogram, and a purity estimate.

Output is a command-line summary plus a self-contained HTML report with inline SVG charts.

```
python -m variant_pipeline demo                                   # synthetic data, both stages
python -m variant_pipeline germline my.vcf.gz --sample S1
python -m variant_pipeline somatic tn.vcf --tumor TUMOR --normal NORMAL
python -m variant_pipeline.evaluate                               # scenario sweep -> out/evaluation.json
python -m pytest tests                                            # 20 tests
```

This **filters candidates a caller already emitted**; it is not a variant caller. It reads plain and
`.gz` VCFs, splits multi-allelic records, and needs AD/DP in the FORMAT column for the somatic stage.

## Measured results (synthetic data with known truth, mean of 20 seeds)

Each scenario regenerates a tumor-normal VCF with 150 clonal and 40 subclonal true somatic SNVs, 1,500
germline sites and 250 noise sites, then scores the calls against the truth.

| scenario | precision | recall, clonal | recall, subclonal | false positives / run |
|---|---|---|---|---|
| baseline (purity 0.7, depth 60/35) | 1.00 | 0.95 | 0.88 | 0.0 |
| low purity 0.3 | 1.00 | 0.95 | 0.33 | 0.0 |
| low purity 0.15 | 1.00 | 0.64 | 0.05 | 0.0 |
| high purity 0.9 | 1.00 | 0.93 | 0.95 | 0.0 |
| low depth 25/15 | 1.00 | 0.82 | 0.32 | 0.3 |
| high depth 120/60 | 1.00 | 0.98 | 0.99 | 0.0 |
| tumor-in-normal leak 0.10 | 1.00 | 0.49 | 0.73 | 0.0 |
| tumor-in-normal leak 0.30 | 1.00 | 0.07 | 0.41 | 0.0 |
| stronger artifacts (up to 6 alt reads) | **0.59** | 0.95 | 0.88 | **123.0** |

Germline stage on the synthetic file: Ti/Tv 2.05 and het/hom 1.48, matching how it was generated
(2.0 and 1.5), with the 6% injected low-confidence sites removed by the QC filters.

## Limits (read these)

- **All numbers are on synthetic data** built by the same author as the filter. They show the logic
  behaves as designed and where it breaks; they are not a benchmark on real tumors. It has not been
  validated against a real caller's output.
- **No strand-bias, base-quality, mapping-quality or sequence-context evidence.** Artifacts that
  reach 4+ alt reads pass, which is the precision collapse in the last row. A real pipeline needs
  those filters (or a panel of normals).
- **Contamination of the normal by tumor DNA** hurts recall sharply (clonal recall 0.49 at a 0.10
  leak), because the "normal support" rule rejects real mutations.
- **Subclonal calls need purity and depth.** At purity 0.3 or low depth, most subclonal mutations
  are missed.
- The purity estimate (2 x VAF mode) assumes a diploid, copy-neutral dominant clone.
- SNV-focused: indels are typed and counted in the germline stage but excluded from the spectrum.
- Positions and genes in the demo are invented; nothing here uses real genomes or patient data.
