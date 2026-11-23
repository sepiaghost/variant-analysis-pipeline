import gzip
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from variant_pipeline import germline as G
from variant_pipeline import somatic as S
from variant_pipeline.cli import main
from variant_pipeline.synth import make_germline, make_tumor_normal, make_genes, CONTIGS
from variant_pipeline.vcf import Variant, read_vcf, write_vcf

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def test_roundtrip(tmp_path):
    vs = make_germline(seed=3, n_sites=50)
    p = tmp_path / "g.vcf"
    write_vcf(p, vs, ["SAMPLE"], CONTIGS)
    names, back = read_vcf(p)
    assert names == ["SAMPLE"] and len(back) == 50
    assert (back[7].chrom, back[7].pos, back[7].ref, back[7].alt) == (vs[7].chrom, vs[7].pos, vs[7].ref, vs[7].alt)
    assert back[7].samples["SAMPLE"]["AD"] == vs[7].samples["SAMPLE"]["AD"]


def test_multiallelic_is_split(tmp_path):
    p = tmp_path / "m.vcf"
    p.write_text(HEADER + "chr1\t100\t.\tA\tG,T\t50\tPASS\t.\tGT:DP:AD:GQ\t1/2:30:5,12,13:99\n")
    _, vs = read_vcf(p)
    assert [v.alt for v in vs] == ["G", "T"]
    assert vs[0].samples["S1"]["AD"] == (5, 12) and vs[1].samples["S1"]["AD"] == (5, 13)


def test_gzip_and_malformed(tmp_path):
    p = tmp_path / "z.vcf.gz"
    with gzip.open(p, "wt") as f:
        f.write(HEADER + "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:AD:GQ\t0/1:30:15,15:99\n")
    assert len(read_vcf(p)[1]) == 1
    bad = tmp_path / "bad.vcf"
    bad.write_text(HEADER + "chr1\t100\tA\n")
    with pytest.raises(ValueError):
        read_vcf(bad)


def test_variant_kinds():
    k = lambda r, a: Variant("c", 1, r, a, 1, "PASS").kind
    assert (k("A", "G"), k("A", "AT"), k("AT", "A"), k("AC", "GT")) == ("SNV", "INS", "DEL", "MNV")


def test_zygosity():
    assert [G.zygosity(x) for x in ("0/1", "1|0", "1/1", "0/0", "./.")] == \
        ["het", "het", "hom_alt", "hom_ref", "missing"]


def test_germline_qc_and_summary():
    vs = make_germline(seed=5)
    kept, reasons = G.qc_filter(vs, "SAMPLE")
    assert len(kept) + sum(reasons.values()) == len(vs)
    assert reasons["low_qual"] > 0
    s = G.summarize(kept, "SAMPLE")
    assert 1.8 < s["ti_tv"] < 2.3          # generated at 2:1
    assert 1.2 < s["het_hom_ratio"] < 1.8  # generated at 60/40


def test_annotation():
    genes = [("chr1", 100, 200, "G1"), ("chr1", 500, 900, "G2")]
    vs = [Variant("chr1", p, "A", "G", 50, "PASS") for p in (150, 600, 300)]
    lab = G.annotate(vs, genes)
    assert [lab[("chr1", p, "A", "G")] for p in (150, 600, 300)] == ["G1", "G2", "intergenic"]


def _site(t_ad, n_ad):
    return Variant("c", 1, "C", "T", 50, "PASS", {
        "TUMOR": {"GT": "0/1", "DP": sum(t_ad), "AD": t_ad, "GQ": 60},
        "NORMAL": {"GT": "0/0", "DP": sum(n_ad), "AD": n_ad, "GQ": 60}})


@pytest.mark.parametrize("t,n,expected", [
    ((40, 20), (35, 0), "SOMATIC"),
    ((5, 3), (35, 0), "low_tumor_depth"),
    ((40, 20), (4, 0), "low_normal_depth"),
    ((30, 30), (18, 17), "germline"),
    ((40, 20), (33, 2), "normal_support"),
    ((55, 3), (35, 0), "low_tumor_alt"),
    ((95, 5), (35, 0), "SOMATIC"),          # VAF 0.05 sits exactly on the cutoff and passes
    ((90, 3), (35, 0), "low_tumor_alt"),
])
def test_classify(t, n, expected):
    assert S.classify(_site(t, n), "TUMOR", "NORMAL") == expected


def test_spectrum_folds_to_pyrimidine_strand():
    vs = [Variant("c", 1, "G", "A", 1, "PASS"), Variant("c", 2, "C", "T", 1, "PASS"),
          Variant("c", 3, "A", "C", 1, "PASS"), Variant("c", 4, "AT", "A", 1, "PASS")]
    sp = S.spectrum(vs)
    assert sp["C>T"] == 2 and sp["T>G"] == 1 and sum(sp.values()) == 3  # indel ignored


def test_somatic_on_synthetic_truth():
    vs, truth = make_tumor_normal(seed=9)
    calls, reasons = S.call_somatic(vs)
    sc = S.score_against_truth(calls, truth)
    assert sc["precision"] > 0.97 and sc["recall_clonal"] > 0.9
    assert reasons["germline"] > 1000           # germline sites are not called somatic
    assert 0.6 < S.clonal_summary(calls)["implied_purity"] < 0.8


def test_known_failure_stronger_artifacts_hurt_precision():
    """Documented limitation: with no strand/context evidence, artifacts that reach 4-6 alt reads pass."""
    vs, truth = make_tumor_normal(seed=9, noise_alt_max=6)
    sc = S.score_against_truth(S.call_somatic(vs)[0], truth)
    assert sc["precision"] < 0.8


def test_cli_demo_end_to_end(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    main(["demo"])
    out = capsys.readouterr().out
    assert "precision" in out
    for f in ("germline_report.html", "somatic_report.html"):
        assert "<svg" in (tmp_path / "out" / f).read_text()


def test_cli_rejects_unknown_sample(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    main(["demo"])
    with pytest.raises(SystemExit):
        main(["somatic", "out/tumor_normal.vcf", "--tumor", "NOPE", "--normal", "NORMAL"])
