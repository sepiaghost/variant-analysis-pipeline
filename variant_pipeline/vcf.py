"""Minimal VCF reader/writer using only the standard library.

Multi-allelic records are split into one record per ALT allele, with AD
narrowed to (ref, that alt). Handles plain and .gz files.
"""
import gzip
from dataclasses import dataclass, field


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float | None
    filter: str
    samples: dict = field(default_factory=dict)  # name -> {"GT","DP","AD":(ref,alt),"GQ"}

    @property
    def kind(self):
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "SNV"
        if len(self.ref) == len(self.alt):
            return "MNV"
        return "INS" if len(self.alt) > len(self.ref) else "DEL"


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _num(x):
    return None if x in (".", "") else float(x)


def read_vcf(path):
    """Return (sample_names, [Variant])."""
    names, out = [], []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = f[9:]
                continue
            if len(f) < 8:
                raise ValueError(f"malformed VCF line (needs 8+ columns): {line[:60]!r}")
            alts = f[4].split(",")
            fmt = f[8].split(":") if len(f) > 8 else []
            for ai, alt in enumerate(alts, start=1):
                v = Variant(f[0], int(f[1]), f[3], alt, _num(f[5]), f[6])
                for name, raw in zip(names, f[9:]):
                    d = dict(zip(fmt, raw.split(":")))
                    rec = {"GT": d.get("GT", "./.")}
                    rec["DP"] = int(d["DP"]) if d.get("DP", ".") != "." else None
                    rec["GQ"] = int(d["GQ"]) if d.get("GQ", ".") != "." else None
                    ad = d.get("AD", ".")
                    if ad != ".":
                        parts = [int(p) for p in ad.split(",")]
                        rec["AD"] = (parts[0], parts[ai]) if ai < len(parts) else (parts[0], 0)
                    else:
                        rec["AD"] = None
                    v.samples[name] = rec
                out.append(v)
    return names, out


def write_vcf(path, variants, sample_names, contigs):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, n in contigs.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        fh.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n')
        fh.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n")
        for v in variants:
            cols = [v.chrom, str(v.pos), ".", v.ref, v.alt,
                    "." if v.qual is None else f"{v.qual:.1f}", v.filter, ".", "GT:DP:AD:GQ"]
            for n in sample_names:
                s = v.samples[n]
                cols.append(f"{s['GT']}:{s['DP']}:{s['AD'][0]},{s['AD'][1]}:{s['GQ']}")
            fh.write("\t".join(cols) + "\n")
