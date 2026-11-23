"""Self-contained HTML report (inline SVG charts, no external files or libraries)."""
import html

CSS = """
body{font-family:'IBM Plex Sans',sans-serif;background:#F6F2EE;color:#2A2530;max-width:860px;margin:36px auto;padding:0 24px;line-height:1.6}
h1,h2{font-family:'IBM Plex Mono',monospace} h1{font-size:25px} h2{font-size:18px;margin-top:32px}
table{border-collapse:collapse;width:100%;background:#fff;margin:10px 0}
th,td{border:1px solid #E8E0D9;padding:7px 10px;text-align:left;font-size:14.5px}
th{background:#FBEAEF;font-family:'IBM Plex Mono',monospace}
.note{background:#fff;border-left:3px solid #C2607F;padding:10px 14px;color:#675F6B}
svg text{font-family:'IBM Plex Mono',monospace;font-size:11px;fill:#675F6B}
"""


def _e(x):
    return html.escape(str(x))


def bar_chart(labels, values, width=640, height=200, color="#C2607F"):
    n, top = len(values), max(values) if values and max(values) else 1
    bw = (width - 50) / max(n, 1)
    parts = [f'<svg viewBox="0 0 {width} {height + 56}" width="100%" role="img">']
    for i, (lab, val) in enumerate(zip(labels, values)):
        h = (val / top) * height
        x = 40 + i * bw
        parts.append(f'<rect x="{x + 2:.1f}" y="{height - h + 22:.1f}" width="{bw - 4:.1f}" '
                     f'height="{h:.1f}" fill="{color}" rx="2"><title>{_e(lab)}: {val}</title></rect>')
        parts.append(f'<text x="{x + bw / 2:.1f}" y="{height + 40}" text-anchor="middle">{_e(lab)}</text>')
        if val:
            parts.append(f'<text x="{x + bw / 2:.1f}" y="{height - h + 18:.1f}" text-anchor="middle">{val}</text>')
    parts.append("</svg>")
    return "".join(parts)


def _table(rows, head=None):
    h = "".join(f"<th>{_e(c)}</th>" for c in head) if head else ""
    body = "".join("<tr>" + "".join(f"<td>{_e(c)}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table>{'<tr>' + h + '</tr>' if h else ''}{body}</table>"


def _fmt(x, nd=2):
    return "n/a" if x is None else f"{x:.{nd}f}"


def page(title, sections):
    return (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>{_e(title)}</title>'
            f"<style>{CSS}</style></head><body><h1>{_e(title)}</h1>{''.join(sections)}</body></html>")


def germline_report(summary, reasons, genes_hit, n_input):
    kinds = summary["by_type"]
    s = [
        '<p class="note">All input here is synthetic or user-supplied; this report describes the file, '
        "it does not interpret variants clinically.</p>",
        "<h2>Filtering</h2>",
        _table([("records in", n_input), ("kept", summary["total"])] +
               [(f"removed: {k}", v) for k, v in sorted(reasons.items())]),
        "<h2>Summary</h2>",
        _table([("SNVs", kinds.get("SNV", 0)), ("insertions", kinds.get("INS", 0)),
                ("deletions", kinds.get("DEL", 0)), ("MNVs", kinds.get("MNV", 0)),
                ("Ti/Tv (SNVs)", _fmt(summary["ti_tv"])),
                ("het / hom-alt", f"{summary['het']} / {summary['hom_alt']} = {_fmt(summary['het_hom_ratio'])}")]),
        "<h2>Variants per chromosome</h2>",
        bar_chart(list(summary["per_chrom"]), list(summary["per_chrom"].values())),
        f"<h2>Annotation</h2><p>{genes_hit} kept variants fall inside the supplied gene intervals; "
        "the rest are labelled intergenic.</p>",
    ]
    return page("Germline variant QC report", s)


def somatic_report(reasons, spectrum, hist, clonal, score, n_input, params):
    bins = len(hist)
    labels = [f"{i / bins:.2f}" for i in range(bins)]
    s = [
        '<p class="note">Somatic filtering of candidate sites from a tumor-normal VCF. This is a filter over '
        "an existing caller's output, not a caller.</p>",
        "<h2>Outcome of every candidate site</h2>",
        _table([("candidates in", n_input)] + sorted(reasons.items(), key=lambda kv: -kv[1])),
        "<h2>Substitution spectrum (called SNVs)</h2>",
        bar_chart(list(spectrum), list(spectrum.values())),
        "<h2>Tumor VAF of called variants</h2>",
        bar_chart(labels, hist, color="#E8A5B8"),
        _table([("VAF histogram mode", _fmt(clonal["mode_vaf"])),
                ("implied purity (2 x mode)", _fmt(clonal["implied_purity"]))]),
        '<p class="note">Implied purity assumes the dominant clone is diploid and copy-neutral at these sites. '
        "Copy-number changes or LOH would bias it.</p>",
    ]
    if score:
        s += ["<h2>Scored against the synthetic truth</h2>",
              _table([("called", score["called"]), ("true positives", score["true_positive"]),
                      ("false positives", score["false_positive"]),
                      ("precision", _fmt(score["precision"])),
                      ("recall, clonal", _fmt(score["recall_clonal"])),
                      ("recall, subclonal", _fmt(score["recall_subclonal"]))])]
    s += ["<h2>Thresholds used</h2>", _table(sorted(params.items()))]
    return page("Somatic (tumor-normal) report", s)
