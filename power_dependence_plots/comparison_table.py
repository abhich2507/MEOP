"""
Interactive side-by-side comparison table:
  Table A  — table.tsv           (reference / previous analysis, uses σsys)
  Table B  — C8_C9_summary.csv  (current notebook analysis, uses σrms)

For each (cell, power, parameter):
  - Two colour-coded sub-cells side by side: [Ref]  [New]
  - Δ = New − Ref shown as a coloured badge
  - Hover tooltip shows all error components
"""

import pandas as pd
import numpy as np

# ── load Table B ─────────────────────────────────────────────────────────────
dfB = pd.read_csv("/Users/snip/Documents/MEOP/power_dependence_plots/C8_C9_summary.csv")
dfB["n"]    = dfB["n"].astype(int)
dfB["cell"] = dfB["cell"].str.upper()
B = dfB.set_index(["cell", "power_W"])

# ── load Table A ─────────────────────────────────────────────────────────────
raw = pd.read_csv("/Users/snip/Documents/MEOP/table.tsv", sep="\t", header=0)
raw = raw.dropna(subset=[raw.columns[0]])
raw = raw[raw.iloc[:, 0].astype(str).str.strip().str.match(r"^\d")]
raw["power_W"] = raw.iloc[:, 0].astype(int)
raw["cell"]    = raw.iloc[:, 1].str.strip().str.upper()
c = list(raw.columns)
raw = raw.rename(columns={
    c[2]:"P_max_mean",  c[3]:"P_max_ssys",  c[4]:"P_max_sfit",  c[5]:"P_max_stot",
    c[6]:"tau_d_mean",  c[7]:"tau_d_ssys",  c[8]:"tau_d_sfit",  c[9]:"tau_d_stot",
    c[10]:"tau_b_mean", c[11]:"tau_b_ssys", c[12]:"tau_b_sfit", c[13]:"tau_b_stot",
})
for col in ["P_max_mean","P_max_stot","tau_d_mean","tau_d_stot","tau_b_mean","tau_b_stot"]:
    raw[col] = pd.to_numeric(raw[col], errors="coerce")
A = raw.set_index(["cell", "power_W"])

powers = sorted(dfB["power_W"].unique())
cells  = ["C8", "C9"]

# ── colour helpers ────────────────────────────────────────────────────────────
def lerp_rgb(val, lo, hi, c0, c1):
    t = float(np.clip((val - lo) / (hi - lo + 1e-12), 0, 1))
    return tuple(int(c0[i] + t * (c1[i] - c0[i])) for i in range(3))

RAMPS = {
    "P_max": ((210,230,255), (0,70,200)),
    "tau_b": ((255,240,210), (190,70,0)),
    "tau_d": ((210,255,220), (0,130,40)),
}
RANGES = {}
for p, Bcol, Acol in [("P_max","P_max_mean","P_max_mean"),
                       ("tau_b","tau_b_mean","tau_b_mean"),
                       ("tau_d","tau_decay_mean","tau_d_mean")]:
    vB = dfB[Bcol].dropna().values
    vA = raw[Acol].dropna().values
    RANGES[p] = (min(vB.min(), vA.min()), max(vB.max(), vA.max()))

def bg_fg(param, val):
    lo, hi = RANGES[param]
    c0, c1 = RAMPS[param]
    r,g,b = lerp_rgb(val, lo, hi, c0, c1)
    lum = 0.299*r + 0.587*g + 0.114*b
    return f"rgb({r},{g},{b})", ("#111" if lum > 140 else "#fff")

def delta_badge(d):
    if abs(d) < 0.005:  col, s = "#888", "≈0"
    elif d > 0:         col, s = "#2ecc71", f"+{d:.3f}"
    else:               col, s = "#e74c3c", f"{d:.3f}"
    return f'<span class="delta" style="color:{col}">Δ {s}</span>'

def pair_cell(param, cell, pw):
    rA = A.loc[(cell,pw)] if (cell,pw) in A.index else None
    rB = B.loc[(cell,pw)] if (cell,pw) in B.index else None

    Bmean = "tau_decay_mean" if param=="tau_d" else f"{param}_mean"
    Bstot = "tau_decay_stot" if param=="tau_d" else f"{param}_stot"
    Bsfit = "tau_decay_sfit" if param=="tau_d" else f"{param}_sfit"
    Bsrms = "tau_decay_srms" if param=="tau_d" else f"{param}_srms"

    vA = float(rA[f"{param}_mean"]) if rA is not None else None
    eA = float(rA[f"{param}_stot"]) if rA is not None else None
    vB = float(rB[Bmean]) if rB is not None else None
    eB = float(rB[Bstot]) if rB is not None else None

    if vA is not None:
        bgA, fgA = bg_fg(param, vA)
        tipA = (f"[LIND] mean={vA:.4f}  σsys={float(rA[f'{param}_ssys']):.4f}"
                f"  σfit={float(rA[f'{param}_sfit']):.5f}  σtot={eA:.4f}")
        cA = (f'<div class="subcell" style="background:{bgA};color:{fgA}" title="{tipA}">'
              f'<span class="src">LIND</span>'
              f'<span class="mean">{vA:.3f}</span>'
              f'<span class="err">±{eA:.3f}</span></div>')
    else:
        cA = '<div class="subcell missing">—</div>'

    if vB is not None:
        bgB, fgB = bg_fg(param, vB)
        tipB = (f"[ABHI] mean={vB:.4f}  σ̄fit={float(rB[Bsfit]):.4f}"
                f"  σrms={float(rB[Bsrms]):.4f}  σtot={eB:.4f}")
        cB = (f'<div class="subcell" style="background:{bgB};color:{fgB}" title="{tipB}">'
              f'<span class="src">ABHI</span>'
              f'<span class="mean">{vB:.3f}</span>'
              f'<span class="err">±{eB:.3f}</span></div>')
    else:
        cB = '<div class="subcell missing">—</div>'

    badge = delta_badge(vB - vA) if (vA is not None and vB is not None) else ""
    return f'<td class="pcell">{cA}{cB}{badge}</td>'

# ── build rows ────────────────────────────────────────────────────────────────
rows_html = []
for cell in cells:
    for pw in powers:
        rB = B.loc[(cell,pw)] if (cell,pw) in B.index else None
        n_str = str(int(rB["n"])) if rB is not None else "?"
        cls = "c8" if cell=="C8" else "c9"
        row = (f'<tr>'
               f'<td class="cell-lbl {cls}">{cell}</td>'
               f'<td class="pow">{pw} W</td>'
               + pair_cell("P_max", cell, pw)
               + pair_cell("tau_b", cell, pw)
               + pair_cell("tau_d", cell, pw)
               + f'<td class="n">{n_str}</td>'
               f'</tr>')
        rows_html.append(row)
    rows_html.append('<tr class="sep"><td colspan="7"></td></tr>')

# ── HTML template ─────────────────────────────────────────────────────────────
html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>MEOP — Reference vs New Analysis</title>
<style>
* { box-sizing:border-box; margin:0; padding:0; }
body { font-family:"Segoe UI",Arial,sans-serif; background:#12121f; color:#ddd;
       display:flex; flex-direction:column; align-items:center; padding:30px 12px; }
h1   { font-size:1.45rem; margin-bottom:4px; letter-spacing:1px; color:#eee; }
.subtitle { font-size:0.82rem; color:#888; margin-bottom:24px; }
.wrapper { overflow-x:auto; width:100%; max-width:1300px; }
table { border-collapse:collapse; width:100%; font-size:0.8rem; }
thead th { background:#0d2137; padding:8px 6px; border-bottom:2px solid #e94560; letter-spacing:0.4px; }
thead tr:nth-child(2) th { background:#0a1a2e; padding:5px 6px; color:#bbb; font-weight:400;
                            border-bottom:2px solid #1e1e3a; font-size:0.73rem; }
td, th { border:1px solid #1e1e3a; text-align:center; vertical-align:middle; padding:3px 5px; }
.cell-lbl { font-weight:700; font-size:0.85rem; width:36px; }
.c8 { color:#ff7675; background:#1a0f0f; }
.c9 { color:#74b9ff; background:#0a0f1a; }
.pow { color:#fdcb6e; background:#1a1500; font-weight:600; white-space:nowrap; width:52px; }
.n   { color:#666; background:#111; font-size:0.73rem; width:28px; }
.pcell { padding:4px 5px; min-width:180px; }
.subcell { display:inline-flex; flex-direction:column; align-items:center;
           border-radius:5px; padding:4px 8px; margin:2px; min-width:72px;
           cursor:default; transition:filter 0.12s; vertical-align:top; }
.subcell:hover { filter:brightness(1.3); outline:2px solid #fff; }
.missing { background:#1e1e1e!important; color:#555!important; border:1px dashed #333; }
.src  { font-size:0.6rem; opacity:0.65; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:1px; }
.mean { font-weight:700; font-size:0.85rem; }
.err  { font-size:0.68rem; opacity:0.82; }
.delta { display:block; font-size:0.68rem; margin-top:4px; font-weight:700; letter-spacing:0.3px; }
.sep td { background:#1e1e3a; height:5px; border:none; }
.hdr-pmax { color:#ff7675; }
.hdr-taub { color:#fdcb6e; }
.hdr-taud { color:#55efc4; }
.legend { margin-top:20px; font-size:0.75rem; color:#666; text-align:center;
          line-height:1.9; max-width:900px; }
.legend b { color:#999; }
.pill { display:inline-block; border-radius:4px; padding:1px 7px; font-size:0.7rem;
        font-weight:600; margin:0 2px; }
</style>
</head>
<body>
<h1>MEOP Fit Parameters — LIND vs ABHI</h1>
<p class="subtitle">
  <span class="pill" style="background:#1a2a1a;color:#2ecc71">LIND</span> table.tsv &nbsp;·&nbsp;
  <span class="pill" style="background:#1a1a2e;color:#74b9ff">ABHI</span> C8_C9_summary.csv &nbsp;·&nbsp;
  Δ = ABHI − LIND
</p>
<div class="wrapper">
<table>
  <thead>
    <tr>
      <th rowspan="2" colspan="2">Cell / Power</th>
      <th class="hdr-pmax">P<sub>max</sub> [%]</th>
      <th class="hdr-taub">τ<sub>b</sub> [s]</th>
      <th class="hdr-taud">τ<sub>decay</sub> [s]</th>
      <th rowspan="2" style="color:#555">n</th>
    </tr>
    <tr>
      <th>mean ± σ<sub>tot</sub> &nbsp;|&nbsp; Δ(ABHI−LIND)</th>
      <th>mean ± σ<sub>tot</sub> &nbsp;|&nbsp; Δ(ABHI−LIND)</th>
      <th>mean ± σ<sub>tot</sub> &nbsp;|&nbsp; Δ(ABHI−LIND)</th>
    </tr>
  </thead>
  <tbody>
    ROWS_PLACEHOLDER
  </tbody>
</table>
</div>
<div class="legend">
  <b>LIND</b> σ<sub>tot</sub> = √(σ<sub>sys</sub>² + σ<sub>fit</sub>²) &nbsp;·&nbsp;
  <b>ABHI</b> σ<sub>tot</sub> = √(σ̄<sub>fit</sub>² + σ<sub>rms</sub>²) &nbsp;·&nbsp;
  Colour intensity encodes the mean value on a shared scale per parameter column
</div>
</body>
</html>
"""

html = html.replace("ROWS_PLACEHOLDER", "\n    ".join(rows_html))
out = "/Users/snip/Documents/MEOP/power_dependence_plots/comparison_table.html"
with open(out, "w") as f:
    f.write(html)
print(f"Saved → {out}")
import subprocess
subprocess.Popen(["open", out])
