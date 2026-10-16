#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig. 1 v19 - TWO panels (user decision 2026-09-22): (a) cross-size
prediction with the published exponents (v18 panel c), (b) direction-
independent within-pair fold discrepancy (v18 panel d).
v18 panels (a) published evidence and (b) controlled validation are
REMOVED from the main figure: (a) moves to SI Fig S1 (gap narrative,
manuscript text will reference it), (b) is dropped entirely -- the
exposed size fractions it visualised will be stated in the Methods
text instead (user: "b는 연구방법에 우리가 노출한 사이즈를 이야기하면
바로 해결이야. 굳이 시각화 할필요 없어").
Computation is IDENTICAL to v18 (same seeds, same pairs): for every
pair D_raw = |log10(EC20 / EC45-75)| (measured disagreement, no model)
and D_norm = |log10((EC20*CF20) / (EC45-75*CF45-75))| under each
published exponent; FD = 10**D, so FD = 1 is perfect agreement and the
absolute value makes the metric direction-independent. Normalization
does NOT reduce the discrepancy: median FD rises from x7.9 (measured)
to x30.5 (alpha 1.60) and x17.6 (alpha 2.07).
Panel (a): parity for the 20 -> 45-75 um transition, TWO exponents
(alpha = 1.60 generic, grey, muted; alpha = 2.07 seawater,
blue, emphasized). Colour = exponent, marker shape =
species, fill = acute(filled)/chronic(open). Predictions are defined
from the smaller to the larger fragment-size fraction.
EC_B = EC_A * CF(A)/CF(B), each fraction's ACTUAL bounds entered into
the power law (point fractions use the canonical +/-0.5 um rule).
FE = EC_hat/EC_obs; AFE = max(FE, 1/FE); per-series cluster bootstrap
(BOOT_N=10000, seed=20260913). No R^2.
"""
import sys, math, json, random
import statistics as st
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from track1_fourpanel_applicability import (
    parse_size, parse_num, norm_shape, geomean, csv_reader)
from mp_pipeline_lib import get_powerlaw_cf

MP = str(Path(__file__).resolve().parents[1])
OUTDIR = f"{MP}/output"
OUT = Path(__file__).resolve().parents[1] / "data"
PNG = f"{OUTDIR}/fig1_v23_empirical_test.png"
AUDIT = f"{OUTDIR}/fig1_v23_audit.json"
# v22 (2026-10-15, user): (1) panel (a) parity KEPT as-is except the
#   grid best-fit series (alpha 4.61) is REMOVED from the figure and its
#   legend entry (kept in the audit only: alpha_best_in_figure=False);
#   (2) alpha 2.07 is EMPHASIZED over alpha 1.60: legend lists 2.07
#   first with a larger marker, 2.07 is drawn last on top with full
#   alpha while 1.60 is muted (alpha 0.45); 1.60 relabelled "generic"
#   (Kooi & Koelmans 2019 multi-environment exponent - NOT a past
#   seawater-only exponent) and 2.07 "seawater". Data 1:1 =
#   v21 (same seeds, same pairs).
# (3) panel (b) REPLACED: the direction-independent fold-discrepancy
#   panel is removed and (b) now shows the power-law size
#   distributions p(L) proportional to L**(-alpha) for the two
#   exponents, both normalized on the SAME interval 20-300 um (shaded
#   band). The FD numbers (median x30.5 for alpha 1.60, x17.6 for
#   alpha 2.07) move to the manuscript text: the improvement is
#   reportable, but normalization accuracy must NOT be claimed as
#   validated (x17.6 remains substantial disagreement).
# v21 (2026-10-14): widened the (a)/(b) gap by widening the canvas
#   (W_IN 10.5 -> 10.9 in; square side UNCHANGED 4.2925 in so panel
#   contents and their fixed-pt fonts keep the 확정 v20 relative layout)
#   -- (b)'s y-axis label was overlapping panel (a) (user: "Fig1 ab 간격을
#   벌려줘 b 축제목이 겹쳐 있어"). Data/numbers 1:1 = v20 (same seeds, same
#   pairs).
# v20 (2026-09-22): removed in-panel interpretation text boxes per user --
#   (a) prediction-direction note + within-factor-2/10 block;
#   (b) top description + italic conclusion. Data, geometry, seeds, numbers 1:1 = v19.

ALPHA = 2.07
STD_MAX = 5000.0
C_MEAS = "#3b4a5a"
C_ASSEM = "#8a8a8a"
C_A207 = "#0072B2"    # alpha 2.07 seawater exponent - MARINE BLUE (통일 리디자인 2026-10-15: 전체 구성색 스타일, user 승인)
C_ABST = "#2ca02c"    # best-fit exponent - NOT drawn since v22 (audit only)
C_A16 = "#7a7a7a"     # alpha 1.60 generic exponent - grey (unified a+b, user 2026-10-15)
C_F2 = "#666666"
C_F10 = "#bbbbbb"
MARKERS = ["o", "s", "^", "v", "D", "P", "X", "p", "h", "d"]

# ---------------- data: pooled series + parity pairs ----------------
rows = []
for cc, d in zip(csv_reader(OUT / "Species_Rep_Candidates_MP.csv"),
                 [r for r in csv_reader(OUT / "Domestic_Long.csv")
                  if "MP" in str(r.get("Sheet", ""))]):
    if str(d.get("Particle_Conc_Basis", "")) != "measured":
        continue
    if cc["_invivo"] != "True" or cc["_cdir"] == "<":
        continue
    shp = norm_shape(cc["Shape"])
    if shp not in ("Fragment", "Fiber") or cc["Environment"] != "Seawater":
        continue
    sz = parse_size(cc["Size_um"])
    if sz is None:
        continue
    ev, fp = parse_num(cc["Endpoint_Value"]), parse_num(d.get("Factor_per_mg"))
    if ev is None or not fp:
        continue
    lo, hi = sz[0], sz[1]
    if hi is None:
        hi = STD_MAX
    if lo is None:
        lo = 1.0
    rows.append(dict(sp=cc["Species_Name"], tier=str(cc["Tier"]).strip(),
                     lo=lo, hi=hi, pc=ev * fp, shp=shp))

series_raw = defaultdict(dict)
for r in rows:
    series_raw[(r["sp"], r["tier"])].setdefault((r["lo"], r["hi"]), []).append(r)
series = {k: {b: geomean([x["pc"] for x in v]) for b, v in bandmap.items()}
          for k, bandmap in series_raw.items() if len(bandmap) >= 3}

def band_class(lo, hi):
    """22707 nominal classes; fiber bands are not canonical pair members."""
    if lo >= 1000:
        return None
    if lo == hi and lo <= 25:
        return "small"
    if hi <= 30:
        return "small"
    if lo <= 45 and hi <= 80:
        return "mid"
    if hi <= 200:
        return "large"
    return None

TMAP = {("small", "mid"): "small2mid", ("small", "large"): "small2large",
        ("mid", "large"): "mid2large"}
# v11: keep only the dominant transition 20 -> 45-75 um (user decision)
KEEP_TT = {"small2mid"}

pairs = []
for (sp, tier), vals in series.items():
    bl = sorted(vals.keys())
    for i in range(len(bl)):
        for j in range(i + 1, len(bl)):
            A, B = bl[i], bl[j]
            tt = TMAP.get((band_class(*A), band_class(*B)))
            if tt is None or tt not in KEEP_TT:
                continue
            cfa = get_powerlaw_cf(A[0], A[1], ALPHA, 1.0, STD_MAX)
            cfb = get_powerlaw_cf(B[0], B[1], ALPHA, 1.0, STD_MAX)
            pred = vals[A] * cfa / cfb
            obs = vals[B]
            fe = pred / obs
            pairs.append(dict(sp=sp, tier=tier, tt=tt,
                              alo=A[0], ahi=A[1], blo=B[0], bhi=B[1],
                              ec_a=vals[A],
                              a_raw=f"{A[0]:g}\u2013{A[1]:g}",
                              b_raw=f"{B[0]:g}\u2013{B[1]:g}",
                              pred=pred, obs=obs, fe=fe,
                              afe=max(fe, 1.0 / fe)))

def preds_at(alpha):
    """Predicted concentrations for every pair under a given exponent."""
    out = []
    for p in pairs:
        cfa = get_powerlaw_cf(p["alo"], p["ahi"], alpha, 1.0, STD_MAX)
        cfb = get_powerlaw_cf(p["blo"], p["bhi"], alpha, 1.0, STD_MAX)
        pr = p["ec_a"] * cfa / cfb
        fe = pr / p["obs"]
        out.append(dict(pred=pr, fe=fe, afe=max(fe, 1.0 / fe)))
    return out

def afe_stats(res):
    med = st.median([q["afe"] for q in res])
    f2 = 100.0 * sum(1 for q in res if q["afe"] <= 2.0) / len(res)
    f10 = 100.0 * sum(1 for q in res if q["afe"] <= 10.0) / len(res)
    return med, f2, f10

# slope sweep: best-fit exponent minimises the median AFE over these pairs
# (grid 0.20-8.00; a coarse cap at 4.00 pinned the minimum on the boundary)
GRID = [0.20 + 0.01 * k for k in range(781)
        if abs(0.20 + 0.01 * k - 1.0) > 1e-6]     # skip degenerate alpha=1
ALPHA_B = round(min(GRID, key=lambda a: st.median(
    [q["afe"] for q in preds_at(a)])), 2)

ALPHA_16 = 1.60                    # pre-update manuscript exponent
res_16 = preds_at(ALPHA_16)
res_207 = preds_at(ALPHA)
res_bst = preds_at(ALPHA_B)
med_16, f2_16, f10_16 = afe_stats(res_16)
med_207, f2_207, f10_207 = afe_stats(res_207)
med_bst, f2_bst, f10_bst = afe_stats(res_bst)

random.seed(20260913)
keys = sorted(series.keys())
def ci(v):
    s = sorted(v)
    return s[250], s[9750]
def boot_ci(res):
    by_s = {k: [q for q, p in zip(res, pairs)
                if (p["sp"], p["tier"]) == k] for k in keys}
    bm, b2, b10 = [], [], []
    for _ in range(10000):
        smp = []
        for _ in range(len(keys)):
            smp += by_s[random.choice(keys)]
        bm.append(st.median([q["afe"] for q in smp]))
        b2.append(100.0 * sum(1 for q in smp if q["afe"] <= 2.0) / len(smp))
        b10.append(100.0 * sum(1 for q in smp if q["afe"] <= 10.0) / len(smp))
    return ci(bm), ci(b2), ci(b10)
ci_16 = boot_ci(res_16)
ci_207 = boot_ci(res_207)
ci_bst = boot_ci(res_bst)

fe_med_207 = st.median([q["fe"] for q in res_207])
fe_med_bst = st.median([q["fe"] for q in res_bst])
print(f"(a) series={len(series)} pairs={len(pairs)}")
print(f"    alpha=1.60 (pre-update): median AFE x{med_16:.2f} "
      f"f2 {f2_16:.1f}% f10 {f10_16:.1f}%")
print(f"    alpha=2.07: median AFE x{med_207:.2f} "
      f"f2 {f2_207:.1f}% f10 {f10_207:.1f}% median FE x{fe_med_207:.2f}")
print(f"      CI AFE ({ci_207[0][0]:.1f},{ci_207[0][1]:.1f}) "
      f"f2 ({ci_207[1][0]:.1f},{ci_207[1][1]:.1f}) "
      f"f10 ({ci_207[2][0]:.1f},{ci_207[2][1]:.1f})")
print(f"    alpha={ALPHA_B:.2f} (best fit): median AFE x{med_bst:.2f} "
      f"f2 {f2_bst:.1f}% f10 {f10_bst:.1f}% median FE x{fe_med_bst:.2f}")
print(f"      CI AFE ({ci_bst[0][0]:.1f},{ci_bst[0][1]:.1f}) "
      f"f2 ({ci_bst[1][0]:.1f},{ci_bst[1][1]:.1f}) "
      f"f10 ({ci_bst[2][0]:.1f},{ci_bst[2][1]:.1f})")

assert len(series) == 16, len(series)
assert len(pairs) == 18, len(pairs)          # v11: 20->45-75 only
assert all(p["tt"] == "small2mid" for p in pairs)

# ---------------- figure ----------------
# v21 geometry: a single row of TWO 1:1 squares. Panel size = the v18
# parity square (4.2925 in, unchanged from 확정 v20); the extra gap comes
# from a wider canvas: W 10.9 in: left margin 0.9975 in (= v18 0.095*10.5)
# + 2 x 4.2925 squares + 0.95 gap => right margin 0.368 in (= v18 0.035
# margin). v20's 0.55 in gap was too narrow for (b)'s y-axis label, which
# overlapped panel (a).
# figsize (10.9, 5.95): top 0.844 in (panel labels + title + stats text,
# label baseline sits 0.115*SQ above the axes top) + square 4.2925 +
# bottom 0.814 in (x ticks + xlabel).
W_IN, H_IN = 10.9, 5.95
SQ = 4.2925                                  # square side, inches (= v18/v20)
SQ_W = SQ / W_IN                             # axes frac
SQ_H = SQ / H_IN
X0 = 0.9975 / W_IN                           # left margin 0.9975 in (= v18)
SQGAP = 0.95 / W_IN                          # widened (v21, user 2026-10-14)
Y0 = 0.814 / H_IN                            # bottom margin (xlabel+ticks)
fig = plt.figure(figsize=(W_IN, H_IN))
axp = fig.add_axes([X0, Y0, SQ_W, SQ_H])     # (a) published-exponent parity
axd = fig.add_axes([X0 + SQ_W + SQGAP, Y0, SQ_W, SQ_H])  # (b) FD (v18 d)
assert axd is not axp
axp.set_aspect("equal", adjustable="box")    # log decades equal
_p = axp.get_position()
assert abs(_p.width * W_IN - _p.height * H_IN) < 1e-9, "1:1 square"

def panel_label(ax, letter, lx=-0.075):
    ax.text(lx, 1.115, letter, transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="bottom", ha="left")

# ---- (a) parity FORWARD: measured 20 -> predicted 45-75 (22707, v16) ----
sps_b = sorted({p["sp"] for p in pairs})
pts_b = []
# v22: alpha 1.60 (generic) drawn first and muted; alpha 2.07 (seawater)
# drawn last on top with larger markers - 2.07 is emphasized.
for res, col, msz, alp, zr in ((res_16, C_A16, 5.4, 0.45, 4),
                               (res_207, C_A207, 6.6, 0.90, 6)):
    for p, q in zip(pairs, res):
        acute = "chronic" not in p["tier"].lower()
        mfc = col if acute else "white"
        axp.plot(p["obs"], q["pred"], ls="none",
                 marker=MARKERS[sps_b.index(p["sp"])],
                 ms=msz if acute else msz + 0.4, color=col, mfc=mfc,
                 mec=col, mew=1.2, alpha=alp, zorder=zr)
        pts_b.append((p["obs"], q["pred"]))

lims = [1e2, 1e12]
axp.plot(lims, lims, color="black", lw=1.0, zorder=2)
axp.plot(lims, [l * 2 for l in lims], color=C_F2, lw=0.8, ls=(0, (4, 3)),
         zorder=1)
axp.plot(lims, [l / 2 for l in lims], color=C_F2, lw=0.8, ls=(0, (4, 3)),
         zorder=1)
axp.plot(lims, [l * 10 for l in lims], color=C_F10, lw=0.8, ls=(0, (2, 3)),
         zorder=1)
axp.plot(lims, [l / 10 for l in lims], color=C_F10, lw=0.8, ls=(0, (2, 3)),
         zorder=1)
axp.set_xscale("log"); axp.set_yscale("log")
axp.set_xlim(*lims); axp.set_ylim(*lims)
axp.set_xlabel("Observed effect concentration (particles L$^{-1}$, log)",
               fontsize=9.5)
axp.set_ylabel("Predicted effect concentration (particles L$^{-1}$, log)",
               fontsize=9.5)
panel_label(axp, "(a)")
axp.text(0.015, 0.968,
         f"Median abs. fold error ×{med_16:.1f} / ×{med_207:.1f}"
         f"\n(α = 1.60 / 2.07)",
         transform=axp.transAxes, fontsize=7.6, va="top", color="black",
         fontweight="bold")

# legend: exponent colours (with n) + species markers + acute/chronic fill
handles_a = []
# 2.07 listed FIRST with a larger marker - emphasized over 1.60 (v22)
handles_a.append(Line2D([0], [0], color=C_A207, marker="o", ls="none",
                        ms=6.0,
                        label=f"α = 2.07 (seawater, n={len(pairs)})"))
handles_a.append(Line2D([0], [0], color=C_A16, marker="o", ls="none",
                        ms=5.0,
                        label=f"α = 1.60 (generic, n={len(pairs)})"))
sp_short = {s: s.split()[0] + " " + s.split()[1][:4] + "." for s in sps_b}
for i, s in enumerate(sps_b):
    handles_a.append(Line2D([0], [0], color="#444444", marker=MARKERS[i],
                            ls="none", ms=5.0, label=sp_short[s]))
handles_a.append(Line2D([0], [0], color="#444444", marker="o", ls="none",
                        ms=5.5, mfc="#444444", label="Acute (filled)"))
handles_a.append(Line2D([0], [0], color="#444444", marker="o", ls="none",
                        ms=5.5, mfc="white", label="Chronic (open)"))
legp = axp.legend(handles=handles_a, loc="lower right", fontsize=6.4,
                  frameon=False, handlelength=1.2, labelspacing=0.34,
                  ncol=2, columnspacing=0.9)

# ---- (b) power-law size distributions (v22; replaces the FD panel) ----
# p(L) proportional to L**(-alpha) for the two exponents; BOTH curves
# are normalized on the SAME interval 20-300 um (the standardization
# target range, shaded) so the comparison is like-for-like. Within the
# band the two models are consistent by construction; outside it they
# diverge - this is what propagates into panel (a)'s exponent-dependent
# correction factors. Median FD numbers (x30.5 for alpha 1.60, x17.6
# for alpha 2.07) move to the manuscript text: reportable, but
# normalization accuracy must NOT be claimed as validated.
L1, L2 = 20.0, 300.0                 # normalization interval (um)
XMIN, XMAX = 1.0, 5000.0             # um; XMAX = pipeline STD_MAX cap

def plaw_c(a):
    """Normalizing constant: integral of C*L**(-a) over [L1, L2] = 1."""
    return (1.0 - a) / (L2 ** (1.0 - a) - L1 ** (1.0 - a))

c_16 = plaw_c(ALPHA_16)
c_207 = plaw_c(ALPHA)
_nlog = 2400
grid = [10.0 ** (math.log10(XMIN) + k * (math.log10(XMAX)
        - math.log10(XMIN)) / _nlog) for k in range(_nlog + 1)]

def integ_plaw(a, c):
    """Exact per-cell integral of c*L**(-a) over [L1, L2] (clipped bounds)."""
    tot = 0.0
    for i in range(len(grid) - 1):
        lo, hi = grid[i], grid[i + 1]
        if hi <= L1 or lo >= L2:
            continue
        lo2, hi2 = max(lo, L1), min(hi, L2)
        tot += c * (hi2 ** (1.0 - a) - lo2 ** (1.0 - a)) / (1.0 - a)
    return tot

int_16 = integ_plaw(ALPHA_16, c_16)
int_207 = integ_plaw(ALPHA, c_207)
assert abs(int_16 - 1.0) < 1e-9, int_16
assert abs(int_207 - 1.0) < 1e-9, int_207
print(f"(b) p(L)~L^-alpha normalized on {L1:g}-{L2:g} um: "
      f"C(1.60)={c_16:.4f}, C(2.07)={c_207:.4f}; numeric integrals "
      f"{int_16:.8f} / {int_207:.8f}")
axd.axvspan(L1, L2, color="#e8e8e8", zorder=0)
ys_16 = [c_16 * L ** (-ALPHA_16) for L in grid]
ys_207 = [c_207 * L ** (-ALPHA) for L in grid]
axd.plot(grid, ys_16, color=C_A16, lw=1.8, zorder=3,
         label="α = 1.60 (generic)")
axd.plot(grid, ys_207, color=C_A207, lw=2.0, zorder=4,
         label="α = 2.07 (seawater)")
axd.set_xscale("log")
axd.set_yscale("log")
axd.set_xlim(XMIN, XMAX)
axd.set_ylim(1e-7, 60.0)
axd.set_xlabel("Particle size (μm, log)", fontsize=9.5)
axd.set_ylabel("Relative particle abundance (log)", fontsize=9.5)
legd = axd.legend(loc="lower left", fontsize=7.4, frameon=False)
axd.text(math.sqrt(L1 * L2), 8.0, "normalization range\n20–300 μm",
         fontsize=7.2, color="#444444", ha="center", va="bottom", zorder=5)
panel_label(axd, "(b)", lx=-0.088)
pts_d = list(zip(grid, ys_16)) + list(zip(grid, ys_207))

# ---------------- QA ----------------
fig.canvas.draw()
rend = fig.canvas.get_renderer()

def colliding_texts(ax):
    out = []
    verts = []
    for ln in ax.lines:
        xy = ln.get_xydata()
        if len(xy):
            verts.extend(ax.transData.transform(xy))
    for coll in ax.collections:
        try:
            segs = coll.get_segments()
        except Exception:
            segs = []
        for s in segs:
            if len(s):
                verts.extend(ax.transData.transform(s))
    for t in ax.texts:
        s = t.get_text().strip()
        if not s:
            continue
        tb = t.get_window_extent(rend)
        inside = sum(1 for (px, py) in verts
                     if tb.x0 <= px <= tb.x1 and tb.y0 <= py <= tb.y1)
        if inside:
            out.append((s[:40], inside))
    return out

def legend_collisions(leg, ax, pts):
    bb = leg.get_window_extent()
    n = 0
    for x, y in pts:
        dx, dy = ax.transData.transform((x, y))
        if bb.x0 <= dx <= bb.x1 and bb.y0 <= dy <= bb.y1:
            n += 1
    return n

qa_text = dict(axa=colliding_texts(axp), axd=colliding_texts(axd))
qa_leg = dict(panel_a=legend_collisions(legp, axp, pts_b),
              panel_b=legend_collisions(legd, axd, pts_d))
qa_misc = dict(pd_gap_overlap=bool(axp.get_position().overlaps(
    axd.get_position())))
print("QA text-vs-data:", qa_text)
print("QA legend-vs-markers:", qa_leg)
print("QA misc:", qa_misc)
assert not any(qa_text.values()), qa_text
_rend2 = legd.figure.canvas.get_renderer()
_legx = {"panel_a": [t.get_text()[:24] for t in axp.texts if
                     legp.get_window_extent(_rend2).overlaps(
                         t.get_window_extent(_rend2))],
         "panel_b": [t.get_text()[:24] for t in axd.texts if
                     legd.get_window_extent(_rend2).overlaps(
                         t.get_window_extent(_rend2))]}
print("QA legend-vs-TEXT:", _legx)
assert all(v == [] for v in _legx.values()), _legx
assert all(v == 0 for v in qa_leg.values()), qa_leg
assert not any(qa_misc.values()), qa_misc
# both panels must be exact 1:1 squares
for _nm, _ax in (("a", axp), ("b", axd)):
    _pb = _ax.get_position()
    _wpx = _pb.width * fig.bbox.width
    _hpx = _pb.height * fig.bbox.height
    print("QA panel-%s box px: %.1f x %.1f  ratio=%.5f"
          % (_nm, _wpx, _hpx, _wpx / _hpx))
    assert abs(_wpx / _hpx - 1.0) < 0.01, (_nm, _wpx, _hpx)
# panel label / top texts must stay inside the canvas
for _nm, _ax in (("a", axp), ("b", axd)):
    for t in _ax.texts:
        tb = t.get_window_extent(rend)
        assert tb.y1 <= fig.bbox.height + 1.0, (_nm, t.get_text()[:30])
        assert tb.x1 <= fig.bbox.width + 1.0, (_nm, t.get_text()[:30])
# stats texts must not overlap each other (text-vs-text, both panels)
def _tvs(boxes):
    bb = [b for b in boxes if b is not None]
    for i in range(len(bb)):
        for j in range(i + 1, len(bb)):
            assert not bb[i].overlaps(bb[j]), (i, j,
                                               bb[i].bounds, bb[j].bounds)
_tvs([t.get_window_extent(rend) for t in axp.texts if t.get_text().strip()])
_tvs([t.get_window_extent(rend) for t in axd.texts if t.get_text().strip()])
print("QA text-vs-text overlaps: 0")
# x axis labels must stay above the canvas bottom
for _nm, _ax in (("a", axp), ("b", axd)):
    _xb = _ax.xaxis.label.get_window_extent(rend)
    assert _xb.y0 >= 2.0, (_nm, _xb.y0)
print("QA xlabel bottom margins OK")
# v21: (b)'s ylabel must clear panel (a)'s right spine
_yb = axd.yaxis.label.get_window_extent(rend)
_pa = axp.get_window_extent(rend)
assert _yb.x0 >= _pa.x1, ("(b) ylabel overlaps panel (a)", _yb.x0, _pa.x1)
print("QA (b) ylabel clears panel (a): x0 %.1f >= %0.1f (+%.2f px)"
      % (_yb.x0, _pa.x1, _yb.x0 - _pa.x1))
# (b) integrity: both curves integrate to 1 over the SAME 20-300 um
# window (asserted at build time, exact per-cell quadrature eps
# 1e-9); the grey/blue curves cross inside the band and diverge outside it, which is what
# drives panel (a)'s exponent-dependent correction factors.
assert c_16 > 0 and c_207 > 0
print("QA (b) same-interval normalization: 20-300 um integrals = 1 "
      "(%.6f / %.6f); curves cross at %.1f um"
      % (int_16, int_207, (c_207 / c_16) ** (1.0 / (ALPHA - ALPHA_16))))

fig.savefig(PNG, dpi=300)
plt.close(fig)

audit = dict(
    panel_a=dict(parity=True, direction="forward 20 -> 45-75 um",
                 n_pairs=len(pairs),
                 alpha=ALPHA, median_afe=med_207, pct_f2=f2_207,
                 pct_f10=f10_207, ci_afe=ci_207[0], ci_f2=ci_207[1],
                 ci_f10=ci_207[2],
                 alpha_16=ALPHA_16, median_afe_16=med_16,
                 pct_f2_16=f2_16, pct_f10_16=f10_16,
                 ci_afe_16=ci_16[0], ci_f2_16=ci_16[1],
                 ci_f10_16=ci_16[2],
                 alpha_best=ALPHA_B, median_afe_best=med_bst,
                 pct_f2_best=f2_bst, pct_f10_best=f10_bst,
                 ci_afe_best=ci_bst[0], ci_f2_best=ci_bst[1],
                 ci_f10_best=ci_bst[2],
                 alpha_best_in_figure=False,
                 emphasis_note=("v22: alpha 2.07 (seawater) emphasized - "
                                "legend first, drawn last on top; "
                                "alpha 1.60 (generic) muted; colours unified "
                                "across panels: 2.07 blue / 1.60 grey; "
                                "legend terms unified: 2.07 seawater / "
                                "1.60 generic (user 2026-10-15)"),
                 pair_list=pairs),
    panel_b=dict(size_distributions=True,
                 model="p(L) proportional to L**(-alpha)",
                 alphas=dict(alpha_16=ALPHA_16, alpha_207=ALPHA),
                 labels=dict(alpha_16="generic",
                             alpha_207="seawater"),
                 normalization_interval_um=[L1, L2],
                 norm_constant_alpha_16=c_16,
                 norm_constant_alpha_207=c_207,
                 numeric_integral_alpha_16=int_16,
                 numeric_integral_alpha_207=int_207,
                 x_range_um=[XMIN, XMAX],
                 band_shaded=True,
                 curves_cross_at_um=(c_207 / c_16)
                 ** (1.0 / (ALPHA - ALPHA_16)),
                 note=("v22 replaces the FD panel; both curves "
                       "normalized on the SAME 20-300 um interval. "
                       "Median FD x30.5 (alpha 1.60) / x17.6 (alpha "
                       "2.07) reportable in text; normalization "
                       "accuracy must NOT be claimed as validated.")),
    provenance=dict(
        removed_from_main=dict(
            published_evidence="moved to SI Fig S1 (user 2026-09-22)",
            controlled_validation=("dropped; exposed size fractions to be "
                                   "stated in Methods text")),
        supersedes=("fig1_empirical_test_v21.py (v22: best-fit series "
                    "removed; 2.07 emphasized; panel (b) replaced by "
                    "power-law size distributions)")),
)
with open(AUDIT, "w") as f:
    json.dump(audit, f, ensure_ascii=False, indent=1)
print("saved", PNG)