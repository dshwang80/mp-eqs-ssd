"""
Track 1 four-panel applicability review figure (v7, seawater-only + shape strata).

Design doc: mp-eqs-derivation skill, references/applicability-figure-design.md
(user-approved 2026-09-02, incl. three review amendments; amended 2026-09-03:
 seawater-only validation + shape strata, per user directive).

User directive (2026-09-03): validate SEAWATER ONLY. Sediment is excluded from
the validation panels: (1) sediment has its own power-law exponent (Kooi 2021:
2.57 +/- 0.20 marine sediment vs 2.07 +/- 0.03 marine surface water) - applying
the seawater slope to sediment points is an error; (2) exposure scenarios differ
strongly. Sediment appears only as the Kooi alpha-span citation in (c). Shape
strata (fragment vs fiber) added per user directive: the seawater evidence had
no shape contrast; fibers are shown as a separate stratum (Track 1 criteria
scoping remains fragment-based; here fibers are a diagnostic stratum, not pooled
with fragments).

Question: does power-law size standardization (Track 1, gamma = 1 - alpha,
1-5000 um) agree with the domestic native size-band evidence in seawater, and
does the answer differ by particle shape?

Panels (v6, 2026-09-08 user verdict on the v5 render - three defects fixed)
  (a) Raw seawater band-level endpoints vs band size: count currency (left),
      mass currency (right). Markers: fragment = circle, fiber = triangle.
      Within-cell link lines (same species x test method, only the size band
      differs) in muted species colors under the markers; legend in the empty
      400-2000 um column (framed). v6 FIX: cell_links stored particle-count
      levels as x, so (a) drew phantom lines at x >= 1e4 and no lines across
      the data region; links now connect band sizes within each cell.
  (b) v7 STORY PANEL (user principle 2026-09-09: an analysis figure must
      explain the claim, not exhibit a diagnostic): pooled observed
      within-cell size slopes, ln(count) vs ln(size-mid), cell-centered,
      with 95% CI, against the assumed single-alpha slope -(alpha-1) =
      -1.07 (dashed) and the flat zero line (dotted). Fragment is flat
      (all bands +0.11 +/- 0.34; <=63 um bands -0.43 +/- 0.41), far from
      -1.07; 63-333 um bands slightly positive (+0.78); fiber steeper
      (-1.33 +/- 0.23, near the assumed slope). Readout: NO observed
      stratum reproduces the assumed slope -> a single alpha = 2.07
      standardization is not supported by the data; this is the empirical
      basis for the (c) band-dependence and the <56 um extrapolation
      limit. v5/v6 standardized-point scatter and its stats box REMOVED;
      beta diagnostics stay in the manuscript/tables.
  (c) Conversion-factor uncertainty BY SOURCE BAND (pure method quantities,
      no data points): CF(alpha) curves for 333-5000 um (adopted band,
      black), 56-333 um (blue), and 20-56 / 5-20 um (dashed = source band
      below the fitted validity lower bound 56 +/- 21 um, Kooi 2021 ->
      extrapolation). Adopted alpha = 2.07 crimson vline with CI whisker
      (x1.4); old default 1.6 gray dashed vline with Koelmans +-10% whisker
      (x~2); alpha-update arrow (CF x13.1) on the adopted band. The alpha
      update changes CF by a different multiple per source band: x13.1
      (333-5000) > x5.2 (56-333) > x2.8 (20-56) > x1.6 (5-20), and sources
      below 56 um are standardization-limited (alpha unfitted there).
  (d) Direct pair test, seawater cells with >= 3 count bands (16 fragment +
      1 fiber): required ln(CF_L/CF_S) vs measured ln(PC_S/PC_L); 1:1 line,
      median-offset guide, residual rug in the right margin.

Conventions (probe2-verified, mp_pipeline_lib / track3 1:1)
  - Pool: measured basis, in vivo, '<' excluded, parseable size, positive
    endpoint -> 490 rows (seawater 362: fragment 245 + fiber 117; bead 0).
  - Cells: species x environment x tier x shape. Bands keyed by size mid
    (track3). CF: get_powerlaw_cf(lo, hi, alpha). Open bounds ('>150') and
    single-value sizes use the narrow +-0.5 um proxy (reproduces the verified
    pair requirement +4.17 for 20 vs '>150' at alpha = 2.07).
  - alpha = 2.07 seawater everywhere (media_alpha retained for provenance).

Verification anchors (asserted on exec)
  pool 490 (SW 362); fragment cells 18 (multi 16, valid 16, bands 50);
  fiber cells 15 (multi 14, valid 1, bands 30);
  beta_count: fragment +0.1074 +/- 0.3427 (n=48), fiber -1.3293 +/- 0.2336
  (n=29); beta_mass: fragment +0.2879 +/- 0.2589, fiber +0.2802 +/- 0.2193;
  (b) points 77; (d) pairs 17; Koelmans bracket x2.34 / x2.24.
All figure text in English (MDPI Toxins).
"""

import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PY_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PY_DIR))
from mp_pipeline_lib import get_powerlaw_cf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "track1_selfrun"
FIG_DIR = ROOT / "03_figures" / "output"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ALPHA_SW, ALPHA_SED, ALPHA_DEF = 2.07, 2.57, 1.6

C_CNT, C_MASS, C_PAIR = "#b2182b", "#2166ac", "#1b7837"
C_FRAG, C_FIB = "#b2182b", "#2166ac"
SHAPES_MK = (("Fragment", "o"), ("Fiber", "^"))
GRAY = "#666666"

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "legend.fontsize": 7, "axes.linewidth": 0.8,
})


# ------------------------------------------------------- pool rebuild (1:1)
def parse_size(s):
    if s is None:
        return None
    t = str(s).strip().replace("\u2013", "-").replace("\u2212", "-")
    if t in ("", "nan", "None", "NA"):
        return None
    t = t.replace("\u03bc", "")
    m = re.match(r"^>\s*([\d.]+)$", t)
    if m:
        v = float(m.group(1))
        return (v, None, None)
    m = re.match(r"^<\s*([\d.]+)$", t)
    if m:
        v = float(m.group(1))
        return (None, v, v / math.sqrt(2))
    m = re.match(r"^([\d.]+)\s*-\s*([\d.]+)$", t)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        return (lo, hi, math.sqrt(lo * hi))
    m = re.match(r"^([\d.]+)$", t)
    if m:
        v = float(m.group(1))
        return (v, v, v)
    return None


def norm_shape(s):
    t = str(s).strip().lower()
    if t.startswith("fr"):
        return "Fragment"
    if t.startswith("f"):
        return "Fiber"
    if t.startswith("b"):
        return "Bead"
    return t.title()


def parse_num(s):
    try:
        v = float(str(s).replace(",", ""))
        return v if v > 0 else None
    except Exception:
        return None


def geomean(vals):
    v = [x for x in vals if x > 0]
    return math.exp(sum(math.log(x) for x in v) / len(v)) if v else None


def build_pool():
    """Measured, in vivo, '<' excluded, fibers KEPT (shape stratum), bead empty.
    490 rows total; seawater subset 362 (fragment 245 + fiber 117)."""
    cand = list(csv_reader(OUT / "Species_Rep_Candidates_MP.csv"))
    dl = list(csv_reader(OUT / "Domestic_Long.csv"))
    dl_mp = [r for r in dl if "MP" in str(r.get("Sheet", ""))]
    rows = []
    for cc, d in zip(cand, dl_mp):
        if str(d.get("Particle_Conc_Basis", "")) != "measured":
            continue
        if cc["_invivo"] != "True":
            continue
        if cc["_cdir"] == "<":
            continue
        sz = parse_size(cc["Size_um"])
        if sz is None:
            continue
        ev = parse_num(cc["Endpoint_Value"])
        if ev is None:
            continue
        fp = parse_num(d.get("Factor_per_mg"))
        szmid = sz[0] if sz[2] is None else sz[2]
        pc = ev * fp if fp else None
        rows.append(dict(sp=cc["Species_Name"], env=cc["Environment"],
                         tier=cc["Tier"], shape=norm_shape(cc["Shape"]),
                         size=szmid, lo=sz[0], hi=sz[1], ev=ev, pc=pc))
    return rows


def csv_reader(path):
    import csv
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


# --------------------------------------------------- band structure (1:1)
def cell_band_pairs(rows, currency):
    """cell -> [(ln size_mid, ln band geomean)] for every cell with >=1 band.
    Cells: species x environment x tier x shape. Bands keyed by size mid
    (track3 convention). currency: 'count'|'mass'."""
    cells = defaultdict(list)
    for r in rows:
        cells[(r["sp"], r["env"], r["tier"], r["shape"])].append(r)
    out = {}
    for k, v in cells.items():
        band = defaultdict(list)
        for r in v:
            val = r["pc"] if currency == "count" else r["ev"]
            if val:
                band[round(r["size"], 3)].append(val)
        pairs = []
        for s, vals in sorted(band.items()):
            g = geomean(vals)
            if g:
                pairs.append((math.log(s), math.log(g)))
        if pairs:
            out[k] = pairs
    return out


def pooled_slope(all_pairs):
    """Species-fixed-effect pooling of within-cell centered rows (track3 1:1).
    SE uses df = n - k with k = number of cells in the dict (incl. 1-band)."""
    dx, dy = [], []
    for pairs in all_pairs.values():
        if len(pairs) < 2:
            continue
        xs = np.array([p[0] for p in pairs])
        ys = np.array([p[1] for p in pairs])
        dx.extend(xs - xs.mean())
        dy.extend(ys - ys.mean())
    n, k = len(dx), len(all_pairs)
    dx = np.array(dx)
    dy = np.array(dy)
    b = float(np.dot(dx, dy) / np.dot(dx, dx))
    resid = dy - b * dx
    se = (math.sqrt(np.sum(resid ** 2) / (n - k) / np.dot(dx, dx))
          if n > k else None)
    return b, se, n


def band_cf(lo, hi, alpha):
    """CF per pipeline convention; open upper bound -> narrow +-0.5 um proxy
    (probe2-verified pair convention)."""
    if hi is None:
        return get_powerlaw_cf(lo - 0.5, lo + 0.5, alpha)
    return get_powerlaw_cf(lo, hi, alpha)


def media_alpha(env):
    return ALPHA_SED if env == "Sediment" else ALPHA_SW


# ---------------------------------------------------------------- panels
def add_label(ax, text):
    ax.text(-0.07, 1.06, text, transform=ax.transAxes, fontsize=10.5,
            fontweight="bold", va="bottom", ha="left")


def scatter_shape(ax, xs, ys, shapes, color, size=26, alpha=0.85,
                  open_mk=False):
    for shape, mk in SHAPES_MK:
        sel = [i for i, s in enumerate(shapes) if s == shape]
        if sel:
            ax.scatter([xs[i] for i in sel], [ys[i] for i in sel],
                       marker=mk, s=size, alpha=alpha,
                       facecolors="none" if open_mk else color,
                       edgecolor=color if open_mk else "white",
                       linewidth=0.9 if open_mk else 0.4, zorder=3)


def cell_links(rows):
    """Within-cell link line segments for panel (a): connect band points
    belonging to the same species x environment x tier x shape cell (same
    test method, only the size band differs). Returns dict cell -> sorted
    [(size, count_level, ln_mass)] with >= 2 bands (v6 fix: x is the band
    size, not the particle-count level)."""
    cells = defaultdict(set)
    pts_by_cell = defaultdict(dict)
    for p in rows:
        k = (p["sp"], p["env"], p["tier"], p["shape"])
        # collect by size-mid (a band may aggregate several rows; store both
        # currencies if present)
        rec = pts_by_cell[k].setdefault(round(p["size"], 3), {})
        rec["size"] = p["size"]
        if p["ln_pc"] is not None:
            rec["ln_pc"] = p["ln_pc"]
        if p.get("ln_mass") is not None:
            rec["ln_mass"] = p["ln_mass"]
        cells[k].add(round(p["size"], 3))
    out = {}
    for k, sizes in cells.items():
        if len(sizes) < 2:
            continue
        recs = pts_by_cell[k]
        seq = []
        for s in sorted(sizes):
            rec = recs.get(s, {})
            if rec.get("ln_pc") is None:
                seq = []
                break
            seq.append((rec["size"], math.exp(rec["ln_pc"]),
                        rec.get("ln_mass")))
        if len(seq) >= 2:
            out[k] = seq
    return out


SPECIES_CYCLE = [
    "#1b9e77", "#7570b3", "#d95f02", "#66a61e", "#e7298a",
    "#a6761d", "#666666", "#56b4e9", "#8c564b", "#c49c94",
]


def panel_a(ax, ax2, pts):
    cnt = [p for p in pts if p["ln_pc"] is not None]
    xs = [p["size"] for p in cnt]
    ys = [math.exp(p["ln_pc"]) for p in cnt]
    shapes = [p["shape"] for p in cnt]
    # v6: within-cell link lines (muted, under markers), colored per species
    sp_list = sorted({p["sp"] for p in pts})
    sp_color = {sp: SPECIES_CYCLE[i % len(SPECIES_CYCLE)]
                for i, sp in enumerate(sp_list)}
    links = cell_links(pts)
    for k in sorted(links):
        seq = links[k]
        col = sp_color.get(k[0], GRAY)
        xs_l = [s[0] for s in seq]
        ys_l = [s[1] for s in seq]
        ax.plot(xs_l, ys_l,
                color=col, linewidth=0.55, alpha=0.30, zorder=1.5,
                solid_capstyle="round")
    scatter_shape(ax, xs, ys, shapes, C_CNT)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Size-band representative diameter ($\\mu$m)")
    ax.set_ylabel("Count endpoint (particles/L)", color=C_CNT)
    ax.tick_params(axis="y", labelcolor=C_CNT)

    xs_m = [p["size"] for p in pts]
    ys_m = [math.exp(p["ln_mass"]) for p in pts]
    shapes_m = [p["shape"] for p in pts]
    scatter_shape(ax2, xs_m, ys_m, shapes_m, C_MASS, size=18, alpha=0.55)
    ax2.set_yscale("log")
    ax2.set_ylabel("Mass endpoint (mg/L)", color=C_MASS)
    ax2.tick_params(axis="y", labelcolor=C_MASS)

    handles = [
        Line2D([0], [0], color=C_CNT, marker="o", linestyle="none",
               markersize=6, label="count $\\cdot$ Fragment (left)"),
        Line2D([0], [0], color=C_CNT, marker="^", linestyle="none",
               markersize=6, label="count $\\cdot$ Fiber (left)"),
        Line2D([0], [0], color=C_MASS, marker="o", linestyle="none",
               markersize=5, alpha=0.7,
               label="mass $\\cdot$ Fragment (right)"),
        Line2D([0], [0], color=C_MASS, marker="^", linestyle="none",
               markersize=5, alpha=0.7,
               label="mass $\\cdot$ Fiber (right)"),
        Line2D([0], [0], color=GRAY, linewidth=0.8, alpha=0.5,
               label="same species $\\times$ test, band series"),
    ]
    # v5: legend INSIDE the empty 400-2000 um column (no marker collision;
    # was upper-right squeeze overlapping the 3000 um fiber cluster)
    ax.legend(handles=handles, loc="center left", frameon=True,
              framealpha=0.9, bbox_to_anchor=(0.52, 0.62), fontsize=6.2,
              handlelength=1.4, labelspacing=0.5)


def panel_b(ax, pts, cc_f, cc_n):
    """v7 (b) STORY PANEL - pooled observed within-cell size slope vs the
    slope assumed by a single-alpha standardization. One sentence: no
    observed stratum reproduces the assumed slope -(alpha-1) = -1.07, so a
    single alpha = 2.07 standardization is not supported by the data (the
    empirical basis for the band-dependence shown in (c) and the <56 um
    extrapolation limit). Fragment split by the 58 um band edge (the pooled
    bands straddle the Kooi validity lower bound): small bands are flat,
    large bands slightly positive, fiber steeper but still not -1.07."""
    # recompute pooled slopes by stratum on the same cell-centered pool
    def pooled_stratum(pairs_dict, cond):
        dx, dy, k = [], [], 0
        for pairs in pairs_dict.values():
            if len(pairs) < 2:
                continue
            xs = np.array([p[0] for p in pairs])
            ys = np.array([p[1] for p in pairs])
            mx = xs.mean()
            sel = [(x - mx, y - ys.mean()) for x, y in zip(xs, ys)
                   if cond(math.exp(x))]
            if not sel:
                continue
            ddx, ddy = zip(*sel)
            dx.extend(ddx)
            dy.extend(ddy)
            k += 1
        dx = np.array(dx)
        dy = np.array(dy)
        b = float(dx @ dy / (dx @ dx))
        resid = dy - b * dx
        n, df = len(dx), len(dx) - k
        se = float(math.sqrt(resid @ resid / df / (dx @ dx)))
        half = float(se * stats.t.ppf(0.975, df))
        return b, half, n

    s_f_small = pooled_stratum(cc_f, lambda m: m <= 58.0)
    s_f_large = pooled_stratum(cc_f, lambda m: m > 58.0)
    s_f_all = pooled_stratum(cc_f, lambda m: True)
    s_n_all = pooled_stratum(cc_n, lambda m: True)
    breq = -(ALPHA_SW - 1.0)  # -1.07: slope a single-alpha standardization assumes

    rows = [
        # (label, b, ci_half, color)
        ("Fragment bands $\\leq$58 $\\mu$m (n=%d)" % s_f_small[2],
         s_f_small[0], s_f_small[1], C_FRAG),
        ("Fragment bands $>$58 $\\mu$m (n=%d)" % s_f_large[2],
         s_f_large[0], s_f_large[1], C_FRAG),
        ("Fragment all bands (n=%d)" % s_f_all[2],
         s_f_all[0], s_f_all[1], C_FRAG),
        ("Fiber bands (n=%d)" % s_n_all[2], s_n_all[0], s_n_all[1], C_FIB),
    ]
    ylab = [r[0] for r in rows][::-1]
    ypos = list(range(len(rows)))[::-1]
    # assume line + flat line FIRST so they sit under the CI bars
    ax.axvline(breq, color="black", linestyle="--", linewidth=1.5, zorder=2,
               label=("slope assumed by single-$\\alpha$ standardization "
                      "($-(\\alpha-1)=-1.07$)"))
    ax.axvline(0, color=GRAY, linestyle=":", linewidth=1.3, zorder=2,
               label="flat trend (no size dependence)")
    for (lab, b, half, col), y in zip(rows, ypos):
        ax.errorbar([b], [y], xerr=[[half], [half]], fmt="o", color=col,
                    markersize=5.5, capsize=4, linewidth=1.6, zorder=5,
                    label=lab + ("" if lab is None else ""))
    ax.set_yticks(ypos)
    ax.set_yticklabels(ylab)
    ax.set_xlabel("Observed within-cell size slope, $\\ln$(count) per "
                  "$\\ln$(size)")
    ax.set_xlim(-7.0, 3.2)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    handles = [
        Line2D([0], [0], color=C_FRAG, marker="o", linestyle="none",
               markersize=5.5, label="Fragment (pooled, 95% CI)"),
        Line2D([0], [0], color=C_FIB, marker="o", linestyle="none",
               markersize=5.5, label="Fiber (pooled, 95% CI)"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=1.5,
               label="assumed slope $-1.07$ ($\\alpha$ = 2.07)"),
        Line2D([0], [0], color=GRAY, linestyle=":", linewidth=1.3,
               label="flat trend (zero)"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=True,
              framealpha=0.92, fontsize=6.0)


def cf_safe(lo, hi, a):
    if abs(a - 1.0) < 0.004:
        return float("nan")
    try:
        return get_powerlaw_cf(lo, hi, a)
    except Exception:
        return float("nan")


def panel_c(ax):
    """v6 (c): CF(alpha) BY SOURCE BAND - the user-agreed message: the alpha
    update changes CF by a different multiple per source band, and source
    bands below the fitted validity lower bound 56 um are extrapolation
    (standardization-limited). Four band curves: 333-5000 (adopted, black),
    56-333 (blue), 20-56 / 5-20 (dashed, below the fitted bound). Markers at
    alpha = 2.07 (crimson, CI whisker) and 1.6 (gray, Koelmans +-10% whisker)
    on the adopted band only; x13 update arrow on the adopted band."""
    grid = np.linspace(1.2, 3.0, 240)
    bands = [
        # (lo, hi, color, linestyle, width, zorder, label)
        (333.0, 5000.0, "black", "-", 2.0, 5,
         "333-5000 $\\mu$m source band (adopted)"),
        (56.0, 333.0, "#2166ac", "-", 1.6, 4,
         "56-333 $\\mu$m source band"),
        (20.0, 56.0, "#1b7837", "--", 1.4, 3.5,
         "20-56 $\\mu$m (below fitted lower bound)"),
        (5.0, 20.0, "#7f2704", "--", 1.4, 3.5,
         "5-20 $\\mu$m (below fitted lower bound)"),
    ]
    for lo, hi, col, ls, lw, z, lab in bands:
        cf = [cf_safe(lo, hi, a) for a in grid]
        ax.plot(grid, cf, color=col, linestyle=ls, linewidth=lw,
                label=lab, zorder=z)

    ax.axvline(ALPHA_DEF, color=GRAY, linestyle=":", linewidth=1.0, zorder=2)
    ax.axvline(ALPHA_SW, color=C_FRAG, linestyle="-", linewidth=1.4, zorder=3)

    # markers on the adopted band: adopted 2.07 (+-0.03) and old 1.6 (+-10%)
    cf16 = get_powerlaw_cf(333.0, 5000.0, ALPHA_DEF)
    cf_dn = get_powerlaw_cf(333.0, 5000.0, ALPHA_DEF * 0.9)
    cf_up = get_powerlaw_cf(333.0, 5000.0, ALPHA_DEF * 1.1)
    ax.errorbar([ALPHA_DEF], [cf16], yerr=[[cf16 - cf_dn], [cf_up - cf16]],
                color=GRAY, fmt="o", markersize=4.5, capsize=3.5,
                linewidth=1.3, zorder=6)
    cf207 = get_powerlaw_cf(333.0, 5000.0, ALPHA_SW)
    cf_lo = get_powerlaw_cf(333.0, 5000.0, ALPHA_SW - 0.03)
    cf_hi = get_powerlaw_cf(333.0, 5000.0, ALPHA_SW + 0.03)
    ax.errorbar([ALPHA_SW], [cf207], yerr=[[cf207 - cf_lo], [cf_hi - cf207]],
                color=C_FRAG, fmt="o", markersize=5.5, capsize=3.5,
                linewidth=1.6, zorder=7)

    # per-band CF-change multiples at the alpha update (computed, asserted)
    mults = {}
    for lo, hi, *_ in bands:
        m = get_powerlaw_cf(lo, hi, ALPHA_SW) / get_powerlaw_cf(lo, hi,
                                                                ALPHA_DEF)
        mults[(lo, hi)] = m
    assert abs(mults[(333.0, 5000.0)] - 13.108) < 0.02, mults
    assert abs(mults[(56.0, 333.0)] - 5.2) < 0.3, mults
    assert abs(mults[(20.0, 56.0)] - 2.8) < 0.3, mults
    assert abs(mults[(5.0, 20.0)] - 1.6) < 0.3, mults

    # x13 arrow on the adopted band only (visual); wording in the legend
    ax.annotate("", xy=(ALPHA_SW, 800), xytext=(ALPHA_DEF, 120),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.1,
                                connectionstyle="arc3,rad=0.15"))
    ax.text(2.24, 150, "adopted\n(Kooi 2021,\n2.07 $\\pm$ 0.03)",
            fontsize=6.2, color=C_FRAG, ha="left", va="bottom", zorder=7)

    ax.set_xlim(1.1, 3.05)
    ax.set_ylim(4, 4e5)
    ax.set_yscale("log")
    ax.set_xlabel("Exponent $\\alpha$")
    ax.set_ylabel("Conversion factor to 1-5000 $\\mu$m band\n"
                  "(by source band)")
    handles = [
        Line2D([0], [0], color="black", linewidth=2.0,
               label="333-5000 $\\mu$m (adopted): CF $\\times$13.1"),
        Line2D([0], [0], color="#2166ac", linewidth=1.6,
               label="56-333 $\\mu$m: CF $\\times$5.2"),
        Line2D([0], [0], color="#1b7837", linestyle="--", linewidth=1.4,
               label="20-56 $\\mu$m: CF $\\times$2.8"),
        Line2D([0], [0], color="#7f2704", linestyle="--", linewidth=1.4,
               label="5-20 $\\mu$m: CF $\\times$1.6"),
        Line2D([0], [0], color=C_FRAG, marker="o", linestyle="none",
               markersize=5.5,
               label="adopted $\\alpha$ = 2.07 (CI $\\times$1.4, Kooi 2021)"),
        Line2D([0], [0], color=GRAY, marker="o", linestyle="none",
               markersize=4.5,
               label="old default 1.6 ($\\pm$10% $\\rightarrow$ $\\times$~2, "
                     "Koelmans)"),
        Line2D([0], [0], color="black", linewidth=1.1, linestyle="-",
               label="$\\alpha$ update 1.6 $\\rightarrow$ 2.07: CF "
                     "$\\times$13.1"),
    ]
    # legend INSIDE axes, lower-right empty block (below the rising curve,
    # right of both reference lines): no curve/line/marker is covered
    ax.legend(handles=handles, loc="lower right", frameon=True,
              framealpha=0.92, fontsize=6.0, handlelength=1.2,
              borderaxespad=1.2, labelspacing=0.45)


def panel_d(ax, d_rows):
    xr = [t["x"] for t in d_rows]
    ym = [t["y"] for t in d_rows]
    shapes = [t["shape"] for t in d_rows]
    lo = min(xr + ym) - 1.0
    hi = max(xr + ym) + 1.0
    ax.plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1.0,
            zorder=2, label="1:1 (requirement met)")
    med = float(np.median(np.array(ym) - np.array(xr)))
    ax.plot([lo, hi], [lo + med, hi + med], color=GRAY, linestyle=":",
            linewidth=1.2, zorder=2, label=f"median offset ({med:+.2f} log)")
    scatter_shape(ax, xr, ym, shapes, C_PAIR, size=30, alpha=0.85,
                  open_mk=True)

    x_rug = hi - 0.45
    ax.plot(np.full(len(d_rows), x_rug),
            [y - x for x, y in zip(xr, ym)],
            "|", color=GRAY, alpha=0.6, markersize=6, zorder=3)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Required pair ratio $\\ln$(CF$_L$/CF$_S$)")
    ax.set_ylabel("Measured pair ratio $\\ln$(PC$_S$/PC$_L$)")
    handles = [
        Line2D([0], [0], color=C_PAIR, marker="o", linestyle="none",
               markersize=6, markerfacecolor="none", label="pair $\\cdot$ Fragment"),
        Line2D([0], [0], color=C_PAIR, marker="^", linestyle="none",
               markersize=6, markerfacecolor="none", label="pair $\\cdot$ Fiber (n=1)"),
        Line2D([0], [0], color=GRAY, linestyle=":", linewidth=1.2,
               label="residuals (right margin)"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=6.2)


# ----------------------------------------------------------------- main
def main():
    rows = build_pool()
    assert len(rows) == 490, f"pool rows = {len(rows)} (expect 490)"
    sw = [r for r in rows if r["env"] == "Seawater"]
    assert len(sw) == 362, f"seawater rows = {len(sw)} (expect 362)"

    frag = [r for r in sw if r["shape"] == "Fragment"]
    fib = [r for r in sw if r["shape"] == "Fiber"]
    cc_f = cell_band_pairs(frag, "count")
    cc_n = cell_band_pairs(fib, "count")
    cm_f = cell_band_pairs(frag, "mass")
    cm_n = cell_band_pairs(fib, "mass")

    multi_f = [k for k, p in cc_f.items() if len(p) >= 2]
    multi_n = [k for k, p in cc_n.items() if len(p) >= 2]
    valid_f = [k for k, p in cc_f.items() if len(p) >= 3]
    valid_n = [k for k, p in cc_n.items() if len(p) >= 3]
    assert len(cc_f) == 18 and len(multi_f) == 16 and len(valid_f) == 16, \
        (len(cc_f), len(multi_f), len(valid_f))
    assert len(cc_n) == 15 and len(multi_n) == 14 and len(valid_n) == 1, \
        (len(cc_n), len(multi_n), len(valid_n))

    b_f, se_f, n_f = pooled_slope(cc_f)
    b_n, se_n, n_n = pooled_slope(cc_n)
    bm_f, sem_f, nm_f = pooled_slope(cm_f)
    bm_n, sem_n, nm_n = pooled_slope(cm_n)
    assert n_f == 48 and n_n == 29, (n_f, n_n)
    assert abs(b_f - 0.1074) < 0.02 and abs(se_f - 0.3427) < 0.02, (b_f, se_f)
    assert abs(b_n - (-1.3293)) < 0.02 and abs(se_n - 0.2336) < 0.02, \
        (b_n, se_n)
    assert abs(bm_f - 0.2879) < 0.02 and abs(sem_f - 0.2589) < 0.02, \
        (bm_f, sem_f)
    assert abs(bm_n - 0.2802) < 0.02 and abs(sem_n - 0.2193) < 0.02, \
        (bm_n, sem_n)

    # panel points: bands of the multi-band cells (count + mass), both shapes
    band_range = {}
    for r in sw:
        band_range.setdefault((r["shape"], r["sp"], r["tier"],
                               round(r["size"], 3)), (r["lo"], r["hi"]))
    pts = []
    for shape, ccnt, cmas in (("Fragment", cc_f, cm_f), ("Fiber", cc_n, cm_n)):
        for k, bands in ccnt.items():
            if len(bands) < 2 or k not in cmas:
                continue
            for (ln_s, ln_pc), (_, ln_mass) in zip(bands, cmas[k]):
                pts.append(dict(sp=k[0], env=k[1], tier=k[2], shape=shape,
                                size=math.exp(ln_s), ln_pc=ln_pc,
                                ln_mass=ln_mass, dev=None, ln_s=ln_s))
    assert len(pts) == 77, f"panel points = {len(pts)} (expect 77)"

    # (b) standardize at alpha = 2.07, within-cell centered
    cells_all = defaultdict(list)
    for p in pts:
        cells_all[(p["sp"], p["env"], p["tier"], p["shape"])].append(p)
    for k, vals in cells_all.items():
        std = []
        for p in vals:
            lo, hi = band_range[(k[3], k[0], k[2], round(p["size"], 3))]
            cf = band_cf(lo, hi, ALPHA_SW)
            std.append(p["ln_pc"] + math.log(cf))
        my = sum(std) / len(std)
        for p, s in zip(vals, std):
            p["dev"] = s - my

    # (d) direct pairs: smallest vs largest band per regression-valid cell
    cells_raw = defaultdict(list)
    for r in sw:
        cells_raw[(r["sp"], r["env"], r["tier"], r["shape"])].append(r)
    d_rows, excluded = [], []
    for (sp, env, tier, shape), v in sorted(cells_raw.items()):
        if (sp, env, tier, shape) not in valid_f + valid_n:
            if len({r["size"] for r in v}) >= 2:
                excluded.append(shape)
            continue
        bysize = defaultdict(list)
        for r in v:
            if r["pc"]:
                bysize[(r["size"], r["lo"], r["hi"])].append(r["pc"])
        bands = sorted(bysize, key=lambda t: t[0])
        slo, shi = bands[0], bands[-1]
        pcS = geomean(bysize[slo])
        pcL = geomean(bysize[shi])
        a = media_alpha(env)
        cfS = band_cf(slo[1], slo[2], a)
        cfL = band_cf(shi[1], shi[2], a)
        d_rows.append(dict(sp=sp, env=env, tier=tier, shape=shape,
                           s=slo[0], l=shi[0],
                           x=math.log(cfL / cfS), y=math.log(pcS / pcL)))
    assert len(d_rows) == 17, f"(d) pairs = {len(d_rows)} (expect 17)"
    n_ex_f = sum(1 for s in excluded if s == "Fragment")
    n_ex_n = sum(1 for s in excluded if s == "Fiber")

    cf16 = get_powerlaw_cf(333.0, 5000.0, ALPHA_DEF)
    up = get_powerlaw_cf(333.0, 5000.0, ALPHA_DEF * 1.1) / cf16
    dn = cf16 / get_powerlaw_cf(333.0, 5000.0, ALPHA_DEF * 0.9)
    assert abs(up - 2.342) < 0.02 and abs(dn - 2.235) < 0.02, (up, dn)

    a_imp_f, a_imp_f_sd = 1.0 - b_f, se_f
    a_imp_n, a_imp_n_sd = 1.0 - b_n, se_n

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.6))
    fig.subplots_adjust(hspace=0.34, wspace=0.30, top=0.95, bottom=0.15,
                        left=0.085, right=0.925)
    panel_a(axes[0][0], axes[0][0].twinx(), pts)
    panel_b(axes[0][1], pts, cc_f, cc_n)
    panel_c(axes[1][0])
    panel_d(axes[1][1], d_rows)

    add_label(axes[0][0], "(a) Raw seawater endpoints (before standardization)")
    add_label(axes[0][1], "(b) Observed size slope vs single-$\\alpha$ assumption")
    add_label(axes[1][0], "(c) Adopted-$\\alpha$ uncertainty (333-5000 $\\mu$m)")
    add_label(axes[1][1], "(d) Required vs measured (seawater direct pairs)")

    out_png = FIG_DIR / "track1_fourpanel_applicability.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    req = -(ALPHA_SW - 1.0)
    gap_f = abs(req - b_f) / se_f
    gap_n = abs(req - b_n) / se_n
    print("== Track 1 four-panel applicability (v7, seawater-only) ==")
    print(f"pool rows: {len(rows)} (SW {len(sw)}: frag {len(frag)} + "
          f"fiber {len(fib)})")
    links_a = cell_links(pts)
    print(f"v5 link lines (same species x test, >=2 bands): {len(links_a)}")
    n_link_rows = sum(len(v) for v in links_a.values())
    print(f"v5 linked band points: {n_link_rows}")
    print(f"fragment: cells {len(cc_f)}, multi {len(multi_f)}, "
          f"valid {len(valid_f)}, bands {n_f}")
    print(f"fiber   : cells {len(cc_n)}, multi {len(multi_n)}, "
          f"valid {len(valid_n)}, bands {n_n}")
    print(f"count: beta frag = {b_f:+.4f} +/- {se_f:.4f} (required {req:+.2f}; "
          f"gap {gap_f:.2f} sigma) | beta fiber = {b_n:+.4f} +/- {se_n:.4f} "
          f"(gap {gap_n:.2f} sigma)")
    print(f"mass : beta frag = {bm_f:+.4f} +/- {sem_f:.4f} | "
          f"beta fiber = {bm_n:+.4f} +/- {sem_n:.4f} (required: 0.00)")
    print(f"implied alpha: fragment {a_imp_f:.2f} +/- {a_imp_f_sd:.2f} | "
          f"fiber {a_imp_n:.2f} +/- {a_imp_n_sd:.2f} (Kooi marine SW 2.07)")
    print(f"Koelmans bracket (333-5000, alpha 1.6 +/-10%): "
          f"x{up:.2f} up / x{dn:.2f} down")
    med = float(np.median([t["y"] - t["x"] for t in d_rows]))
    print(f"(d) pairs: {len(d_rows)} ({sum(1 for t in d_rows if t['shape'] == 'Fragment')} frag + "
          f"{sum(1 for t in d_rows if t['shape'] == 'Fiber')} fiber); "
          f"median offset (meas - req) = {med:+.2f} log")
    print(f"(d) excluded (<3 bands): {n_ex_f} fragment, {n_ex_n} fiber cells")
    for t in d_rows[:3]:
        print(f"    e.g. {t['sp'][:26]:26s} {t['shape'][:9]:9s} "
              f"{t['tier'][:6]:6s} S={t['s']:g} L={t['l']:g} "
              f"lnReq={t['x']:+.2f} lnMeas={t['y']:+.2f}")
    for t in d_rows:
        if t["shape"] == "Fiber":
            print(f"    fiber {t['sp'][:26]:26s} {t['tier'][:6]:6s} "
                  f"S={t['s']:g} L={t['l']:g} "
                  f"lnReq={t['x']:+.2f} lnMeas={t['y']:+.2f}")
    print(f"== saved: {out_png} ==")


if __name__ == "__main__":
    main()