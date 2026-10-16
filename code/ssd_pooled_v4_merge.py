# -*- coding: utf-8 -*-
# ssd_pooled_v4_merge.py — v4 merge (v3 engine + two fixes, user directive: no acute/chronic inversion)
# FIX-1: exclude ToMEx rows with Metric == "HONEC" (right-censored top-of-design
#        no-effect dose, not an effect threshold; sole cause of Crassostrea
#        chronic 1e9 vs acute 1.27e5 = x7889 inversion). Uniform rule, both
#        tiers, no per-species selection.
# FIX-2: v2 MIN track now uses species_rep_v2.csv min_particle (ours). The v3
#        merge fell back to GM because ssd_pooled_reps_v2.csv has no MIN rows;
#        the v3 header spec ("species_rep_v2.csv; min_particle -> MIN track")
#        was never implemented. Jung2021 species: single NOEC -> mn = gm.
# HC5   : hc5_modelavg canonical (lnorm MLE ddof=0 + llogis, w=exp(ll-2) norm).
import math
import numpy as np
import pandas as pd
import scipy.stats as st

from pathlib import Path
MP = str(Path(__file__).resolve().parents[1])
OUT = f"{MP}/output/py_out"


def P(*a):
    print(*a)


def gm(vals):
    return float(np.exp(np.mean(np.log(np.asarray(vals, float)))))


def hc5_modelavg(vals):
    x = np.log(np.asarray(vals, float))
    mu, sd = x.mean(), x.std(ddof=0)
    ll_n = float(np.sum(st.norm.logpdf(x, mu, sd)))
    a_l, b_l, *_ = st.logistic.fit(x)
    ll_l = float(np.sum(st.logistic.logpdf(x, a_l, b_l)))
    en, el = math.exp(ll_n - 2.0), math.exp(ll_l - 2.0)
    w_n = en / (en + el)
    h_n = float(st.norm.ppf(0.05, mu, sd))
    h_l = float(st.logistic.ppf(0.05, a_l, b_l))
    return math.exp(w_n * h_n + (1.0 - w_n) * h_l), w_n, (mu, sd), (a_l, b_l)


ALGAE = ["Dunaliella salina", "Karenia mikimotoi", "Skeletonema costatum"]

v3 = pd.read_csv(f"{OUT}/v3_extract_reps.csv", encoding="utf-8-sig")
v3a = v3[~v3["Species"].isin(ALGAE)].copy()
v3a["Species"] = v3a["Species"].str.strip()
honec_rows = int((v3a["Metric"] == "HONEC").sum())
v3h = v3a[v3a["Metric"] != "HONEC"].copy()
P(f"ToMEx reps: non-algae {len(v3a)} rows / {v3a['Species'].nunique()} species; "
  f"HONEC removed {honec_rows} rows -> {len(v3h)} rows / {v3h['Species'].nunique()} species")
lost = sorted(set(zip(v3a["Study_Tier"], v3a["Species"])) - set(zip(v3h["Study_Tier"], v3h["Species"])))
P(f"ToMEx (tier,species) lost to HONEC-only: {len(lost)} -> {[f'{t}/{s}' for t, s in lost]}")

v3map = {}
for (tier, sp), g in v3h.groupby(["Study_Tier", "Species"]):
    v3map[(tier, sp)] = dict(gm=gm(g["gm_p"]), mn=float(g["min_p"].min()),
                             metrics="|".join(sorted(g["Metric"])), reps=int(len(g)))

v2 = pd.read_csv(f"{OUT}/ssd_pooled_reps_v2.csv", encoding="utf-8-sig")
sp2 = pd.read_csv(f"{OUT}/species_rep_v2.csv", encoding="utf-8-sig")
sp2["학명"] = sp2["학명"].astype(str).str.strip()
TMAP = {"급성": "Acute", "만성": "Chronic"}
mn_lookup = {(TMAP[str(t)], s): float(m) for t, s, m in zip(sp2["급만성"], sp2["학명"], sp2["min_particle"])}

v2map = {}
for tier, g in v2.groupby("Tier"):
    for sp, gg in g.groupby("Species"):
        spc = str(sp).strip()
        gm_r = gg[gg["Agg"] == "GM"]
        mn_r = gg[gg["Agg"] == "MIN"]
        gm_v = float(gm_r["Rep"].iloc[0]) if len(gm_r) else float(gg["Rep"].iloc[0])
        mn_v = float(mn_r["Rep"].iloc[0]) if len(mn_r) else mn_lookup.get((tier, spc), gm_v)
        v2map[(tier, spc)] = dict(gm=gm_v, mn=mn_v, src=str(gg["Source"].iloc[0]))

n_hit = sum(1 for k in v2map if (k[0], k[1]) in mn_lookup)
P(f"v2 pool: {len(v2)} rows / {v2['Species'].nunique()} species; MIN restored from species_rep_v2: {n_hit}/{len(v2map)}")

merged = {}
for key, info in v2map.items():
    ov = v3map.get(key)
    if ov is None:
        merged[key] = dict(gm=info["gm"], mn=info["mn"], src=info["src"])
    else:
        merged[key] = dict(gm=gm([info["gm"], ov["gm"]]), mn=min(info["mn"], ov["mn"]),
                           src=info["src"] + "+ToMEx")

added = []
for (tier, sp), info in v3map.items():
    if (tier, sp) not in v2map:
        merged[(tier, sp)] = dict(info, src="ToMEx")
        added.append(f"{tier}/{sp}")
P(f"ToMEx-only species added: {len(added)}")

final_rows = []
for tier in sorted({k[0] for k in merged}):
    for track, fld in (("GM", "gm"), ("MIN", "mn")):
        vals = [d[fld] for k, d in merged.items() if k[0] == tier]
        hc, w_n, (mu, sd), _ = hc5_modelavg(vals)
        final_rows.append(dict(Tier=tier, Track=track, n_species=len(vals), gm_pool=gm(vals),
                               log_mu=mu, log_sd=sd, w_lnorm=w_n, hc5=hc))
        P(f"{tier}-{track}: n={len(vals)} poolGM={gm(vals):.1f} HC5={hc:.1f} (w_lnorm={w_n:.3f})")

# ---- acute vs chronic inversion audit (species level) -----------------------
by_tier = {}
for (tier, sp), d in merged.items():
    by_tier.setdefault(tier, {})[sp] = d
shared = sorted(set(by_tier.get("Acute", {})) & set(by_tier.get("Chronic", {})))
P(f"shared species (Acute&Chronic): {len(shared)}")
inv = 0
for sp in shared:
    ra, rc = by_tier["Acute"][sp], by_tier["Chronic"][sp]
    r_gm = rc["gm"] / ra["gm"]
    r_mn = rc["mn"] / ra["mn"]
    flag = " <-- CHRONIC>ACUTE" if (r_gm > 1.0 or r_mn > 1.0) else ""
    if r_gm > 1.0 or r_mn > 1.0:
        inv += 1
    P(f"  {sp:32s} gm {ra['gm']:.3e} -> {rc['gm']:.3e} (x{r_gm:.4g}) | "
      f"mn {ra['mn']:.3e} -> {rc['mn']:.3e} (x{r_mn:.4g}){flag}")
P(f"species-level inversions remaining: {inv}")

pd.DataFrame(final_rows).to_csv(f"{OUT}/ssd_pooled_hc5_v4.csv", index=False, encoding="utf-8-sig")
pd.DataFrame([dict(Tier=k[0], Species=k[1], **d) for k, d in merged.items()]).to_csv(
    f"{OUT}/ssd_pooled_reps_v4.csv", index=False, encoding="utf-8-sig")
P("saved: ssd_pooled_hc5_v4.csv / ssd_pooled_reps_v4.csv")

try:
    v3hc = pd.read_csv(f"{OUT}/ssd_pooled_hc5_v3.csv", encoding="utf-8-sig")
    v4hc = pd.read_csv(f"{OUT}/ssd_pooled_hc5_v4.csv", encoding="utf-8-sig")
    P("== HC5 v3 -> v4 ==")
    for _, r4 in v4hc.iterrows():
        m = v3hc[(v3hc["Tier"] == r4["Tier"]) & (v3hc["Track"] == r4["Track"])]
        old = float(m["hc5"].iloc[0]) if len(m) else float("nan")
        P(f"  {r4['Tier']:7s}-{r4['Track']:3s} n={int(r4['n_species']):2d} HC5 {old:12.1f} -> {float(r4['hc5']):12.1f}")
except Exception as e:
    P(f"(v3 comparison skipped: {e})")