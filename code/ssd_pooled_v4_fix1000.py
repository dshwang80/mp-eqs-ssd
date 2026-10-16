# -*- coding: utf-8 -*-
# ssd_pooled_v4_fix1000.py (v2) — x1000 defect correction for the conventional-approach pool ([B])
# Defect (recorded line, ssd_pooled_v3_extract.py: mass_g = vol*dens*1e3, actual 1e6)
#   -> particle count x1000 inflated in the ToMEx-side components of the v4 pool.
# Affected: 16 ToMEx-containing cells of locked ssd_pooled_reps_v4.csv (14 pure + 2 merged).
# Fix: recompute ONLY those 16 cells with ToMEx components /1000; the 21 non-ToMEx cells
#      (ours 17 + Jung2021 4) are carried over VERBATIM from the locked v4 pool
#      (incl. Cyprinodon variegatus Chronic, which no longer exists in current source files):
#      minimal intervention, locked originals untouched; outputs are new files.
# Engine: hc5_modelavg canonical identical to ssd_pooled_v4_merge.py (lnorm MLE + llogis).
import math
import numpy as np
import pandas as pd
import scipy.stats as st

from pathlib import Path
MP = str(Path(__file__).resolve().parents[1])
OUT = f"{MP}/output/py_out"
LOG = []
def P(*a):
    print(*a); LOG.append(" ".join(str(x) for x in a))

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

# --- locked v4 pool
lock = pd.read_csv(f"{OUT}/ssd_pooled_reps_v4.csv", encoding="utf-8-sig")
lock["Species"] = lock["Species"].astype(str).str.strip()
lockmap = {(r["Tier"], r["Species"]): r for _, r in lock.iterrows()}
P(f"locked v4 pool: {len(lock)} cells "
  f"(ToMEx-containing {sum('ToMEx' in str(s) for s in lock['src'])})")

# --- defective ToMEx side from the v3 extract chain (same gates as v4 merge)
vx = pd.read_csv(f"{OUT}/v3_extract_reps.csv", encoding="utf-8-sig")
raw_sum = float(vx["gm_p"].sum())
vxa = vx[~vx["Species"].isin(ALGAE)].copy()
vxa["Species"] = vxa["Species"].str.strip()
v3h = vxa[vxa["Metric"] != "HONEC"].copy()
v3map = {}
for (tier, sp), g in v3h.groupby(["Study_Tier", "Species"]):
    v3map[(tier, sp)] = dict(gm=gm(g["gm_p"]), mn=float(g["min_p"].min()),
                             metrics="|".join(sorted(g["Metric"])), reps=int(len(g)))

# --- ours side for the 2 merged cells (drift-verified: recompute == locked within 1e-9)
v2 = pd.read_csv(f"{OUT}/ssd_pooled_reps_v2.csv", encoding="utf-8-sig")
sp2 = pd.read_csv(f"{OUT}/species_rep_v2.csv", encoding="utf-8-sig")
sp2["학명"] = sp2["학명"].astype(str).str.strip()
TMAP = {"급성": "Acute", "만성": "Chronic"}
mn_lookup = {(TMAP[str(t)], s): float(m) for t, s, m in zip(sp2["급만성"], sp2["학명"], sp2["min_particle"])}
v2map = {}
for tier, g in v2.groupby("Tier"):
    for sp, gg in g.groupby("Species"):
        spc = str(sp).strip()
        gm_r = gg[gg["Agg"] == "GM"]; mn_r = gg[gg["Agg"] == "MIN"]
        gm_v = float(gm_r["Rep"].iloc[0]) if len(gm_r) else float(gg["Rep"].iloc[0])
        mn_v = float(mn_r["Rep"].iloc[0]) if len(mn_r) else mn_lookup.get((tier, spc), gm_v)
        v2map[(tier, spc)] = dict(gm=gm_v, mn=mn_v, src=str(gg["Source"].iloc[0]))

# --- corrected pool
merged = {}
n_keep = 0
for _, r in lock.iterrows():
    key = (r["Tier"], r["Species"])
    src = str(r["src"])
    if "ToMEx" not in src:
        merged[key] = dict(gm=float(r["gm"]), mn=float(r["mn"]), src=src,
                           metrics=r.get("metrics"), reps=r.get("reps"))
        n_keep += 1
n_fixed, n_merged = 0, 0
for key, info in v3map.items():
    lk = lockmap.get(key)
    if lk is None:
        P(f"[AUDIT-UNEXPECTED] v3map key not in locked pool: {key}")
        merged[key] = dict(gm=info["gm"] / 1000.0, mn=info["mn"] / 1000.0,
                           src="ToMEx_fix1000", metrics=info["metrics"], reps=info["reps"])
        n_fixed += 1
        continue
    if "ours" in str(lk["src"]):        # merged cell: ours untouched, ToMEx component /1000
        o = v2map.get(key)
        if o is None:
            P(f"[AUDIT-ABORT] merged cell missing ours side: {key}")
            raise SystemExit(2)
        merged[key] = dict(gm=gm([o["gm"], info["gm"] / 1000.0]),
                           mn=min(o["mn"], info["mn"] / 1000.0),
                           src=str(lk["src"]) + "+ToMEx_fix1000",
                           metrics=lk.get("metrics"), reps=lk.get("reps"))
        n_merged += 1
    else:                                # pure ToMEx cell
        merged[key] = dict(gm=info["gm"] / 1000.0, mn=info["mn"] / 1000.0,
                           src="ToMEx_fix1000", metrics=info["metrics"], reps=info["reps"])
        n_fixed += 1
P(f"pool: verbatim non-ToMEx {n_keep} + merged-fixed {n_merged} + pure ToMEx-fixed {n_fixed} "
  f"= {len(merged)} cells")
rows = []
for tier in sorted({k[0] for k in merged}):
    for track, fld in (("GM", "gm"), ("MIN", "mn")):
        vals = [d[fld] for k, d in merged.items() if k[0] == tier]
        hc, w_n, (mu, sd), _ = hc5_modelavg(vals)
        rows.append(dict(Tier=tier, Track=track, n_species=len(vals), gm_pool=gm(vals),
                         log_mu=mu, log_sd=sd, w_lnorm=w_n, hc5=hc))
        P(f"{tier}-{track}: n={len(vals)} poolGM={gm(vals):.3f} HC5={hc:.3f} (w_lnorm={w_n:.3f})")

fix = pd.DataFrame(rows)
oldhc = pd.read_csv(f"{OUT}/ssd_pooled_hc5_v4.csv", encoding="utf-8-sig")
P("== HC5 v4(defective) -> v4_fix1000 ==")
for _, r in fix.iterrows():
    m = oldhc[(oldhc["Tier"] == r["Tier"]) & (oldhc["Track"] == r["Track"])]
    o = float(m["hc5"].iloc[0]); n = float(r["hc5"])
    P(f"  {r['Tier']:7s}-{r['Track']:3s} n={int(r['n_species']):2d} HC5 {o:12.1f} -> {n:12.1f} (x{n/o:.4g})")

# --- save (no overwriting of locked v4 files)
fix.to_csv(f"{OUT}/ssd_pooled_hc5_v4_fix1000.csv", index=False, encoding="utf-8-sig")
out_reps = pd.DataFrame([dict(Tier=k[0], Species=k[1], src=d["src"], gm=d["gm"], mn=d["mn"],
                              metrics=d.get("metrics"), reps=d.get("reps"))
                         for k, d in merged.items()])
out_reps = out_reps.sort_values(["Tier", "Species"])
out_reps.to_csv(f"{OUT}/ssd_pooled_reps_v4_fix1000.csv", index=False, encoding="utf-8-sig")
with open(f"{OUT}/ssd_pooled_v4_fix1000_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG))
chk = pd.read_csv(f"{OUT}/ssd_pooled_reps_v4_fix1000.csv", encoding="utf-8-sig")
P(f"[SAVED+READBACK] hc5 rows={len(fix)} reps cells={len(chk)} "
  f"({sum('fix1000' in str(s) for s in chk['src'])} corrected / "
  f"{sum('fix1000' not in str(s) for s in chk['src'])} verbatim) | raw v3 gm_p sum={raw_sum:.6f}")