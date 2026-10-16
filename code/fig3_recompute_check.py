#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig3 render-input verification (2026-10-13).

Recompute the engine HC5/HC10 + parametric-bootstrap CIs directly from the
ssdtools input CSVs (fig3_merge/ssdtools_inputs/) and compare with
TR_abcd_merged_T12_20_300.csv (the user-confirmed a,b,c,d set).
Engine formulas (ssd_tomex.py, spec alpha=2.07, window 20-300 um):
  fit_lognormal_hc: mu=mean(log v), sd=std(log v, ddof=0), HC_p=exp(mu+z_p*sd)
  bootstrap_hc    : parametric, seed=20261011, nboot=1000,
                    x=exp(N(mu,sigma,(nboot,len))), lcl/ucl = 2.5%/97.5% quantiles
  TR mapping      : TR1=HC5 lcl (q1), TR2=HC5 (q1), TR3=HC5 (median), TR4=HC10 (median)
"""
import numpy as np
import pandas as pd

from pathlib import Path
IN = str(Path(__file__).resolve().parents[1] / "data")

FILES = {
    ("food_dilution", "q1_all"): "ssdtools_input_marine_food_dilution_q1_all.csv",
    ("food_dilution", "median_organism_population"): "ssdtools_input_marine_food_dilution_median_organism_population.csv",
    ("translocation", "q1_all"): "ssdtools_input_marine_translocation_q1_all.csv",
    ("translocation", "median_organism_population"): "ssdtools_input_marine_translocation_median_organism_population.csv",
}


def engine_hc(values):
    v = np.array(values, dtype=float)
    lx = np.log(v)
    mu = float(lx.mean())
    sigma = float(lx.std(ddof=0))
    hc5 = float(np.exp(mu + norm_ppf(0.05) * sigma))
    hc10 = float(np.exp(mu + norm_ppf(0.10) * sigma))
    # parametric bootstrap, exact engine recipe
    nboot, seed = 1000, 20261011
    rng = np.random.default_rng(seed)
    x = np.exp(rng.normal(mu, sigma, size=(nboot, len(v))))
    mus = np.log(x).mean(axis=1)
    sigs = np.log(x).std(axis=1, ddof=0)
    b5 = np.exp(mus + norm_ppf(0.05) * sigs)
    b10 = np.exp(mus + norm_ppf(0.10) * sigs)
    return {
        "n": len(v), "mu": mu, "sigma": sigma,
        "HC5": hc5, "HC10": hc10,
        "HC5_lcl": float(np.quantile(b5, 0.025)),
        "HC5_ucl": float(np.quantile(b5, 0.975)),
        "HC10_lcl": float(np.quantile(b10, 0.025)),
        "HC10_ucl": float(np.quantile(b10, 0.975)),
    }


def norm_ppf(p):
    from scipy.stats import norm
    return float(norm.ppf(p))


def main():
    tr = pd.read_csv(f"{IN}/TR_abcd_merged_T12_20_300.csv", encoding="utf-8-sig")
    calc = {}
    for (mech, collapse), fn in FILES.items():
        df = pd.read_csv(f"{IN}/{fn}", encoding="utf-8-sig")
        calc[(mech, collapse)] = {
            "vals": df["Conc"].to_numpy(float),
            "n_eff": int(df["n_effects"].sum()),
            "hc": engine_hc(df["Conc"].to_numpy(float)),
        }

    def close(a, b, tol=1e-6):
        if a == b:
            return True
        return abs(a - b) / max(abs(b), 1e-300) <= tol

    print("=== cohort summary (from input CSVs) ===")
    for k in FILES:
        hc = calc[k]["hc"]
        print(f"{k[0]:>13s} {k[1]:>28s}: n_species={hc['n']} n_effects={calc[k]['n_eff']} "
              f"HC5={hc['HC5']:.6g} HC10={hc['HC10']:.6g} lcl5={hc['HC5_lcl']:.6g}")

    print()
    print("=== TR verification vs TR_abcd_merged_T12_20_300.csv ===")
    rows = []
    nfail = 0
    for _, r in tr.iterrows():
        mech = str(r["Mechanism"])
        collapse = str(r["SSD_Collapse"])
        trn = str(r["TR"])
        val = float(r["particles_L"])
        hc = calc[(mech, collapse)]["hc"]
        if trn == "TR1":
            tgt = hc["HC5_lcl"]
        elif trn == "TR2":
            tgt = hc["HC5"]
        elif trn == "TR3":
            tgt = hc["HC5"]
        elif trn == "TR4":
            tgt = hc["HC10"]
        else:
            continue
        ok = close(tgt, val)
        nfail += 0 if ok else 1
        rows.append((mech, trn, collapse, val, tgt, ok))
        print(f"{mech:>13s} {trn} [{collapse:>28s}] csv={val:>16.10g} recomputed={tgt:>16.10g} "
              f"{'MATCH' if ok else 'MISMATCH'}")
    print()
    print(f"VERDICT: {'ALL_MATCH' if nfail == 0 else f'{nfail} MISMATCH(ES)'}")


if __name__ == "__main__":
    main()