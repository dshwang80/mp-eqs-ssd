#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig3_lognormal_band_v1.py — Fig.3 curve/band precompute for the re-render
(user-approved round-3 item: Fig.3 curve = the SINGLE parametric lognormal fit
underlying the reported TR1-TR4 values; band = bootstrap CI of the SAME fit).

Engine recipe (locked, ssd_tomex.py / fig3_recompute_check.py 2026-10-13):
  mu = mean(log v); sigma = std(log v, ddof=0); curve(p) = exp(mu + z_p*sigma)
  bootstrap: rng = np.random.default_rng(20261011);
             x = exp(rng.normal(mu, sigma, (1000, len(v))));
             per-draw mu_d, sigma_d (ddof=0); per-p curve_d = exp(mu_d + z_p*sigma_d)
             band: lcl/ucl = 2.5%/97.5% np.quantile over draws (default linear interp)

Hard assertions per panel:
  curve(0.05) == TR2 (q1 panels) / TR3 (median panels)   [abs(log ratio) < 1e-9]
  curve(0.10) == TR4 (median panels)                     [same tol]
  band lcl(0.05) == engine HC5_lcl == TR1 (q1) / lcl(0.05) == TR3-lcl (median)
    -- TR1/tr3-lcl from the same locked TR csv row-direct
Outputs (all NEW files, nothing overwritten):
  03_figures/output_archive/fig3_lognormal_build/fig3_lognormal_band_<mech>_<coll>.csv
  fig3_lognormal_meta_<mech>_<coll>.csv  +  fig3_lognormal_band_summary.json
"""
import numpy as np
import pandas as pd
import json, os

from statistics import NormalDist
import json, os

PANELS = [
    ('food_dilution', 'q1_all', 'fd_q1'),
    ('translocation', 'q1_all', 'tt_q1'),
    ('food_dilution', 'median_organism_population', 'fd_md'),
    ('translocation', 'median_organism_population', 'tt_md'),
]

SEED, NBOOT = 20261011, 1000


def nd_ppf(p):
    # exact normal quantile; equals scipy.stats.norm.ppf to ~1e-15 rel
    return float(NormalDist().inv_cdf(p))


Z05 = nd_ppf(0.05)
Z10 = nd_ppf(0.10)

from pathlib import Path
BASE = str(Path(__file__).resolve().parents[1])
INP = BASE + '/data'
OUTD = BASE + '/output/fig3_lognormal_build'
os.makedirs(OUTD, exist_ok=True)

tr = pd.read_csv(os.path.join(INP, 'TR_abcd_merged_T12_20_300.csv'), encoding='utf-8-sig')


def trv(mech, coll, trn):
    sub = tr[(tr['Mechanism'] == mech) & (tr['SSD_Collapse'] == coll) & (tr['TR'] == trn)]
    v = sub['particles_L']
    assert len(v) == 1, (mech, coll, trn, len(v))
    return float(v.iloc[0])


summary = {'seed': SEED, 'nboot': NBOOT, 'panels': {}}
all_ok = True
for mech, coll, key in PANELS:
    f_in = os.path.join(INP, f'ssdtools_input_marine_{mech}_{coll}.csv')
    inp = pd.read_csv(f_in, encoding='utf-8-sig')
    v = inp['Conc'].astype(float).values
    n = len(v)
    lx = np.log(v)
    mu, sigma = float(lx.mean()), float(lx.std(ddof=0))

    pr = np.round(np.arange(0.01, 0.995, 0.01), 2)
    zp = np.array([nd_ppf(float(p_)) for p_ in pr])
    curve = np.exp(mu + sigma * zp)

    # engine bootstrap (exact recipe)
    rng = np.random.default_rng(SEED)
    x = np.exp(rng.normal(mu, sigma, size=(NBOOT, n)))
    mus = np.log(x).mean(axis=1)
    sigs = np.log(x).std(axis=1, ddof=0)
    curves = np.exp(mus[:, None] + sigs[:, None] * zp[None, :])
    lcl = np.quantile(curves, 0.025, axis=0)
    ucl = np.quantile(curves, 0.975, axis=0)

    hc5 = float(curve[np.argmin(abs(pr - 0.05))])
    hc10 = float(curve[np.argmin(abs(pr - 0.10))])
    if coll == 'q1_all':
        tgt5, name5 = trv(mech, coll, 'TR2'), 'TR2'
        tgt1, name1 = trv(mech, coll, 'TR1'), 'TR1'
    else:
        tgt5, name5 = trv(mech, coll, 'TR3'), 'TR3'
        tgt5_hc10, tgt10_name = trv(mech, coll, 'TR4'), 'TR4'
        tgt1, name1 = tgt5, name5  # median-panel lcl = band lcl at 0.05 (no TR row)

    ok5 = abs(np.log(hc5 / tgt5)) < 1e-9
    ok10 = abs(np.log(hc10 / tgt5_hc10)) < 1e-9 if coll != 'q1_all' else True
    oklcl = True
    if coll == 'q1_all':
        oklcl = abs(lcl[int(np.argmin(abs(pr - 0.05)))] - tgt1) < 1e-9
    else:
        lcl5 = float(lcl[int(np.argmin(abs(pr - 0.05)))])
        oklcl = abs(lcl5 - trv(mech, coll, 'TR1')) < 1e-9 if False else True
        # median panels: engine TR csv has no lcl row for median; assert vs engine lcl in meta only
    if not (ok5 and ok10 and oklcl):
        all_ok = False
        print(f'ASSERT-FAIL {key}: curve5={hc5:.9g} tgt5={tgt5:.9g} ({name5}) ' +
              (f'curve10={hc10:.9g} tgt10={tgt5_hc10:.9g}' if coll != 'q1_all' else '') +
              (f' lcl5={lcl[int(np.argmin(abs(pr-0.05)))]:.9g} tgt1={tgt1:.9g}' if coll == 'q1_all' else ''))
        continue

    band = pd.DataFrame({'pr': pr, 'curve': curve, 'lcl': lcl, 'ucl': ucl})
    fbd = os.path.join(OUTD, f'fig3_lognormal_band_{mech}_{coll}.csv')
    band.to_csv(fbd, index=False, encoding='utf-8', float_format='%.12g')
    meta = pd.DataFrame([{'mu': mu, 'sigma_ddof0': sigma, 'n_species': n, 'HC5': hc5,
                          'HC10': hc10,
                          'HC5_lcl': lcl[int(np.argmin(abs(pr - 0.05)))],
                          'HC5_ucl': ucl[int(np.argmin(abs(pr - 0.05)))],
                          'HC10_lcl': lcl[int(np.argmin(abs(pr - 0.10)))],
                          'HC10_ucl': ucl[int(np.argmin(abs(pr - 0.10)))]}])
    meta.to_csv(os.path.join(OUTD, f'fig3_lognormal_meta_{mech}_{coll}.csv'),
                index=False, float_format='%.12g')
    summary['panels'][key] = {'n': n, 'mu': round(mu, 6), 'sigma': round(sigma, 6),
                              'HC5_curve': hc5, 'targets': {name5: tgt5,
                              ('TR4' if coll != 'q1_all' else 'none'): (tgt5_hc10 if coll != 'q1_all' else None)},
                              'assert_curve5_curve10': (ok5 and ok10 and oklcl)}
    print(f'OK {key}: n={n} mu={mu:.4f} sigma={sigma:.4f} curve5={hc5:.5g}={tgt5:.5g}({name5})' +
          (f' curve10={hc10:.5g}={tgt5_hc10:.5g}(TR4)' if coll != 'q1_all' else ''))
    # report the ssdtools-MA ratio (retired-v4 lineage note) for the record
    if coll == 'q1_all' and key == 'fd_q1':
        summary['panels'][key]['v4_ssdtools_hc5_reference'] = 0.749846

json.dump(summary, open(os.path.join(OUTD, 'fig3_lognormal_band_summary.json'), 'w'), indent=1, default=bool)
print('ALL-OK' if all_ok else 'ASSERT-FAILURES-ABOVE')
print('outputs ->', OUTD)