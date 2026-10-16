# mp-eqs-ssd

**Toxicity data processing and SSD analysis reproduction code (Figs. 1–3)** — microplastic
marine water-quality criteria (MP_EQS study).

This repository reproduces Figures 1–3 of the manuscript, i.e. the toxicity-data
processing chain and the species-sensitivity-distribution (SSD) analyses behind every
HC5 / TR (threshold ratio) value used in those figures. It is **not** a full-reproduction
package for the entire manuscript: manuscript text, SI figures (Figs. 4/S1), and
monitoring-data analyses are out of scope.

## Scope

| Output | What it shows | Reproduced by |
|---|---|---|
| Fig. 1 | Size-resolved cross-comparison of microplastic effect data (parity of acute/chronic pool sensitivities; power-law shape check `p(L) ∝ L^-α`) | `code/fig1_empirical_test_v23.py` (+ helper `code/track1_fourpanel_applicability.py`) |
| Fig. 2 | Pooled SSD curves (GM vs. minimum-value pools, acute/chronic) → HC5 | `code/fig2_v9_pools.py`, `code/ssd_pooled_v4_merge.py`, `code/ssd_pooled_v4_fix1000.py` → `code/fig2_ssdtools_render_v12.R` |
| Fig. 3 | Single-lognormal SSDs per uptake mechanism (food dilution / translocation) and nominal particle size (T = 12/20/300 µm) → TR1–TR4 | `code/ssd_tomex_alpha207.py`, `code/fig3_lognormal_band_v1.py`, `code/fig3_recompute_check.py`, `code/fig3_panel_sizeclass_v2.py` → `code/fig3_lognormal_render_v3.R` |

Toxicity-data processing chain (`code/step1_rescale.py` → `step2_align.py` →
`step3_thresholds.py`, built on `code/mp_pipeline_lib.py`) rescales and aligns the
experimental endpoints and derives the SSD input database.

## Repository layout

```
code/          13 Python + 2 R scripts (tables above)
data/          derived, figure-ready CSV/JSON used as locked inputs
data/raw/      3 self-compiled Excel workbooks of the study toxicity database
verification/  locked value tables, render audit logs, published figures
```

## Data

- `data/raw/` contains the three **self-compiled** Excel workbooks used to assemble the
  toxicity database (curated by the authors; safe to redistribute).
- `data/` contains the derived CSV/JSON inputs for every figure script, so Figs. 1–3 can be
  regenerated from this repository alone, except where noted below.
- **ToMEx** (Toxicity of Microplastics Exploration database; Koelmans *et al.*): the raw
  ToMEx export is **not redistributed** (third-party database). Two entry points expect it
  as `data/ToMEx_Exploration2026-03-24.csv`:
  - `code/ssd_tomex_alpha207.py` (Fig. 3 TR source engine),
  - `code/step2_align.py` (unified DB build).
  Users with their own ToMEx export can place it at that path; all ToMEx-derived
  intermediates needed for the Fig. 2–3 renders are already included in `data/`.
- `code/ssd_pooled_v4_merge.py` / `_fix1000.py` (Fig. 2 HC5 bootstrap engine) additionally
  expect upstream size-representativeness intermediates not shipped here; the locked
  `data/ssd_pooled_hc5_v9.csv` / `ssd_pooled_reps_v9.csv` are what the renderer consumes.

## Reproduction (run from the repository root)

```bash
# Fig. 1 (out-of-the-box)
python3 code/fig1_empirical_test_v23.py          # -> output/fig1_v23_empirical_test.png

# Fig. 2 (out-of-the-box render; HC5 engine needs upstream chain, see Data notes)
python3 code/fig2_v9_pools.py                    # -> output/fig2_pools/
Rscript code/fig2_ssdtools_render_v12.R          # -> output/fig2_ssdtools_v12.png

# Fig. 3 band/TR computation (out-of-the-box from locked inputs)
python3 code/fig3_lognormal_band_v1.py           # -> output/fig3_lognormal_build/
python3 code/fig3_recompute_check.py             # cross-check vs ssdtools inputs
python3 code/fig3_panel_sizeclass_v2.py          # -> output/fig3_sizeclass/ (needs ToMEx export)
Rscript code/fig3_lognormal_render_v3.R          # -> output/fig3_lognormal_v3.png

# Toxicity-data processing chain
python3 code/step1_rescale.py                    # -> output/py_out/ (from data/raw xlsx)
python3 code/step2_align.py                      # needs ToMEx export + search-DB xlsx
python3 code/step3_thresholds.py                 # -> output/py_out/ criteria tables
```

## Locked values (verbatim from `verification/`)

HC5 (pooled SSD, Fig. 2; `HC5_lock_ssd_pooled_hc5_v9.csv`):

| Tier | Track | n species | Pool GM | HC5 |
|---|---|---|---|---|
| Acute | GM | 13 | 17,224,128.83 | 161,305.9 |
| Acute | MIN | 13 | — | 2,435.6 |
| Chronic | GM | 15 | 436,199.85 | 156.008 |
| Chronic | MIN | 15 | — | 60.492 |

TR (Fig. 3; `TR_lock_abcd_merged_T12_20_300.csv`, quoted `particles_L` column; a
mass proxy `mg_L_proxy` is included in the file — definitions and units are given in the
manuscript):

| Mechanism | TR1 HC5 95% CI lower (q1_all) | TR2 HC5 point (q1_all) | TR3 HC5 point (median) | TR4 HC10 point (median) |
|---|---|---|---|---|
| Food dilution | 0.025480689 | 1.069964859 | 11.52677653 | 69.43953824 |
| Translocation | 0.691025531 | 39.14279487 | 144.8661936 | 1400.640967 |

## Requirements

- Python ≥ 3.10: `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl`
- R with `ssdtools`, `ggplot2`, `patchwork`, `ggrepel`
  (internal runs: R 4.5.2 + ssdtools 2.6.0; older R/ssdtools combinations were unstable —
  use a recent pair)
- `output/` holds regenerated intermediates and is gitignored; figure scripts write
  rendered PNGs + `*_audit.json` logs there.

## Verification

`verification/` bundles the exact audited artifacts behind the published figures:
`*_audit.json` render audits, the Fig. 3 render log, the band-summary JSON, the HC5 and
TR lock CSVs, the three published figure PNGs, and `md5_manifest.txt` with per-file
checksums of the released tree. `code/fig3_recompute_check.py` independently re-derives
the Fig. 3 lognormal fits from the shipped SSD inputs and reports agreement
in `verification/fig3_verification_report.json` (tolerance 1e-9).

## License

MIT — see `LICENSE`.