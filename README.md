# mp-eqs-ssd

**Toxicity data processing and SSD analysis reproduction code (Figs. 1–3)** — microplastic
marine water-quality criteria (MP_EQS study).

This repository provides the analysis code and the **processed (derived) input tables**
behind Figures 1–3 of the manuscript (every HC5 / TR (threshold realization) value used
in those figures). With the tables supplied here, the three figures can be **regenerated**
and selected SSD computations can be **recomputed** (see Reproduction). It is **not** a
full-reproduction package from raw records: the third-party ToMEx export and some upstream
size-representativeness intermediates (consumed by the Fig. 2 HC5 bootstrap engine) are not
included. Also out of scope: the manuscript text, **Fig. 4 (national monitoring application)
and the monitoring-data analyses**, and the Supplementary Materials.

## Scope

| Output | What it shows | Reproduced by |
|---|---|---|
| Fig. 1 | Particle-size normalization paired assessment: (a) predicted vs. observed effect concentrations (particles/L) under two exponent settings (α = 2.07 seawater-associated, α = 1.60 generic); (b) empirical size–effect power-law description and reference (environmental) particle-size distributions, `p(L) ∝ L^-α` (20–300 µm window) | `code/fig1_empirical_test_v23.py` (+ helper `code/track1_fourpanel_applicability.py`) |
| Fig. 2 | Conventional pooled SSDs (GM vs. minimum-value pools; acute/chronic) with HC5 point estimates and 95% CIs, rendered from the locked ssdtools model-averaged fit values | `code/fig2_v9_pools.py`, `code/ssd_pooled_v4_merge.py`, `code/ssd_pooled_v4_fix1000.py` → `code/fig2_ssdtools_render_v12.R` |
| Fig. 3 | Normalized, mechanism-specific single-lognormal SSDs with derived threshold realizations TR1–TR4 — four panels by uptake mechanism × realization level: (a) FD–Q1, (b) TT–Q1, (c) FD–median, (d) TT–median; size-class records (T = 12, 20, 300 µm) pooled within each panel | `code/ssd_tomex_alpha207.py`, `code/fig3_lognormal_band_v1.py`, `code/fig3_recompute_check.py`, `code/fig3_panel_sizeclass_v2.py` → `code/fig3_lognormal_render_v3.R` |

Toxicity-data processing chain (`code/step1_rescale.py` → `step2_align.py` →
`step3_thresholds.py`, built on `code/mp_pipeline_lib.py`) rescales and aligns the
experimental endpoints and derives the SSD input database.

## Repository layout

```
code/          13 Python + 2 R scripts (tables above)
data/          processed, figure-ready CSV/JSON used as locked inputs
data/raw/      3 self-compiled Excel workbooks of the study toxicity database
verification/  locked value tables, render audit logs, published figures
```

## Data

- `data/raw/` contains the three **self-compiled** Excel workbooks used to assemble the
  toxicity database (curated by the authors; safe to redistribute).
- `data/` contains the **processed input tables** for every figure script, so Figs. 1–3 can
  be regenerated from this repository alone (tier A below), except where noted.
- **ToMEx** (Toxicity of Microplastics Exploration database; Thornton Hampton *et al.* 2022,
  2025): the raw export is **not redistributed** (third-party database). Two entry points
  expect it as `data/ToMEx_Exploration2026-03-24.csv`:
  - `code/ssd_tomex_alpha207.py` (Fig. 3 TR source engine — recomputes the size-class-level
    thresholds from raw records),
  - `code/step2_align.py` (unified DB build from raw records).
  Users with their own ToMEx export can place it at that path; all ToMEx-derived
  intermediates needed for the Fig. 2–3 **renders** are already included in `data/`.
- `code/ssd_pooled_v4_merge.py` / `_fix1000.py` (Fig. 2 HC5 bootstrap engine) additionally
  expect upstream size-representativeness intermediates **not shipped here**; the locked
  `data/ssd_pooled_hc5_v9.csv` / `ssd_pooled_reps_v9.csv` are what the renderer consumes.
  Fig. 2 regeneration therefore means rendering from the locked tables, not re-running the
  HC5 bootstrap engine itself.

## Reproduction (run from the repository root)

Two levels — figure regeneration from the shipped processed tables, and recomputation from
raw records that requires data not included here:

```bash
# --- A) Regenerate figures / check from shipped processed inputs (runs as-is) ---
python3 code/fig1_empirical_test_v23.py          # Fig. 1 -> output/fig1_v23_empirical_test.png
python3 code/fig2_v9_pools.py                    # Fig. 2 pool assembly -> output/fig2_pools/
Rscript code/fig2_ssdtools_render_v12.R          # Fig. 2 render (reads locked ssd_pooled_hc5_v9/reps)
python3 code/fig3_lognormal_band_v1.py           # Fig. 3 band/curve build -> output/fig3_lognormal_build/
python3 code/fig3_recompute_check.py             # Fig. 3 SSD recomputation check (agreement 1e-9)
Rscript code/fig3_lognormal_render_v3.R          # Fig. 3 render -> output/fig3_lognormal_v3.png

# --- B) Recompute-level engines (require data NOT included here) ---
python3 code/ssd_pooled_v4_merge.py              # Fig. 2 HC5 bootstrap engine -> upstream intermediates (not shipped)
python3 code/ssd_pooled_v4_fix1000.py            # Fig. 2 HC5 fix -> same upstream intermediates
python3 code/ssd_tomex_alpha207.py               # Fig. 3 TR source engine -> needs ToMEx export
python3 code/fig3_panel_sizeclass_v2.py          # Fig. 3 size-class adjudication -> needs ToMEx export

# --- C) Toxicity-data processing chain (full rebuild) ---
python3 code/step1_rescale.py                    # runs as-is (reads data/raw xlsx) -> output/py_out/
python3 code/step2_align.py                      # unified DB build -> needs ToMEx export + search-DB xlsx
python3 code/step3_thresholds.py                 # criteria tables -> needs step2 output
```

The shipped `verification/` tables let every figure be regenerated without tier B, while
tier B/C document honestly which recomputation steps need the (third-party or upstream)
files that are not redistributed.

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