# -*- coding: utf-8 -*-
"""Fig.2 v9 acute pool rebuild (2026-10-13, LC50/EC50-only rule).

User-confirmed rule (verbatim basis):
  - Endpoint: LC50/EC50 only; NOEC/LOEC excluded (ToMEx 'Acute' label is NOT sufficient)
  - Size: 20-300 um, applied row-wise
  - Concentration: original particle-number concentration (particles/L)
  - No acute->chronic AF
  - Species aggregation: existing GM/MIN convention
  - Existing agreed eligibility (In vivo, exposure duration <=96h) applied
  - Left-censor '<' rows excluded (established policy)

Outputs (fig2_v9_tmp/):
  reps_C7fix.csv   : legacy 11-species list - censor-species(3) - Perinereis(>96h) = 7 spp,
                     with the Haliotis <=96h fix (20-min EC50 only; 14 d rows dropped)
  reps_P13keep.csv : all qualifying species (13), censor '>' rows kept as reported values
  reps_P13drop.csv : all qualifying species (13), censor '>' rows excluded
Prints per-species reps + bottom-5 for QA; exits non-zero if any pool empty.
"""
import numpy as np
import pandas as pd

from pathlib import Path
SRC = str(Path(__file__).resolve().parents[1] / "data" / "Domestic_Long.csv")
OUT = str(Path(__file__).resolve().parents[1] / "output" / "fig2_pools")

# Legacy 11-species acute pool (v7 era, pre-size-rule construct)
LEGACY11 = {
    "Aiptasiogeton eruptaurantia", "Artemia franciscana", "Brachionus koreanus",
    "Haliotis discus hannai", "Macrophthalmus japonicus", "Mesocentrotus nudus",
    "Monocorophium uenoi", "Mytilus galloprovincialis", "Neomysis awatschensis",
    "Perinereis linea", "Ruditapes philippinarum", "Tigriopus japonicus",
}
# legacy censor-species bucket excluded from C7 (v8B-style policy)
CENSOR_SP = {"Ruditapes philippinarum", "Monocorophium uenoi", "Mytilus galloprovincialis"}


def dur_days(s: str):
    s = str(s).replace(" ", "")
    if s.endswith("d"):
        return float(s[:-1])
    if s.endswith("h"):
        return float(s[:-1]) / 24.0
    if s.endswith("min"):
        return float(s[:-3]) / 1440.0
    return None


def main():
    dl = pd.read_csv(SRC, encoding="utf-8-sig")
    b = dl[(dl["Acute"].astype(str).str.strip() == "o")
           & dl["Endpoint_Type"].isin(["LC50", "EC50"])].copy()
    b["cen"] = b["Censor_Dir"].fillna("").astype(str).str.strip()
    b = b[b["cen"] != "<"]                                   # left-censor excluded (policy)
    b["Size_num"] = pd.to_numeric(b["Size_um"], errors="coerce")
    b = b[(b["Size_num"] >= 20) & (b["Size_num"] <= 300)]      # 20-300 um row-wise
    b["dur_d"] = b["Test_duration"].apply(dur_days)
    n_bad_dur = int(b["dur_d"].isna().sum())
    b = b[b["dur_d"] <= 4.0]                                  # <=96h eligibility
    assert not b.empty
    b["Particle_Conc"] = b["Particle_Conc"].astype(float)

    def reps(pool: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for sp, g in pool.groupby("Species_Name"):
            v = g["Particle_Conc"].values
            rows.append(dict(Tier="Acute", Species=sp, src="ours",
                             gm=float(np.exp(np.log(v).mean())),
                             mn=float(v.min())))
        out = pd.DataFrame(rows).sort_values("gm").reset_index(drop=True)
        out["metrics"] = np.nan
        out["reps"] = np.nan
        return out[["Tier", "Species", "src", "gm", "mn", "metrics", "reps"]]

    # (C7fix) legacy list - censor species - duration-fails
    c7 = b[b["Species_Name"].isin(LEGACY11 - CENSOR_SP)]
    c7 = c7[~c7["Species_Name"].isin(["Perinereis linea", "Monocorophium uenoi"])]
    r_c7 = reps(c7)

    r_keep = reps(b)                      # '>' rows kept as reported
    r_drop = reps(b[b["cen"] != ">"])     # '>' rows excluded

    for name, d in [("reps_C7fix", r_c7), ("reps_P13keep", r_keep), ("reps_P13drop", r_drop)]:
        path = f"{OUT}/{name}.csv"
        d.to_csv(path, index=False, encoding="utf-8-sig")
        chk = pd.read_csv(path)          # disk read-back verification
        assert len(chk) == len(d)

    print("rows after filters:", len(b), "| missing-duration rows (dropped):", n_bad_dur)
    print("endpoint classes:", sorted(b["Endpoint_classification"].dropna().astype(str).unique()))
    for name, d in [("C7fix", r_c7), ("P13keep", r_keep), ("P13drop", r_drop)]:
        print(f"\n=== {name} ({len(d)} spp) ===")
        show = d[["Species", "gm", "mn"]].copy()
        show["gm"] = show["gm"].map("{:,.2f}".format)
        show["mn"] = show["mn"].map("{:,.2f}".format)
        print(show.sort_values("mn").to_string(index=False))
        b5 = d.sort_values("mn").head(5)
        print("bottom5(MIN):", "; ".join(f"{r.Species}: {r.mn:,.2f}" for r in b5.itertuples()))


if __name__ == "__main__":
    main()