# -*- coding: utf-8 -*-
"""
step3_thresholds.py — [STEP 3] Tier(1~5)별 임계값 SSD 도출 + 3계층 기준값 + 최종 2x2
R 원본: 20260410_MP_criteria_2.Rmd 청크3(run_threshold_ssd_v4) + 청크6(Final 2x2)
산출물: Final_SSD_Criteria_1_300.csv / SSD_Management_Summary_with_N.csv
       + SSD 곡선용 예측 데이터( fig 코드용 )
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mp_pipeline_lib import *
from mp_pipeline_lib import _cdf
import numpy as np
import pandas as pd

INVITRO = "sod|cat|enzyme|biomarker|ros|molecular|gene|cellular|tissue"


def run_threshold_ssd(df_master, media, mech_type, tier, rep_cf):
    """R run_threshold_ssd_v4 1:1. tier 1~5 → 3계층 기준값(1-5000 Bio/Std·1-300)."""
    data = df_master.copy()
    if tier == 5:
        data = data[data["Duration_Type"] == "Acute"].copy()
        data["Conc_Food"] = data["Conc_Food"] * 10
        data["Conc_Trans"] = data["Conc_Trans"] * 10
    if tier in (3, 4):
        data = data[~data["Endpoint_classification"].astype(str)
                    .str.lower().str.contains(INVITRO, regex=True, na=False)]
    df_media = data[data["Media_Group"] == media]
    target_col = "Conc_Food" if mech_type == "Food" else "Conc_Trans"
    if len(df_media) == 0:
        return None
    # 종 대표값: T1/T2/T5 = Q1(25%), T3/T4 = median (R 1:1)
    if tier in (1, 2, 5):
        rep = df_media.groupby("학명")[target_col].quantile(0.25)
    else:
        rep = df_media.groupby("학명")[target_col].median()
    df_input = rep.reset_index()
    df_input.columns = ["학명", "Rep_Conc"]
    df_input["Alpha"] = df_media.groupby("학명")["Alpha"].first().values
    df_input = df_input[df_input["Rep_Conc"].notna() & (df_input["Rep_Conc"] > 0)]
    if len(df_input) < 5:  # R: n<5 → skip
        return None
    fits = ssd_fit_dists(df_input["Rep_Conc"], dists=("lnorm", "llogis"))
    if not fits:
        return None
    prop_val = 0.10 if tier == 4 else 0.05
    if tier == 1:
        hc = ssd_hc_modelavg(fits, percent=prop_val * 100, ci=True, nboot=100)
        raw_hc = hc.get("lcl", hc["est"])  # R: T1은 lcl(보수적 하한) 사용
    else:
        hc = ssd_hc_modelavg(fits, percent=prop_val * 100, ci=False)
        raw_hc = hc["est"]
    # 스케일 재조정: r_300 = 1-5000 중 1-300 구간 분율
    alpha_val = float(df_input["Alpha"].iloc[0])
    g = 1.0 - alpha_val
    r_300 = (300.0 ** g - 1.0 ** g) / (5000.0 ** g - 1.0 ** g)
    hc_bio_5000 = raw_hc
    hc_std_5000 = raw_hc / rep_cf if mech_type == "Food" else raw_hc
    hc_std_300 = hc_std_5000 * r_300
    tier_name = "Threshold 5 (Acute)" if tier == 5 else f"Threshold {tier}"
    return {"Media": media, "Mechanism": mech_type, "Threshold_Tier": tier_name,
            "N_Species": len(df_input), "CF_bio_Used": rep_cf,
            "HC5_Bio_1_5000": hc_bio_5000,
            "HC5_Std_Monitoring_1_5000": hc_std_5000,
            "HC5_Std_Monitoring_1_300": hc_std_300}


def main():
    print("[1/2] Unified_Master 로드 + 매체별 대표 CF 산출...")
    df = pd.read_csv(f"{OUT}/Unified_Master_Aligned_DB.csv", encoding="utf-8-sig")
    df["CF"] = df["Conc_Food"] / df["Conc_Trans"]
    cf_medians = df.groupby("Media_Group")["CF"].median()
    print("  media CF:", {k: round(v, 4) for k, v in cf_medians.items()})

    print("[2/2] Tier 1~5 × 매체 × 메커니즘 SSD 임계값 도출...")
    results = []
    for m in ["Seawater", "Sediment"]:
        m_cf = float(cf_medians.get(m, 1.0))
        for mech in ["Food", "Trans"]:
            for tier in range(1, 6):
                try:
                    res = run_threshold_ssd(df, m, mech, tier, m_cf)
                except Exception as e:
                    print(f"  ! {m}/{mech}/T{tier}: {e}")
                    res = None
                if res:
                    results.append(res)
    fw = pd.DataFrame(results).sort_values(["Media", "Mechanism", "Threshold_Tier"],
                                           ascending=[False, True, True])
    fw.to_csv(f"{OUT}/Final_SSD_Criteria_1_300.csv", index=False, encoding="utf-8-sig")
    print(fw.to_string(index=False))

    print("\n[3/3] 최종 관리 프레임워크 2x2 (Seawater/Sediment × Chronic/Acute, n 포함)...")
    summ_rows, curves = [], []
    for media_name in ["Seawater", "Sediment"]:
        for set_type in ["Chronic", "Acute"]:
            sub = df[df["Media_Group"] == media_name]
            if set_type == "Acute":
                sub = sub[~sub["Endpoint_classification"].astype(str)
                          .str.lower().str.contains(INVITRO, regex=True, na=False)]
            rep = sub.groupby("학명").agg(
                Conc=("Conc_Food", lambda s: float(s.quantile(0.25))),
                Alpha=("Alpha", "first")).reset_index()
            rep = rep[rep["Conc"] > 0]
            n_sp = len(rep)
            if n_sp < 5:
                continue
            fit = ssd_fit_dists(rep["Conc"], dists=("lnorm", "llogis"))
            hc5_bio = ssd_hc_modelavg(fit, percent=5.0)["est"]
            m_cf = float(cf_medians.get(media_name, 1.0))
            alpha_val = float(rep["Alpha"].iloc[0])
            g = 1.0 - alpha_val
            r_300 = (300.0 ** g - 1.0 ** g) / (5000.0 ** g - 1.0 ** g)
            hc5_std_300 = (hc5_bio / m_cf) * r_300
            unit = "particles/L" if media_name == "Seawater" else "particles/kg"
            summ = {"Media": media_name, "Type": set_type, "N_Species": n_sp,
                    "HC5_Bio_1_5000": hc5_bio, "HC5_Std_1_300": hc5_std_300}
            summ_rows.append(summ)
            # fig용: 종 포인트 + 모델평균 예측곡선
            rep_sorted = rep.sort_values("Conc").reset_index(drop=True)
            rep_sorted["y_prob"] = (np.arange(len(rep_sorted)) + 0.5) / len(rep_sorted)
            for _, rr in rep_sorted.iterrows():
                curves.append({"Media": media_name, "Type": set_type,
                               "Conc": rr["Conc"], "y_prob": rr["y_prob"], "kind": "point"})
            xs = np.geomspace(rep["Conc"].min() / 10, rep["Conc"].max() * 10, 120)
            w = {d: f for d, f in fit.items()}
            aic = {d: 2 * f["nll"] + 2 * f["k"] for d, f in fit.items()}
            mn = min(aic.values())
            wt = {d: math.exp(-0.5 * (aic[d] - mn)) for d in aic}
            tw = sum(wt.values())
            wt = {d: v / tw for d, v in wt.items()}
            for xv in xs:
                yv = sum(wt[d] * _cdf(d, xv, fit[d]["params"]) for d in wt)
                curves.append({"Media": media_name, "Type": set_type,
                               "Conc": xv, "y_prob": yv, "kind": "curve"})
            print(f"  {media_name} {set_type}: n={n_sp} HC5_bio={hc5_bio:.4f} "
                  f"HC5_std(1-300)={hc5_std_300:.4f} {unit}")
    pd.DataFrame(summ_rows).to_csv(f"{OUT}/SSD_Management_Summary_with_N.csv",
                                   index=False, encoding="utf-8-sig")
    pd.DataFrame(curves).to_csv(f"{OUT}/SSD_Curve_Points_forFig.csv",
                                index=False, encoding="utf-8-sig")
    print("\n✔ STEP3 완료 — Final_SSD_Criteria_1_300 / Summary / fig용 곡선 산출")


if __name__ == "__main__":
    main()