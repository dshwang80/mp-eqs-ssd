# -*- coding: utf-8 -*-
"""
step1_rescale.py — [STEP 1] 최종 취합본 → 엔드포인트 long 변환 → 입자 농도 표준화
R 원본: 20260410_MP_criteria.Rmd 청크 1(마스터 파이프라인) + 5(동기화 마스터) + 6(최종 통합)
산출물:
  MP_Acute_Endpoints_Rescaled_Final.csv / MP_Chronic_Endpoints_Rescaled_Final.csv
  MP_Cleaned_Final.csv (+ Acute/Chronic)
  SSD_Master_Acute.csv / SSD_Master_Chronic.csv          (This Study + ToMEx, 종별 최저 유효값)
  SSD_Master_Acute_Comparison.csv / _Chronic_Comparison.csv (Lowest/Q1/Median 3 대푯값)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mp_pipeline_lib import *
import numpy as np
import pandas as pd

os.makedirs(OUT, exist_ok=True)


def main():
    print("[1/4] 최종 취합본(xlsx) → long 변환 + 입자 환산/크기·Bio 보정...")
    lab = load_lab_final()
    print(f"  long rows: {len(lab)}  (종 수: {lab['학명'].nunique()})")

    # ---- 최종 데이터셋 저장 (Rmd#1 청크1 산출물 대응: Acute/Chronic 'o' 표시 기준)
    lab["CF_bio_alpha1.6"] = get_bio_correction(1.6)  # R 원본: 랩데이터 최종 보정에 alpha 1.6 사용
    lab["Final_Val"] = lab["Rescaled_Endpoint_Value"] * lab["CF_bio_alpha1.6"]

    keep = ["No", "Category", "Species_Name", "학명", "Life_stage", "Environment", "Stage",
            "Type", "Shape", "Size_um", "Size_Low", "Size_High", "Size_um_Midpoint",
            "Test_duration", "Acute", "Chronic", "Endpoint_classification",
            "Endpoint_Type", "Endpoint_Unit", "Endpoint_Value", "Media_Group", "Alpha",
            "CF_size", "CF_bio", "Rescaled_Endpoint_Value", "Final_Val", "Origin"]
    lab_out = lab[keep]
    lab_out[lab_out["Acute"] == "o"].to_csv(f"{OUT}/MP_Acute_Endpoints_Rescaled_Final.csv",
                                            index=False, encoding="utf-8-sig")
    lab_out[lab_out["Chronic"] == "o"].to_csv(f"{OUT}/MP_Chronic_Endpoints_Rescaled_Final.csv",
                                              index=False, encoding="utf-8-sig")
    lab_out.to_csv(f"{OUT}/MP_Cleaned_Final.csv", index=False, encoding="utf-8-sig")
    n_ac = int((lab_out["Acute"] == "o").sum())
    n_ch = int((lab_out["Chronic"] == "o").sum())
    print(f"  Acute('o'): {n_ac} / Chronic('o'): {n_ch} rows 저장")

    print("[2/4] ToMEx(3,924행) 국제표준 보정(질량→입자, alpha 1.6/퇴적물 3.25)...")
    tomex = pd.read_csv(TOMEX_CSV)
    tomex["학명"] = tomex["Species"].astype(str).str.strip()
    tomex["Dose"] = pd.to_numeric(
        tomex["Unaligned Dose Values"].astype(str).str.replace(r"[^0-9.\-]", "", regex=True),
        errors="coerce")
    tomex["Size"] = pd.to_numeric(
        tomex["Particle Length (μm)"].astype(str).str.replace(r"[^0-9.\-]", "", regex=True),
        errors="coerce")
    tomex["Dens"] = pd.to_numeric(tomex["Density (g/cm^3)"], errors="coerce")
    tomex["Alpha_Val"] = np.where(tomex["Environment"].astype(str).str.lower().str.contains("sediment"),
                                  3.25, 1.6)
    tomex = tomex[tomex["Dose"].notna() & (tomex["Dose"] > 0)
                  & tomex["Size"].notna() & tomex["Dens"].notna()].copy()
    tomex["Final_Val"] = [convert_tomex(m, s, d, a) for m, s, d, a in
                          zip(tomex["Dose"], tomex["Size"], tomex["Dens"], tomex["Alpha_Val"])]
    tomex = tomex[tomex["Final_Val"].notna() & (tomex["Final_Val"] > 0)]
    tomex["Source"] = "ToMEx Database"
    tomex["Origin"] = np.where(tomex["학명"].str.lower().isin([s.lower() for s in DOMESTIC_SPECIES]),
                               "Domestic (국내종)", "Foreign (국외종)")
    tomex["Environment"] = np.where(tomex["Environment"].astype(str).str.lower().str.contains("fresh"),
                                    "Freshwater (담수)", "Marine (해수)")
    tomex["Duration_Type"] = np.where(tomex["Acute/Chronic"].astype(str).str.lower().str.contains("acute"),
                                      "급성", "만성")
    print(f"  ToMEx 유효 rows: {len(tomex)} (종 {tomex['학명'].nunique()})")

    print("[3/4] This Study + ToMEx 통합 → 종별 최저 유효값(0 제외) SSD 마스터 구축...")
    lab_ss = lab_out[["Origin", "Environment", "Stage", "Acute", "Chronic",
                      "학명", "Final_Val"]].copy()
    lab_ss["Source"] = "This Study"
    lab_ss["Duration_Type"] = np.where(lab_ss["Acute"] == "o", "급성",
                                       np.where(lab_ss["Chronic"] == "o", "만성", "기타"))
    tomex_ss = tomex[["Source", "Origin", "Environment", "Duration_Type", "학명", "Final_Val"]].copy()
    comb = pd.concat([lab_ss[["Source", "Origin", "Environment", "Duration_Type", "학명", "Final_Val"]],
                      tomex_ss], ignore_index=True)
    comb = comb[comb["Final_Val"] > 0]
    master = (comb.groupby(["Source", "Origin", "Environment", "Duration_Type", "학명"], as_index=False)
              .agg(Lowest_Value=("Final_Val", "min"),
                   Q1_Value=("Final_Val", lambda s: float(np.quantile(s, 0.25))),
                   Median_Value=("Final_Val", "median")))
    master[master["Duration_Type"] == "급성"].to_csv(f"{OUT}/SSD_Master_Acute.csv",
                                                     index=False, encoding="utf-8-sig")
    master[master["Duration_Type"] == "만성"].to_csv(f"{OUT}/SSD_Master_Chronic.csv",
                                                     index=False, encoding="utf-8-sig")
    master[master["Duration_Type"] == "급성"].to_csv(f"{OUT}/SSD_Master_Acute_Comparison.csv",
                                                     index=False, encoding="utf-8-sig")
    master[master["Duration_Type"] == "만성"].to_csv(f"{OUT}/SSD_Master_Chronic_Comparison.csv",
                                                     index=False, encoding="utf-8-sig")
    for d in ("급성", "만성"):
        sub = master[master["Duration_Type"] == d]
        dom = sub[sub["Origin"].str.contains("Domestic")]
        print(f"  SSD_Master {d}: {len(sub)}종 (국내종 {len(dom)})")

    print("[4/4] Median/Q1 대푯값 검증표 (종별 N>1 그룹)...")
    grp = (comb.groupby(["학명", "Duration_Type"], as_index=False)
           .agg(N=("Final_Val", "size"), Min=("Final_Val", "min"),
                Q1=("Final_Val", lambda s: float(np.quantile(s, 0.25))),
                Median=("Final_Val", "median"), Max=("Final_Val", "max")))
    grp[grp["N"] > 1].to_csv(f"{OUT}/Representative_Value_Comparison_Check.csv",
                             index=False, encoding="utf-8-sig")
    print("\n✔ STEP1 완료 — Rescaled/Cleaned/SSD_Master 산출")


if __name__ == "__main__":
    main()