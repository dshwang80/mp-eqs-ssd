# -*- coding: utf-8 -*-
"""
step2_align.py — [STEP 2] 3소스(Lab/검색/ToMEx) AF 적용 + Koelmans 메커니즘 정렬 통합DB
R 원본: 20260410_MP_criteria_2.Rmd 청크1(수정 마스터) — Unified_Master_Aligned_DB.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mp_pipeline_lib import *
import numpy as np
import pandas as pd

SEARCH_SHEETS = ["2022년(1차년도)", "2023년(2차년도)", "2024년(3차년도)", "2025년(4차년도)"]
EP_ORDER = [("NOEC", "NOEC"), ("EC10", "EC10"), ("LC10", "LC10"),
            ("LOEC", "LOEC"), ("EC50", "EC50"), ("LC50", "LC50")]


def load_search_long():
    """검색DB 4개 연도 시트 → long. (particles/L 열, ug/mL 열) 쌍 — particles/L 우선."""
    frames = []
    s_cols = ["No", "Category", "Species_Name", "Life_stage", "Environment",
              "Type", "Shape", "Size_um", "Test_duration", "Target_Cell",
              "Acute", "Chronic", "Endpoint_classification",
              "EC50_p", "EC50_u", "LC50_p", "LC50_u", "NOEC_p", "NOEC_u",
              "LOEC_p", "LOEC_u", "Data_source", "DOI"]
    for sh in SEARCH_SHEETS:
        try:
            df = pd.read_excel(SEARCH_XLSX, sheet_name=sh, skiprows=3, header=None)
        except Exception:
            continue
        df = df.iloc[:, :len(s_cols)]
        df.columns = s_cols[: df.shape[1]]
        frames.append(df)
    raw = pd.concat(frames, ignore_index=True)
    ep_pairs = [("NOEC", "NOEC_p", "NOEC_u"), ("EC10", "EC10_p", "EC10_u"),
                ("LC10", "LC10_p", "LC10_u"), ("LOEC", "LOEC_p", "LOEC_u"),
                ("EC50", "EC50_p", "EC50_u"), ("LC50", "LC50_p", "LC50_u")]
    rows = []
    for _, r in raw.iterrows():
        for ep, cp, cu in ep_pairs:
            # R coalesce 순서: particles/L 우선, 부재 시 ug/mL
            v = clean_tox_value(r.get(cp))
            unit = "particles/L"
            if not np.isfinite(v):
                v = clean_tox_value(r.get(cu))
                unit = "ug/mL"
            if np.isfinite(v):
                rows.append({
                    "Species_Name": str(r["Species_Name"]).strip() if pd.notna(r["Species_Name"]) else "",
                    "Environment": r["Environment"], "Type": r["Type"], "Shape": r["Shape"],
                    "Size_um": r["Size_um"], "Test_duration": r["Test_duration"],
                    "Acute": r["Acute"], "Chronic": r["Chronic"],
                    "Endpoint_classification": r["Endpoint_classification"],
                    "Endpoint_Type": ep, "Endpoint_Value": v, "Endpoint_Unit": unit,
                    "Reference": r.get("Data_source"), "DOI": r.get("DOI")})
    return pd.DataFrame(rows)


def apply_paper_cleanup(df, src):
    """STEP2 국외DB 정제 — Mehinto 2022·Hampton 2025(ToMEx 2.0) SSD 투입 전
    데이터 수준 제외 규칙 (논문 Methods 그대로):
      ① HONEC 제외   — "Highest observed no effect concentrations (HONEC) were
         excluded from our analyses due to their limited reliability"
      ② 조류(Algae) 제외 — "Studies using algal species were excluded as
         ingestion is not a plausible mechanism" (food dilution 메커니즘 기준)
      ※ 문헌 수준 제외(>5mm·현장·독물동태·in vitro 배제)는 원 DB 스크리닝에서
         이미 적용된 사항 — ToMEx 2.0 풀 export에는 해당 없음(max 3,000㎛, In
         vitro 157행은 threshold 분석 계층에 사용되므로 유지).
    """
    n0 = len(df)
    if src == "tomex":
        em = df["Effect Metric"].astype(str).str.strip()
        grp = df["Organism Group"].astype(str).str.strip()
        keep = (em != "HONEC") & (grp != "Algae")
    else:  # search: HONEC 없음(실측 LOEC/NOEC/LC50/EC50뿐), 조류 2종만 대상
        sp = df["학명"].astype(str).str.strip()
        keep = ~sp.isin({"Skeletonema costatum", "Karenia mikimotoi"})
    out = df[keep].copy()
    print(f"  [{src}] 논문 정제 규칙 적용: {n0} → {len(out)}행 "
          f"(HONEC/조류 제외 {n0 - len(out)}행)")
    return out


def main():
    print("[1/4] 입 크기(Mouth Size) 마스터: ToMEx 'Estimated Maximum Ingestible Size' + 수동 보강")
    tomex_raw = pd.read_csv(TOMEX_CSV)
    col_ms = "Estimated Maximum Ingestible Size (mm)"
    ms = tomex_raw[tomex_raw[col_ms].notna()].copy()
    ms["학명"] = ms["Species"].astype(str).str.strip()
    ms["Mouth_Size_um"] = pd.to_numeric(
        ms[col_ms].astype(str).str.replace(r"[^0-9.\-]", "", regex=True), errors="coerce") * 1000
    mouth_db = (ms[ms["Mouth_Size_um"].notna() & (ms["Mouth_Size_um"] > 0)]
                .groupby("학명", as_index=False)["Mouth_Size_um"].median())
    manual = pd.DataFrame({"학명": ["Haliotis discus hannai", "Paralichthys olivaceus", "Mytilus edulis"],
                           "Mouth_Size_um": [500.0, 1000.0, 150.0]})
    mouth_db = pd.concat([mouth_db, manual]).drop_duplicates(subset="학명", keep="first")
    print(f"  mouth size DB: {len(mouth_db)}종")

    print("[2/4] Lab(This Study): AF 적용 + Koelmans food/trans 정렬...")
    lab = pd.read_csv(f"{OUT}/MP_Cleaned_Final.csv", encoding="utf-8-sig")
    domestic_list = set(lab["학명"].str.strip().str.lower())
    lab = lab.merge(mouth_db, on="학명", how="left")
    lab["Mouth_Size_um"] = lab["Mouth_Size_um"].fillna(100.0)
    lab["Duration_Type"] = np.where(lab["Acute"].astype(str).str.strip() == "o", "Acute", "Chronic")
    lab["Chronic_NOEC_Val"] = [apply_af(v, ep, dt) for v, ep, dt in
                               zip(lab["Endpoint_Value"], lab["Endpoint_Type"], lab["Duration_Type"])]
    lab["Conc_Food"] = [align_koelmans(c, l, h, a, "food", m) for c, l, h, a, m in
                        zip(lab["Chronic_NOEC_Val"], lab["Size_Low"], lab["Size_High"],
                            lab["Alpha"], lab["Mouth_Size_um"])]
    lab["Conc_Trans"] = [align_koelmans(c, l, h, a, "trans", m) for c, l, h, a, m in
                         zip(lab["Chronic_NOEC_Val"], lab["Size_Low"], lab["Size_High"],
                             lab["Alpha"], lab["Mouth_Size_um"])]
    lab = lab[lab["Conc_Food"].notna() & (lab["Conc_Food"] > 0)].copy()
    lab["Source_Tag"] = "This Study"
    lab["Origin"] = "Domestic (국내종)"
    lab["Media_Group"] = np.where(lab["Environment"].astype(str).str.lower().str.contains("sediment"),
                                  "Sediment", "Seawater")
    print(f"  lab rows: {len(lab)}")

    print("[3/4] 검색DB(Literature) + ToMEx 처리...")
    search = load_search_long()
    # 국외DB 정제 — 논문 규칙(조류 제외, HONEC은 원본에 없음) 적용
    search_tmp = search.copy()
    search_tmp["학명"] = search_tmp["Species_Name"].astype(str).str.strip()
    search = apply_paper_cleanup(search_tmp, "search")
    search["학명"] = search["Species_Name"].str.strip()
    search = search.merge(mouth_db, on="학명", how="left")
    search["Mouth_Size_um"] = search["Mouth_Size_um"].fillna(100.0)
    search["Media_Group"] = search["Environment"].map(lambda e: media_group(e))
    search["Alpha"] = search["Media_Group"].map(media_alpha)
    search["Duration_Type"] = np.where(search["Acute"].astype(str).str.strip() == "o", "Acute", "Chronic")
    search["Chronic_NOEC_Val"] = [apply_af(v, ep, dt) for v, ep, dt in
                                  zip(search["Endpoint_Value"], search["Endpoint_Type"], search["Duration_Type"])]
    si = search["Size_um"].map(parse_size_info)
    search["Size_Low"] = [x[0] for x in si]
    search["Size_High"] = [x[1] for x in si]
    search["Conc_Food"] = [align_koelmans(c, l, h, a, "food", m) for c, l, h, a, m in
                           zip(search["Chronic_NOEC_Val"], search["Size_Low"], search["Size_High"],
                               search["Alpha"], search["Mouth_Size_um"])]
    search["Conc_Trans"] = [align_koelmans(c, l, h, a, "trans", m) for c, l, h, a, m in
                            zip(search["Chronic_NOEC_Val"], search["Size_Low"], search["Size_High"],
                                search["Alpha"], search["Mouth_Size_um"])]
    search = search[search["Conc_Food"].notna() & (search["Conc_Food"] > 0)].copy()
    search["Source_Tag"] = "Literature/Search"
    search["Origin"] = np.where(search["학명"].str.strip().str.lower().isin(domestic_list),
                                "Domestic (국내종)", "Foreign (국외종)")
    print(f"  search rows: {len(search)}")

    t = tomex_raw.copy()
    # 국외DB 정제 — 논문 규칙(HONEC·조류 제외) 적용
    t = apply_paper_cleanup(t, "tomex")
    t["학명"] = t["Species"].astype(str).str.strip()
    t = t.merge(mouth_db, on="학명", how="left")
    t["Mouth_Size_um"] = t["Mouth_Size_um"].fillna(100.0)
    env = t["Environment"].astype(str).str.lower()
    t["Media_Group"] = np.where(env.str.contains("sediment")
                                | t["Exposure Route"].astype(str).str.lower().str.contains("sediment"),
                                "Sediment", "Seawater")
    t["Alpha"] = t["Media_Group"].map(media_alpha)
    t["raw_val"] = pd.to_numeric(
        t["Unaligned Dose Values"].astype(str).str.replace(r"[^0-9.eE+\-]", "", regex=True),
        errors="coerce")
    t["is_mass"] = t["Dose Metric"].astype(str).str.lower().str.contains("mg|mass|g")
    t["part_size"] = pd.to_numeric(
        t["Particle Length (μm)"].astype(str).str.replace(r"[^0-9.\-]", "", regex=True), errors="coerce")
    t["dens"] = pd.to_numeric(t["Density (g/cm^3)"], errors="coerce")
    # R: raw_val / (dens * 4/3π((size/2)/10000)^3 * 1000) — size μm에서 mm 반경 환산 등가식
    def tomex_count(row):
        if (not row["is_mass"]) or pd.isna(row["part_size"]) or pd.isna(row["dens"]):
            return row["raw_val"]
        if row["part_size"] <= 0 or row["dens"] <= 0 or row["raw_val"] <= 0:
            return np.nan
        radius = (row["part_size"] / 2) / 10000  # R 원본 수식 그대로 (μm→cm)
        return row["raw_val"] / (row["dens"] * (4 / 3) * math.pi * radius ** 3 * 1000)
    t["raw_count"] = t.apply(tomex_count, axis=1)
    t["Duration_Type"] = np.where(t["Acute/Chronic"].astype(str).str.lower().str.contains("acute"),
                                  "Acute", "Chronic")
    t["EP_Type"] = t["Effect Metric"].astype(str)
    t["Chronic_NOEC_Val"] = [apply_af(v, ep, dt) for v, ep, dt in
                             zip(t["raw_count"], t["EP_Type"], t["Duration_Type"])]
    t["Conc_Food"] = [align_koelmans(c, s * 0.9, s * 1.1, a, "food", m) for c, s, a, m in
                      zip(t["Chronic_NOEC_Val"], t["part_size"], t["Alpha"], t["Mouth_Size_um"])]
    t["Conc_Trans"] = [align_koelmans(c, s * 0.9, s * 1.1, a, "trans", m) for c, s, a, m in
                       zip(t["Chronic_NOEC_Val"], t["part_size"], t["Alpha"], t["Mouth_Size_um"])]
    t = t[t["Conc_Food"].notna() & (t["Conc_Food"] > 0)].copy()
    t["Endpoint_classification"] = t["Endpoint"].astype(str)
    t["Source_Tag"] = "ToMEx DB"
    t["Origin"] = np.where(t["학명"].str.strip().str.lower().isin(domestic_list),
                           "Domestic (국내종)", "Foreign (국외종)")
    print(f"  tomex rows: {len(t)}")

    print("[4/4] 통합 저장 → Unified_Master_Aligned_DB.csv")
    cols = ["학명", "Origin", "Media_Group", "Alpha", "Endpoint_classification",
            "Duration_Type", "Mouth_Size_um", "Conc_Food", "Conc_Trans", "Source_Tag"]
    unified = pd.concat([lab[cols], search[cols], t[cols]], ignore_index=True)
    unified.to_csv(f"{OUT}/Unified_Master_Aligned_DB.csv", index=False, encoding="utf-8-sig")
    print(f"  unified rows: {len(unified)} (종 {unified['학명'].nunique()})")
    print("\n✔ STEP2 완료 — Unified_Master_Aligned_DB 산출")


if __name__ == "__main__":
    main()