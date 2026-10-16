# -*- coding: utf-8 -*-
"""
mp_pipeline_lib.py — MP_EQS 공용 라이브러리 (R 코드 파이썬 변환)
원본: 20260410_MP_criteria.Rmd / _2.Rmd (R → Python 1:1 변환, 2026-09-01)
데이터 소스: [실험독성자료 1] 최종 취합본-260805-f.xlsx (V5→최종본 교체)

R 함수 → 파이썬 함수 대응표:
  clean_tox_value()            -> clean_tox_value()
  parse_size_info()            -> parse_size_info()
  get_powerlaw_cf()            -> get_powerlaw_cf()
  get_bio_correction()         -> get_bio_correction()
  convert_tomex_accurate()     -> convert_and_rescale_standard() / convert_tomex()
  apply_af()                   -> apply_af()
  align_koelmans()             -> align_koelmans()
  R quantile(type=7)           -> numpy.quantile (기본 linear = type 7 동일)
  ssdtools ssd_fit_dists/ssd_hc-> fit_ssd() / ssd_hc_modelavg() (scipy MLE + AICc 가중 평균)
"""
import re
import math
import warnings
import numpy as np
import pandas as pd
from scipy import stats, optimize

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------- 경로 상수
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / "data")
OUT = str(ROOT / "output" / "py_out")
CODE = str(ROOT / "code")

FINAL_XLSX = f"{DATA}/raw/[실험독성자료 1] MP 생태 독성 및 거동 DB (최종 취합본)-260805-f.xlsx"
SEARCH_XLSX = f"{DATA}/raw/[검색데이터] MP 생태독성 및 거동 DB-2025-f.xlsx"
TOMEX_CSV = f"{DATA}/ToMEx_Exploration2026-03-24.csv"
MP_SIZE_CSV = f"{DATA}/MP_size.csv"

ENDPOINT_PAIRS = [  # (이름, ㎎/L 열, ㎎/kg 열) — 1단계 MP 시트 28열 기준
    ("EC50", "EC50_L", "EC50_k"), ("EC10", "EC10_L", "EC10_k"),
    ("LC50", "LC50_L", "LC50_k"), ("LC10", "LC10_L", "LC10_k"),
    ("NOEC", "NOEC_L", "NOEC_k"), ("LOEC", "LOEC_L", "LOEC_k"),
    ("MATC", "MATC_L", "MATC_k"),
]
EP_REGEX = "sod|cat|enzyme|biomarker|ros|oxidative|molecular|gene|cellular|tissue|expression|세포독성|효소"

# [R 원본 유지] 국내종 하드코드 리스트 (Rmd#1 '동기화 마스터' 청크)
DOMESTIC_SPECIES = [
    "Brachionus koreanus", "Tigriopus japonicus", "Oryzias latipes",
    "Sebastes schlegeli", "Venerupis philippinarum", "Monodonta labio",
    "Skeletonema costatum", "Neomysis awatschensis", "Paracyclopina nana",
]

# ================================================================ [함수 0]
# 독성 수치 정밀 클리닝 (Rmd#2 clean_tox_value 1:1)
# ================================================================
def clean_tox_value(val_str):
    """한글·괄호·부등호·공백 제거 → 범위형(–,~) 첫 숫자 추출. NA/빈값 → NaN."""
    if val_str is None or (isinstance(val_str, float) and math.isnan(val_str)):
        return np.nan
    s = str(val_str)
    s = re.sub(r"[가-힣]+|[()]|[><]|\s", "", s)
    s = s.replace("–", "-").replace("~", "-")
    nums = re.findall(r"[0-9.]+", s)
    if nums:
        try:
            v = float(nums[0])
            return v if v > 0 else np.nan
        except ValueError:
            return np.nan
    return np.nan


# ================================================================ [함수 2]
# 크기 텍스트 파싱 (Rmd#1 parse_size_info 1:1)
# ================================================================
SKIP_TOKENS = {"-", "mixture", "particle", "leachate", "mixed", "nan", ""}

def parse_size_info(size_str):
    """'45-75', '>150', 3000 등 → (low, high, mid=sqrt(l*h)). mm 단위 자동 환산."""
    if size_str is None or str(size_str).strip().lower() in SKIP_TOKENS:
        return (np.nan, np.nan, np.nan)
    s = str(size_str).replace("–", "-").replace(" ", "")
    multi = 1000.0 if re.search(r"㎜|mm", s, re.I) else 1.0
    vals = re.findall(r"[0-9.]+", s)
    if not vals:
        return (np.nan, np.nan, np.nan)
    try:
        nums = [float(v) for v in vals]
    except ValueError:
        return (np.nan, np.nan, np.nan)
    if "-" in s and len(nums) >= 2:
        low, high = nums[0] * multi, nums[1] * multi
    elif ">" in s:
        low, high = nums[0] * multi, 5000.0
    else:
        low = high = nums[0] * multi
    if low <= 0:
        low = 0.1
    mid = math.sqrt(low * high)
    return (low, high, mid)


# ================================================================ [함수 1]
# 파워-로 크기 보정 계수 (Rmd#1 get_powerlaw_cf 1:1)
# ================================================================
def get_powerlaw_cf(l_low, l_high, alpha, std_min=1.0, std_max=5000.0):
    """질량·입자 농도를 표준 크기대(1–5000㎛)로 재조정하는 적분 기반 CF."""
    if (l_low is None or l_high is None or (isinstance(l_low, float) and math.isnan(l_low))
            or (isinstance(l_high, float) and math.isnan(l_high))):
        return 1.0
    if l_low == l_high:
        l_low = max(0.1, l_low - 0.5)
        l_high = l_high + 0.5
    gamma = 1.0 - alpha
    num = std_max ** gamma - std_min ** gamma
    den = l_high ** gamma - l_low ** gamma
    cf = num / den
    if math.isinf(cf) or math.isnan(cf):
        return 1.0
    return cf


# ================================================================ [신규 함수]
# 생물 이용 가능 부피 보정 (Rmd#1 get_bio_correction 1:1)
# ================================================================
def get_bio_correction(alpha, std_min=1.0, std_max=5000.0, bio_max=100.0):
    """1–5000㎛ 전체 부피 중 1–100㎛(bio-available) 분율. p = 4 - alpha."""
    p = 4.0 - alpha
    return float((bio_max ** p - std_min ** p) / (std_max ** p - std_min ** p))


# ================================================================ [ToMEx 변환]
# Rmd#1 convert_and_rescale_standard 1:1 (질량 mg/L → 입자수/L + 크기·bio 보정)
# ================================================================
def convert_tomex(mass, size, dens, alpha=1.6):
    if mass is None or size is None or dens is None:
        return np.nan
    if any(isinstance(v, float) and math.isnan(v) for v in (mass, size, dens)) or mass <= 0:
        return np.nan
    radius_mm = (size / 2.0) / 1000.0
    vol_mm3 = (4.0 / 3.0) * math.pi * radius_mm ** 3
    part_mass_mg = dens * vol_mm3
    if part_mass_mg <= 0:
        return np.nan
    count_l = mass / part_mass_mg
    std_min, std_max, gamma = 1.0, 5000.0, 1.0 - alpha
    cf_size = (std_max ** gamma - std_min ** gamma) / ((size * 1.1) ** gamma - (size * 0.9) ** gamma)
    cf_bio = get_bio_correction(alpha)
    res = count_l * cf_size * cf_bio
    return res if (math.isfinite(res) and res > 0) else np.nan


# ================================================================ [AF 적용]
# Rmd#2 apply_af 1:1 (EC50/LC50 ×10, LOEC ×2, Acute ×10)
# ================================================================
def apply_af(val, endpoint, duration_type):
    if val is None or (isinstance(val, float) and (math.isnan(val) or val <= 0)):
        return np.nan
    af = 1.0
    if endpoint is not None and not (isinstance(endpoint, float) and math.isnan(endpoint)):
        ep = str(endpoint).strip().lower()
        if re.search(r"ec50|lc50", ep):
            af *= 10.0
        elif re.search(r"loec", ep):
            af *= 2.0
    if duration_type is not None and duration_type == "Acute":
        af *= 10.0
    return val / af


# ================================================================ [Koelmans 정렬]
# Rmd#2 align_koelmans 1:1 (섭취=food / 이동=trans 메커니즘별 표준 크기 정렬)
# ================================================================
def align_koelmans(conc, size_low, size_high, alpha, mech="food", mouth_size=100.0):
    if any(v is None or (isinstance(v, float) and math.isnan(v))
           for v in (conc, size_low, size_high)) or size_high <= 0 or conc <= 0:
        return np.nan
    if size_low == size_high:
        size_low = max(0.1, size_low * 0.9)
        size_high = size_high * 1.1
    gamma = 1.0 - alpha
    cf_count = (5000.0 ** gamma - 1.0 ** gamma) / (size_high ** gamma - size_low ** gamma)
    if mech == "food":
        p = 4.0 - alpha
        limit_upper = mouth_size
    else:
        p = 3.0 - alpha
        limit_upper = 83.0
    limit_upper = max(limit_upper, 1.1)
    cf_metric = (size_high ** p - size_low ** p) / (limit_upper ** p - 1.0 ** p)
    return float(conc * cf_count * cf_metric)


# ================================================================ [SSD fitting]
# ssdtools 호환: scipy MLE 피팅 + AICc 가중 모델 평균 HC 산출
# 분포: lnorm / llogis / gamma / weibull / gompertz / lgumbel (Rmd 사용 6종 전부)
# ================================================================
DIST_LIB = ("lnorm", "llogis", "gamma", "weibull", "gompertz", "lgumbel")

def _loglik_logspace(dist_name, logx):
    """log-농도에서 각 분포의 음의 로그우도(최소화 대상)."""
    if dist_name == "lnorm":
        mu, sd = logx.mean(), logx.std(ddof=0)
        return -np.sum(stats.norm.logpdf(logx, mu, max(sd, 1e-8))), (mu, max(sd, 1e-8))
    if dist_name == "llogis":
        m = logx.mean()
        sh = np.sqrt(3.0) * logx.std(ddof=0) / math.pi  # 분산=π²σ²/3 초기값
        sh = max(sh, 1e-3)
        return -np.sum(stats.fisk.logpdf(np.exp(logx), sh, scale=np.exp(m))), (sh, np.exp(m))
    if dist_name == "lgumbel":
        mu, be = logx.mean(), max(logx.std(ddof=0) * math.sqrt(6) / math.pi, 1e-8)
        return -np.sum(stats.gumbel_r.logpdf(logx, mu, be)), (mu, be)
    return None, None

def _cdf(dist_name, x, params):
    """분포 CDF. params는 피팅에서 반환한 형태."""
    if dist_name == "lnorm":
        mu, sd = params
        return stats.norm.cdf(np.log(x), mu, sd)
    if dist_name == "llogis":
        sh, sc = params
        return stats.fisk.cdf(x, sh, scale=sc)
    if dist_name == "gamma":
        a, _, sc = params
        return stats.gamma.cdf(x, a, loc=0, scale=sc)
    if dist_name == "weibull":
        a, _, sc = params
        return stats.weibull_min.cdf(x, a, loc=0, scale=sc)
    if dist_name == "gompertz":
        a, _, sc = params
        return stats.gompertz.cdf(x, a, loc=0, scale=sc)
    if dist_name == "lgumbel":
        mu, be = params
        return stats.gumbel_r.cdf(np.log(x), mu, be)
    raise ValueError(dist_name)

def _fit_one(dist_name, x):
    """단일 분포 MLE 피팅 → (nll, k, params). log-공간 절친값 보호."""
    lx = np.log(x)
    if dist_name in ("lnorm", "llogis", "lgumbel"):
        nll, params = _loglik_logspace(dist_name, lx)
        return nll, params
    if dist_name == "gamma":
        a, loc, sc = stats.gamma.fit(x, floc=0)
        return -np.sum(stats.gamma.logpdf(x, a, loc=0, scale=sc)), (a, 0.0, sc)
    if dist_name == "weibull":
        a, loc, sc = stats.weibull_min.fit(x, floc=0)
        return -np.sum(stats.weibull_min.logpdf(x, a, loc=0, scale=sc)), (a, 0.0, sc)
    if dist_name == "gompertz":
        try:
            a, loc, sc = stats.gompertz.fit(x, floc=0)
        except Exception:
            a = 1.0
            sc = x.mean()
        return -np.sum(stats.gompertz.logpdf(x, a, loc=0, scale=sc)), (a, 0.0, sc)
    raise ValueError(dist_name)

def ssd_fit_dists(conc, dists=("lnorm", "llogis", "gamma")):
    """종 대표 농도 배열에 다중 분포 MLE 피팅. R ssd_fit_dists 대응.
    returns {name: {'nll', 'k', 'params', 'n', '_x'}}"""
    x = np.asarray([c for c in conc if np.isfinite(c) and c > 0], float)
    fits = {}
    for d in dists:
        try:
            nll, params = _fit_one(d, x)
            if not np.isfinite(nll):
                continue
            fits[d] = {"nll": float(nll), "k": len(params),
                       "params": params, "n": len(x), "_x": x}
        except Exception:
            continue
    return fits

def ssd_hc_modelavg(fits, percent=5.0, nboot=0, seed=42, ci=False):
    """AICc 가중 모델 평균(mixture CDF) HC 산출. R ssd_hc(average=TRUE) 대응.
    returns dict {est, lcl?, ucl?} — 단위: particles/L"""
    if not fits:
        return {"est": np.nan}
    n = next(iter(fits.values()))["n"]

    def weights(fdict):
        aic = {d: 2 * f["nll"] + 2 * f["k"] for d, f in fdict.items()}
        mn = min(aic.values())
        w = {d: math.exp(-0.5 * (aic[d] - mn)) for d in aic}
        tw = sum(w.values())
        return {d: v / tw for d, v in w.items()}

    w = weights(fits)

    def mixture_cdf(xv):
        return sum(w[d] * _cdf(d, xv, fits[d]["params"]) for d in w)

    lo, hi = _bounds(fits)
    est = _solve_hc(mixture_cdf, percent / 100.0, lo, hi)
    out = {"est": est}
    if ci and nboot > 0:
        all_x = next(iter(fits.values()))["_x"]
        rng = np.random.default_rng(seed)
        boots = []
        for _ in range(nboot):
            bx = rng.choice(all_x, size=len(all_x), replace=True)
            bfit = ssd_fit_dists(bx, dists=tuple(fits.keys()))
            if not bfit:
                continue
            bw = weights(bfit)

            def bmix(xv, bw=bw, bfit=bfit):
                return sum(bw[d] * _cdf(d, xv, bfit[d]["params"]) for d in bw)

            try:
                boots.append(_solve_hc(bmix, percent / 100.0, lo, hi))
            except Exception:
                continue
        if len(boots) >= max(10, nboot // 5):
            out["lcl"] = float(np.percentile(boots, 2.5))
            out["ucl"] = float(np.percentile(boots, 97.5))
    return out

def _bounds(fits):
    xs = None
    for f in fits.values():
        if "_x" in f:
            xs = f["_x"]
            break
    if xs is not None:
        return max(1e-6, xs.min() / 100), xs.max() * 100
    return 1e-4, 1e8

def _solve_hc(cdf_func, p, lo, hi):
    hh = hi
    while cdf_func(hh) < p and hh < 1e12:
        hh *= 10
    r = optimize.brentq(lambda t: cdf_func(t) - p, lo, hh, maxiter=200)
    return float(r)

def model_weights(fits):
    """모델별 AICc 가중치 표 (모델 비교/보고용)."""
    if not fits:
        return pd.DataFrame()
    n = next(iter(fits.values()))["n"]
    rows = []
    for d, f in fits.items():
        aic = 2 * f["nll"] + 2 * f["k"]
        aicc = aic + (2 * f["k"] * (f["k"] + 1)) / max(n - f["k"] - 1, 1)
        rows.append({"Model": d, "AICc": aicc, "n": n})
    df = pd.DataFrame(rows)
    mn = df["AICc"].min()
    df["Weight"] = np.exp(-0.5 * (df["AICc"] - mn))
    df["Weight"] = df["Weight"] / df["Weight"].sum()
    return df.sort_values("AICc").reset_index(drop=True)

def bfit_ok(x, fits):
    return bool(fits)

# ================================================================ [유틸]
def is_apical(endpoint_classification, extra=None):
    """SOD/CAT 등 in-vitro·바이오마커 제외 판정 (R 필터 'sod|cat|...' 1:1)."""
    if endpoint_classification is None:
        return True
    pat = EP_REGEX if extra is None else EP_REGEX + "|" + extra
    return not bool(re.search(pat, str(endpoint_classification), re.I))


def media_group(environment, route=None):
    """Seawater / Freshwater / Sediment 분류 (R Media_Group 로직 1:1)."""
    e = str(environment or "").lower()
    r = str(route or "").lower()
    if "sediment" in e or "sediment" in r:
        return "Sediment"
    if "fresh" in e:
        return "Freshwater"
    return "Seawater"


def media_alpha(mg):
    """퇴적물 2.57 / 해수 2.07 (R Alpha 로직 1:1)."""
    return 2.57 if mg == "Sediment" else 2.07


def load_lab_final():
    """최종 취합본 1단계(MP) + 2단계(Paint) → 긴 형식 DataFrame.
    신규 레이아웃: 엔드포인트가 (㎎/L, ㎎/kg) 쌍열 — 매체별 단위 열 선택."""
    frames = []
    for sheet, stage in [("1단계 (MP)", "MP"), ("2단계 (Paint Particles)", "Paint")]:
        df = pd.read_excel(FINAL_XLSX, sheet_name=sheet, skiprows=3, header=None)
        df = df.iloc[:, :28]
        df.columns = (["No", "Category", "Species_Name", "Life_stage", "Environment",
                       "Type", "Shape", "Size_um", "Test_duration", "Target_Cell",
                       "Acute", "Chronic", "Endpoint_classification",
                       "EC50_L", "EC50_k", "EC10_L", "EC10_k", "LC50_L", "LC50_k",
                       "LC10_L", "LC10_k", "NOEC_L", "NOEC_k", "LOEC_L", "LOEC_k",
                       "MATC_L", "MATC_k", "Etc"][: df.shape[1]])
        df["Stage"] = stage
        frames.append(df)
    raw = pd.concat(frames, ignore_index=True)

    rows = []
    for _, r in raw.iterrows():
        for ep, cl_name, ck_name in ENDPOINT_PAIRS:
            for suffix, col_name in (("_L", cl_name), ("_k", ck_name)):
                v = clean_tox_value(r[col_name])
                if np.isfinite(v):
                    rows.append({
                        "No": r["No"], "Category": r["Category"],
                        "Species_Name": str(r["Species_Name"]).strip() if pd.notna(r["Species_Name"]) else None,
                        "Life_stage": r["Life_stage"], "Environment": r["Environment"],
                        "Type": r["Type"], "Shape": r["Shape"], "Size_um": r["Size_um"],
                        "Test_duration": r["Test_duration"], "Stage": r["Stage"],
                        "Acute": r["Acute"], "Chronic": r["Chronic"],
                        "Endpoint_classification": r["Endpoint_classification"],
                        "Endpoint_Type": ep, "Endpoint_Unit": "mg/L" if suffix == "_L" else "mg/kg",
                        "Endpoint_Value": v,
                    })
    long = pd.DataFrame(rows)
    long["Species_Name"] = long["Species_Name"].fillna("")
    long["Acute"] = long["Acute"].astype(str).str.strip().str.lower()
    long["Chronic"] = long["Chronic"].astype(str).str.strip().str.lower()
    long["Media_Group"] = long["Environment"].map(media_group)
    long["Alpha"] = long["Media_Group"].map(media_alpha)
    si = long["Size_um"].map(parse_size_info)
    long["Size_Low"] = [s[0] for s in si]
    long["Size_High"] = [s[1] for s in si]
    long["Size_um_Midpoint"] = [s[2] for s in si]
    long["CF_size"] = [get_powerlaw_cf(l, h, a) for l, h, a in
                       zip(long["Size_Low"], long["Size_High"], long["Alpha"])]
    long["CF_bio"] = long["Alpha"].map(lambda a: get_bio_correction(a))
    long["Rescaled_Endpoint_Value"] = long["Endpoint_Value"] * long["CF_size"]
    long["학명"] = long["Species_Name"].str.strip()
    long["Origin"] = "Domestic (국내종)"  # 본 연구 실험 데이터 = 국내종 (R 원본 동일)
    return long