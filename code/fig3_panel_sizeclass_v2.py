#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig3_panel_sizeclass_v2.py — Fig.3 4패널 연구계열 단위 크기구성 재판정 (v2).

[기준: 승인된 '연구계열 단위 분류 v2'(2026-10-14 INFO verbatim)]
- 동일한 생물종에 대해 비교 가능한 시험조건에서 두 가지 이상의 독립적인 입자
  크기 독성자료가 확보되었다면, 시험 날짜가 달라도 Multiple-size data로 분류하되
  Within-run과 Across-run을 구분하여 기록한다.
- 분류 단위 = 종 × SSD 패널(FD-Q1/TT-Q1/FD-Median/TT-Median).
- 판정 근거 = 해당 패널 실제 투입 적격 원자료 레코드(엔진 load_and_screen → align_one 집합).
  패널 적격 크기조건 1개면 Single-size(=data 명칭의 이유; 원논문 전 실험조건 아님).
- 연구계열(study series) 정의:
  * ours(DOM*) : 동일 연구실 전체가 한 연구계열(species-level series) — 날짜 무관 통합.
  * 문헌       : DOI 1차 토큰(study) = 연구계열. 독립 연구를 임의로 묶지 않는다.
- 4클래스: Multiple-size data(연구계열 ≥2 크기 조건만) / Single-size data(전 연구계열
  1 크기 조건) / Both size-data types(다중크기 연구계열 + 별도 단일크기 연구계열 병존) /
  Unclassified(적격 레코드의 크기 조건 자체 확인 불가 NaN).
- run_type(보충자료 기록열): ours 연구계열 기준
  Within-run = 동일 DOM{date} 런 안에 ≥2 크기 조건이 존재
  Across-run = 다른 날짜 런들에 걸쳐 ≥2 크기 조건 (동일 날짜 다중크기 없음)
  Within+Across-run = 둘 다 (동일 날짜 다중 + 날짜 간 추가 크기)
  single = 적격 크기조건 1개뿐. 날짜는 분류 기준이 아니라 run 기록용 보조 변수.
- 불변: Conc·n·TR·HC5 일절 재계산 금지(재현 게이트 1e-6로 배전 입력 동일성 확인).
[계보] v1(fig3_panel_sizeclass_v1.py, DOM{date} 캠페인 단위)은 이전 계보로 보존 —
이번 v1 재생산은 배전 재현 게이트 재확인용이며 계산상 결함이 아니라 기준 변경이다.
[산출] 02_data/py_out/fig3_sizeclass_v2.csv(+_summary.csv), 콘솔 v1→v2 diff.
"""
from __future__ import annotations
import importlib.util, json, os
import numpy as np
import pandas as pd

from pathlib import Path
DELIV = str(Path(__file__).resolve().parents[1] / 'data')
ENGINE = str(Path(__file__).resolve().parents[1] / 'code' / 'ssd_tomex_alpha207.py')
MERGED = os.path.join(DELIV, 'merged_T12_20_300.csv')
OUT = str(Path(__file__).resolve().parents[1] / 'output' / 'fig3_sizeclass')
OUT_V1 = str(Path(__file__).resolve().parents[1] / 'output' / 'fig3_sizeclass_v1check')
SIZE_MIN, SIZE_MAX, ALPHA = 20.0, 300.0, 2.07

PANELS = {
    ('food_dilution', 'q1_all'): 'FD-Q1',
    ('translocation', 'q1_all'): 'TT-Q1',
    ('food_dilution', 'median_organism_population'): 'FD-Median',
    ('translocation', 'median_organism_population'): 'TT-Median',
}


def series_of(doi: str) -> str:
    """연구계열 ID: ours=OURS(기관 통합), 문헌=DOI 1차 토큰(study)."""
    s = str(doi)
    return 'OURS' if s.startswith('DOM') else s.split('_')[0]


def srun_type(group: pd.DataFrame) -> str:
    """ours 연구계열의 run 기록: 동일 날짜 다중크기 존재 여부 + 날짜 수."""
    dts = sorted({str(d).split('_')[0][3:] for d in group['DOI']})  # DOM{yyyymmdd}
    per_date = group.assign(date=[str(d).split('_')[0][3:] for d in group['DOI']]) \
                    .groupby('date')['len_um'].nunique()
    n_sizes = int(group['len_um'].nunique())
    if n_sizes < 2:
        return 'single'
    if (per_date >= 2).any() and len(dts) >= 2:
        return 'Within+Across-run'
    if (per_date >= 2).any():
        return 'Within-run'
    if len(dts) >= 2:
        return 'Across-run'
    return 'single'  # 방어: 단일 날짜 1크기


def classify(group: pd.DataFrame) -> tuple[str, str]:
    """연구계열별 적격 크기조건 수로 4클래스 판정 + ours run_type 문자열."""
    g = group.assign(series=group['DOI'].map(series_of))
    size_info = g[g['len_um'].notna()]
    n_unclass = int(g['len_um'].isna().sum())
    per = size_info.groupby('series')['len_um'].nunique()
    types = {'multiple' if v >= 2 else 'single' for v in per}
    if 'multiple' in types and 'single' in types:
        cls = 'Both size-data types'
    elif 'multiple' in types:
        cls = 'Multiple-size data'
    else:
        cls = 'Single-size data'
    if n_unclass > 0:
        cls = 'Unclassified'
    if group['is_dom'].all():
        ours = group
    else:
        ours = group[group['is_dom']]
    run = srun_type(ours) if len(ours) else ''
    return cls, run


def camp_detail(group: pd.DataFrame) -> str:
    """근거 문자열: 연구계열|ours|sizes(...)|dates(...) 형식."""
    g = group.assign(series=group['DOI'].map(series_of))
    parts = []
    for s, sub in g.groupby('series', sort=True):
        ours = bool(sub['is_dom'].all())
        sizes = ','.join('%g' % v for v in sorted(set(np.round(sub['len_um'].dropna(), 4))))
        if ours:
            dts = sorted({str(d).split('_')[0][3:] for d in sub['DOI']})
            parts.append('%s|ours|sizes=%s|dates=%s' % (s, sizes, '+'.join(dts)))
        else:
            parts.append('%s|ToMEx|sizes=%s' % (s, sizes))
    return '; '.join(parts)


def main() -> None:
    spec = importlib.util.spec_from_file_location('ssd_engine207', ENGINE)
    eng = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(eng)

    raw, screen, audit = eng.load_and_screen(MERGED, None, 3)
    print('[1] screening: raw=%d screened=%d species=%d studies=%d' % (
        audit['raw_rows'], audit['screened_rows'], audit['species'], audit['studies']))
    env = screen[screen['Environment'] == 'Marine']

    aligned: dict[str, pd.DataFrame] = {}
    for mech in ('food_dilution', 'translocation'):
        recs = []
        for _, row in env.iterrows():
            out = eng.align_one(row, mech, SIZE_MIN, SIZE_MAX, ALPHA)
            if out:
                recs.append({
                    'DOI': row['DOI'], 'is_dom': str(row['DOI']).startswith('DOM'),
                    'Species': row['Species'],
                    'LBO': row['Level of Biological Organization'],
                    'len_um': float(row['_length_um']),
                    'aligned_particles_L': out['aligned_particles_L'],
                })
        aligned[mech] = pd.DataFrame(recs)
        print('[2] %s aligned rows=%d NaN-size=%d' % (
            mech, len(aligned[mech]), int(aligned[mech]['len_um'].isna().sum())))

    # ---- 재현 게이트(불변 확인): 종별 Conc == 배전 ssdtools_input ----
    gate = {}
    for mech in ('food_dilution', 'translocation'):
        dfa = aligned[mech]
        for collapse, mode in (('q1_all', 0.25), ('median_organism_population', 0.50)):
            inc = dfa if collapse == 'q1_all' else dfa[dfa['LBO'].str.lower().isin(
                ['organism', 'population'])]
            g = inc.groupby('Species').agg(
                Conc=('aligned_particles_L', lambda x: float(np.quantile(x, mode))),
                n_effects=('aligned_particles_L', 'size')).reset_index()
            f = os.path.join(DELIV, 'ssdtools_input_marine_%s_%s.csv' % (mech, collapse))
            d = pd.read_csv(f)
            mm = d.merge(g, on='Species', suffixes=('_deliv', '_recomp'))
            rel = (mm['Conc_deliv'] - mm['Conc_recomp']).abs() / mm['Conc_deliv'].abs()
            pid = PANELS[(mech, collapse)]
            gate[pid] = {'delivered': int(len(d)), 'recomp': int(len(g)),
                         'matched': int(len(mm)),
                         'max_rel_err': float(rel.max()) if len(rel) else None}
    print('[3] reproduction gate:', json.dumps(gate))
    ok = all(v['matched'] == v['delivered'] == v['recomp']
             and v['max_rel_err'] is not None and v['max_rel_err'] <= 1e-6
             for v in gate.values())
    print('[3] gate PASS=%s' % ok)
    if not ok:
        raise SystemExit('reproduction gate FAILED — 배전 입력과 불일치, 분류 중단')

    # ---- v2 분류(연구계열 단위) ----
    rows = []
    for (mech, collapse), pname in PANELS.items():
        dfa = aligned[mech]
        inc = dfa if collapse == 'q1_all' else dfa[dfa['LBO'].str.lower().isin(
            ['organism', 'population'])]
        for sp, sub in inc.groupby('Species'):
            g = sub.assign(series=sub['DOI'].map(series_of))
            size_info = g[g['len_um'].notna()]
            cls, run = classify(sub)
            srcs = 'ours' if sub['is_dom'].all() else (
                'ToMEx' if sub['is_dom'].sum() == 0 else 'ours+ToMEx')
            rows.append({
                'panel': pname, 'Species': sp, 'source_mix': srcs,
                'class_v2': cls, 'run_type_ours': run,
                'n_effects': int(len(sub)),
                'n_series': int(g['series'].nunique()),
                'n_series_multi': int((size_info.groupby('series')['len_um']
                                       .nunique() >= 2).sum()),
                'n_series_single': int((size_info.groupby('series')['len_um']
                                        .nunique() == 1).sum()),
                'size_union': ','.join('%g' % v for v in sorted(set(
                    np.round(g['len_um'].dropna(), 4)))),
                'series_detail': camp_detail(sub),
            })
    cls_df = pd.DataFrame(rows)

    # ---- v1 재생산본과 대조(diff) ----
    v1 = pd.read_csv(os.path.join(OUT_V1, 'fig3_panel_sizeclass_v1.csv'),
                     encoding='utf-8-sig')
    m = cls_df.merge(v1[['panel', 'Species', 'class']],
                     on=['panel', 'Species'], how='left')
    m['class_v1'] = m['class']
    diff = m[m['class'].str.replace(' ', '') != m['class'].str.replace(' ', '')]
    diff = m[m['class_v2'] != m['class_v1']]
    cls_df_out = cls_df.copy()
    cls_df_out['class_v1'] = m['class_v1'].values
    cls_df_out['changed'] = (cls_df_out['class_v2'] != cls_df_out['class_v1'])

    # ---- 패널 집계 ----
    summ = cls_df.groupby('panel').apply(
        lambda d: pd.Series({'n_species': len(d),
                             **d['class_v2'].value_counts().to_dict()}),
        include_groups=False).reset_index()
    summ['both_species'] = [
        '; '.join(sorted(cls_df.loc[(cls_df['panel'] == p) &
                                    (cls_df['class_v2'] == 'Both size-data types'),
                                    'Species']))
        for p in summ['panel']]

    cls_df_out.to_csv(os.path.join(OUT, 'fig3_sizeclass_v2.csv'),
                      index=False, encoding='utf-8-sig')
    summ.to_csv(os.path.join(OUT, 'fig3_sizeclass_v2_summary.csv'),
                index=False, encoding='utf-8-sig')
    print('[4] v2 summary:\n%s' % summ.to_string(index=False))
    ch = cls_df_out[cls_df_out['changed']]
    print('[5] v1->v2 changed species=%d:' % len(ch))
    if len(ch):
        print(ch[['panel', 'Species', 'source_mix', 'class_v1', 'class_v2',
                  'run_type_ours']].to_string(index=False))
    print('[6] ours rows(v2 판정근거):\n%s' % cls_df_out[
        cls_df_out['source_mix'].str.contains('ours')][
        ['panel', 'Species', 'class_v2', 'run_type_ours', 'size_union',
         'series_detail']].to_string(index=False))
    print('[7] files: fig3_sizeclass_v2.csv(%d)/_summary.csv -> %s' % (
        len(cls_df_out), OUT))


if __name__ == '__main__':
    main()