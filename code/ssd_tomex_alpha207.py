#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ToMEx -> exploratory SSD -> TR1, TR2, TR3, TR4.

Mehinto et al. (2022) and their 2025 erratum: https://doi.org/10.1186/s43591-025-00149-2
Hampton et al. (2025): https://doi.org/10.1186/s43591-025-00145-6

IMPORTANT: NO claim to replicate the published SSD numerical results.
- The input CSV has no paper-specific 14-item acceptability decisions nor already aligned doses.
- The environmental joint size/shape distribution and corrected supplementary raw
  SSD dataset are missing. The spherical environmental geometry and analytical
  power-law integration below are *explicit proxies*, NOT the full Koelmans model.
- For publication-grade inference use the generated 'ssdtools_input_*.csv' with
  run_ssdtools.R after expert screening/actual alignment has been validated.
- The native Python exploratory model is a *single lognormal* (not the multi-
  distribution model averaging used by R ssdtools).

Usage:
  python ssd_tomex.py --input ToMEx_Exploration2026-03-24.csv --output ssd_results
  python ssd_tomex.py --input ToMEx_Exploration2026-03-24.csv --output ssd_20_300 --size-min 20 --size-max 300
  python ssd_tomex.py --help
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

EXCLUDED_DOI = {'10.1007/s00128-017-2186-3'}  # Mehinto 2025 correction
SUPPORTED = {'NOEC': 1.0, 'LOEC': 2.0, 'EC50': 10.0, 'LC50': 10.0}
COLS = {
    'dose': 'Unaligned Dose Values',
    'length': 'Particle Length (μm)',
    'width': 'Particle Width (μm)',
    'mass': 'Particle Mass (mg)',
    'rho': 'Density (g/cm^3)',
    'body': 'Estimated Body Length (cm)',
    'mouth': 'Estimated Maximum Ingestible Size (mm)',
}


def number(x):
    try:
        n = float(str(x).strip().replace(',', ''))
        return n if np.isfinite(n) else np.nan
    except (TypeError, ValueError):
        return np.nan


def ellipsoid_metrics_um(length, width, shape='Sphere'):
    """Updated in 2025 erratum: SEMI-axes in volume AND surface area.
    The original ToMEx table lacks height: for non-spheres c=0.67*b
    is an explicit approximation, with fiber cross-section equal to width.
    """
    L, W = float(length), float(width)
    if L <= 0 or W <= 0:
        return (np.nan, np.nan)
    a, b = L / 2.0, W / 2.0
    if str(shape).lower() == 'sphere':
        c = a
        b = a
    else:
        c = 0.67 * b
    volume = (4 * np.pi / 3) * a * b * c
    power = 1.6
    area = 4 * np.pi * (((a*b)**power + (a*c)**power + (b*c)**power) / 3)**(1/power)
    return volume, area


def integrate_power(lo, hi, power):
    """Integral of x^power over [lo,hi], stable at power=-1."""
    if hi <= lo:
        return 0.0
    if abs(power + 1) < 1e-11:
        return np.log(hi / lo)
    return (hi**(power + 1) - lo**(power + 1)) / (power + 1)


def ambient_moment_times_availability(k, avail_upper_um, lower, upper, alpha):
    """E[L^k * I(L<=avail_upper)] under p(L) proportional to L^-alpha.
    Avoids numerical underflow from independently calculating f_avail * mean.
    Environmental geometry is equivalent-sphere, documented as proxy.
    """
    hi = min(upper, avail_upper_um)
    return integrate_power(lower, hi, k - alpha) / integrate_power(lower, upper, -alpha)


def corrected_mouth_um(body_cm, fallback_mm):
    # Mehinto (2025) erratum: body_cm * 10 INSIDE log10.
    if np.isfinite(body_cm) and body_cm > 0:
        return 1000 * 10**(0.9341 * np.log10(body_cm * 10) - 1.1200)
    if np.isfinite(fallback_mm) and fallback_mm > 0:
        return 1000 * fallback_mm
    return np.nan


def load_and_screen(p, quality_csv=None, min_doses=3):
    df = pd.read_csv(p, encoding='utf-8-sig', dtype=str, keep_default_na=False)
    df['source_csv_row'] = np.arange(len(df)) + 2
    audit = {'raw_rows': len(df), 'reason_counts': {}}
    flags = pd.Series(True, index=df.index)
    def criterion(reason, selection):
        nonlocal flags
        selection = selection.fillna(False).astype(bool)
        audit['reason_counts'][reason] = int((flags & ~selection).sum())
        flags = flags & selection

    df['_dose_mg_L'] = df[COLS['dose']].map(number) # 1 ug/mL == 1 mg/L
    df['_length_um'] = df[COLS['length']].map(number)
    df['_width_um'] = df[COLS['width']].map(number)
    df['_rho_g_cm3'] = df[COLS['rho']].map(number)
    df['_body_cm'] = df[COLS['body']].map(number)
    df['_mouth_mm'] = df[COLS['mouth']].map(number)
    df['_number_doses'] = df['Number of Doses'].map(number)
    metric = df['Effect Metric'].astype(str).str.upper().str.strip()
    df['_metric'] = metric
    criterion('not_noec_or_loec', metric.isin(SUPPORTED.keys()))
    criterion('honec_excluded_and_inconsistent_effect',
              ((metric == 'NOEC') & (df['Effect'].str.lower() == 'no')) |
              ((metric != 'NOEC') & (df['Effect'].str.lower() == 'yes')))
    criterion('not_water_particle_only',
              (df['Experiment Type'].str.lower() == 'particle only') &
              (df['Exposure Route'].str.lower() == 'water'))
    criterion('negative_control_not_confirmed',
              df['Negative Control'].str.lower().isin(['y','yes','true']))
    criterion('insufficient_dose_levels', df['_number_doses'] >= min_doses)
    criterion('invalid_effect_dose', df['_dose_mg_L'] > 0)
    criterion('lab_particle_outside_1_5000_um',
              df['_length_um'].between(1,5000, inclusive='both'))
    criterion('invalid_particle_density', df['_rho_g_cm3'] > 0)
    criterion('size_or_geometry_unresolved',
              ((df['_width_um'] > 0) | (df['Shape'].str.lower() == 'sphere')))
    criterion('erratum_excluded_study', ~df['DOI'].str.lower().isin(EXCLUDED_DOI))

    # Optional expert-reviewed whitelist. Row-level accept lists more precise
    # than DOI-level; select rows by (source_csv_row, approved = 1).
    if quality_csv:
        quality = pd.read_csv(quality_csv)
        required = {'source_csv_row','approved'}
        if not required.issubset(quality.columns):
            raise ValueError(f'quality CSV must contain {required}')
        approved = quality.loc[quality['approved'].astype(str).str.lower().isin(['1','true','yes','y']), 'source_csv_row']
        criterion('not_expert_approved', df['source_csv_row'].isin(approved))
        audit['expert_screening_applied'] = True
    else:
        audit['expert_screening_applied'] = False
    subset = df.loc[flags].copy()
    # Remove only exact analytic duplicates, preserving distinct endpoint/size/method rows.
    key = ['DOI','Species','Endpoint','Life Stage','Acute/Chronic',
           'Exposure Duration (days)','_metric','_dose_mg_L','_length_um',
           '_width_um','Shape','Polymer','Level of Biological Organization']
    before = len(subset)
    subset = subset.drop_duplicates(subset=key, keep='first')
    audit['exact_duplicate_rows_removed'] = before - len(subset)
    audit['screened_rows'] = len(subset)
    audit['species'] = int(subset['Species'].nunique())
    audit['studies'] = int(subset['DOI'].nunique())
    return df, subset, audit


def align_one(row, mechanism, lower, upper, alpha):
    """Effect concentration conversion by preservation of lab particle volume
    or surface area in accessible environmental particles. Source CSV dose is
    mass concentration (mg/L), so FIRST convert to lab particles/L.
    """
    width = row['_width_um'] if np.isfinite(row['_width_um']) else row['_length_um']
    vlab, alab = ellipsoid_metrics_um(row['_length_um'], width, row['Shape'])
    mass_lab_mg = vlab * row['_rho_g_cm3'] * 1e-9
    if not np.isfinite(mass_lab_mg) or mass_lab_mg <= 0:
        return None
    n_lab = row['_dose_mg_L'] / mass_lab_mg
    metric = row['_metric']
    AF = SUPPORTED[metric] * (10.0 if row['Acute/Chronic'].lower()=='acute' else 1.0)
    n_lab_chronic = n_lab / AF

    if mechanism == 'food_dilution':
        if row['Organism Group'].lower() in ('algae','plant','bacterium'):
            return None
        available_size = corrected_mouth_um(row['_body_cm'],row['_mouth_mm'])
        k, shape_factor, lab_metric = 3, np.pi/6, vlab
    elif mechanism == 'translocation':
        available_size = 83.0
        # Lab particle >83 um is not taken to imply translocation; exclude.
        if row['_length_um'] > 83:
            return None
        k, shape_factor, lab_metric = 2, np.pi, alab
    else:
        raise ValueError(mechanism)

    if not np.isfinite(available_size):
        return None
    avail_moment = ambient_moment_times_availability(k, available_size, lower, upper, alpha)
    available_metric = shape_factor * avail_moment
    if not np.isfinite(available_metric) or available_metric <= 0:
        return None
    # number of all ambient particles in defined size range [lower,upper]
    n_env = n_lab_chronic * lab_metric / available_metric
    if not (np.isfinite(n_env) and n_env > 0):
        return None
    # Reference environmental mass of an equivalent sphere, not lab mass.
    mean_env_mass = (np.pi/6) * ambient_moment_times_availability(3,upper,lower,upper,alpha) * 1.07e-9
    return {
        'lab_mass_conc_mg_L': row['_dose_mg_L'],
        'lab_particle_mass_mg_corrected': mass_lab_mg,
        'lab_particles_L': n_lab,
        'descriptor_AF': SUPPORTED[metric],
        'acute_to_chronic_AF': 10.0 if row['Acute/Chronic'].lower()=='acute' else 1.0,
        'chronic_lab_particles_L': n_lab_chronic,
        'available_upper_um': available_size,
        'ambient_accessible_metric_per_total_particle': available_metric,
        'aligned_particles_L': n_env,
        'aligned_mg_L_ambient_proxy': n_env*mean_env_mass,
        'ambient_mean_particle_mass_mg_proxy': mean_env_mass,
    }


def fit_lognormal_hc(sample, p):
    """Python exploratory parametric lognormal SSD (not R ssdtools)."""
    z = np.log(np.asarray(sample, dtype=float))
    if len(z)<5 or not np.all(np.isfinite(z)):
        return np.nan
    sd = float(np.std(z,ddof=0))
    return float(np.exp(np.mean(z) + norm.ppf(p)*sd))


def bootstrap_hc(species_values, p_values=(0.05,0.10), nboot=1000, seed=20261011):
    values=np.asarray(species_values,dtype=float)
    if len(values)<5:
        return {}
    rng = np.random.default_rng(seed)
    mu = np.log(values).mean()
    sigma = np.log(values).std(ddof=0)
    out={}
    for p in p_values:
        # Parametric bootstrap of fitted lognormal, analogous in spirit (not
        # numerically identical) to ssdtools confidence interval estimation.
        x=np.exp(rng.normal(mu,sigma,size=(nboot,len(values))))
        mus=np.log(x).mean(axis=1)
        sigmas=np.log(x).std(axis=1,ddof=0)
        boots=np.exp(mus+norm.ppf(p)*sigmas)
        out[str(p)] = {'est':fit_lognormal_hc(values,p),
                       'lcl':float(np.quantile(boots,.025)),
                       'ucl':float(np.quantile(boots,.975))}
    return out


def plot_ssd(samples, title, output):
    fig,ax = plt.subplots(figsize=(8,5.4),dpi=160)
    xmin,xmax=np.inf,0
    for collapse, df in samples.items():
        vals=np.sort(df['Conc'].to_numpy(float))
        if len(vals)<5: continue
        xs=np.geomspace(vals.min()/2,vals.max()*2,250)
        zs=np.log(vals); mu=zs.mean(); sd=zs.std(ddof=0)
        ys=norm.cdf((np.log(xs)-mu)/sd)
        ax.plot(xs,ys,label=f'{collapse}: fitted lognormal (n={len(vals)})')
        ax.scatter(vals, np.arange(1,len(vals)+1)/(len(vals)+1),s=14,alpha=.6)
        xmin=min(xmin,vals.min());xmax=max(xmax,vals.max())
    if np.isfinite(xmin):
        ax.set_xscale('log');ax.set_xlim(xmin/4,xmax*4)
    ax.set_ylim(0,1)
    ax.set_xlabel('Aligned microplastic concentration (particles/L) [proxy]')
    ax.set_ylabel('Cumulative fraction of species')
    ax.set_title(title)
    ax.grid(True,alpha=.15)
    ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(output,bbox_inches='tight');plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True,type=Path)
    ap.add_argument('--output',default='ssd_results',type=Path)
    ap.add_argument('--quality-csv',type=Path,help='Expert row whitelist source_csv_row,approved')
    ap.add_argument('--size-min',type=float,default=1)
    ap.add_argument('--size-max',type=float,default=5000)
    ap.add_argument('--nboot',type=int,default=1000)
    ap.add_argument('--min-doses',type=int,default=3)
    args=ap.parse_args()
    if not (args.size_min>=1 and args.size_max>args.size_min and args.size_max<=5000):
        ap.error('size limits must satisfy 1 <= min < max <= 5000')
    args.output.mkdir(parents=True,exist_ok=True)
    raw,screen,audit=load_and_screen(args.input,args.quality_csv,args.min_doses)
    audit['ambient_size_um']=[args.size_min,args.size_max]
    audit['ambient_geometry']='equivalent sphere (ASSUMPTION)'
    audit['environment_alpha']={'Marine':2.07,'Freshwater':2.07}
    audit['model']='single lognormal Python proxy; R SSDTools output is authoritative only after inputs are validated'
    audit['excluded_doi_erratum']=list(EXCLUDED_DOI)
    (args.output/'screening_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    screen.to_csv(args.output/'screened_candidate_rows.csv',index=False,encoding='utf-8-sig')
    quality_cols = ['source_csv_row','DOI','Species','Organism Group','Endpoint',
                    'Level of Biological Organization','Effect Metric','Negative Control',
                    'Number of Doses','Size Validated?','Polymer Validated?',
                    'Concentration Validated?','Particle Length (μm)',
                    'Particle Width (μm)','Unaligned Dose Values']
    review = screen[quality_cols].copy()
    review['approved'] = ''
    review['review_notes'] = ''
    review.to_csv(args.output/'expert_quality_review_template.csv',index=False,encoding='utf-8-sig')

    results=[]
    for environment,alpha in [('Marine',2.07),('Freshwater',2.07)]:
        env=screen[screen['Environment']==environment]
        for mechanism in ['food_dilution','translocation']:
            aligned=[]
            for _,row in env.iterrows():
                out=align_one(row, mechanism,args.size_min,args.size_max,alpha)
                if out:
                    aligned.append({
                        'source_csv_row':int(row['source_csv_row']),'DOI':row['DOI'],
                        'Species':row['Species'],'Organism Group':row['Organism Group'],
                        'Level of Biological Organization':row['Level of Biological Organization'],
                        'Endpoint':row['Endpoint'],'Effect Metric':row['Effect Metric'],
                        'Acute/Chronic':row['Acute/Chronic'], 'Shape':row['Shape'],
                        'Particle Length (μm)':row['Particle Length (μm)'], **out
                    })
            if not aligned:
                warnings.warn(f'{environment}/{mechanism}: no aligned rows')
                continue
            dfa=pd.DataFrame(aligned)
            prefix=f'{environment.lower()}_{mechanism}'
            dfa.to_csv(args.output/f'aligned_{prefix}.csv',index=False,encoding='utf-8-sig')
            datasets={}
            for collapse,mode in [('q1_all',.25),('median_organism_population',.50)]:
                inc=dfa if collapse=='q1_all' else dfa[dfa['Level of Biological Organization'].str.lower().isin(['organism','population'])]
                by_species=inc.groupby('Species',as_index=False).agg(
                    Conc=('aligned_particles_L',lambda x:float(np.quantile(x,mode))),
                    n_effects=('aligned_particles_L','size'),
                    n_studies=('DOI','nunique'),
                    min_effect=('aligned_particles_L','min'),
                    max_effect=('aligned_particles_L','max')
                )
                datasets[collapse]=by_species
                by_species.to_csv(args.output/f'ssdtools_input_{prefix}_{collapse}.csv',index=False,encoding='utf-8-sig')
            for collapse,df in datasets.items():
                if len(df)<5:
                    warnings.warn(f'{environment}/{mechanism}/{collapse} has only {len(df)} species')
                    continue
                hc=bootstrap_hc(df['Conc'].to_numpy(),nboot=args.nboot)
                # TR1 = lower CI of TR2, NOT a separately fitted fourth SSD.
                if collapse=='q1_all':
                    rows=[('TR1','HC5 95% CI lower',hc['0.05']['lcl'],hc['0.05']),
                          ('TR2','HC5 point',hc['0.05']['est'],hc['0.05'])]
                else:
                    rows=[('TR3','HC5 point',hc['0.05']['est'],hc['0.05']),
                          ('TR4','HC10 point',hc['0.1']['est'],hc['0.1'])]
                for tr,method,val,ci in rows:
                    results.append({'Environment':environment,'Mechanism':mechanism,
                     'TR':tr,'Definition':method,'SSD_Collapse':collapse,'N_species':len(df),
                     'N_input_effects':int(df['n_effects'].sum()),
                     'particles_L':val,'mg_L_proxy':val*dfa['ambient_mean_particle_mass_mg_proxy'].iloc[0],
                     'hc_ci_lcl':ci['lcl'],'hc_ci_ucl':ci['ucl'],
                     '95CI':f"{ci['lcl']:.6g} - {ci['ucl']:.6g}",
                     'classification':'EXPLORATORY_ONLY_NOT_PUBLISHED_REPRODUCTION'})
            plot_ssd(datasets, f'{environment} / {mechanism} / {args.size_min:g}-{args.size_max:g} µm',args.output/f'ssd_{prefix}.png')
    out=pd.DataFrame(results)
    out.to_csv(args.output/'TR1_TR4_exploratory.csv',index=False,encoding='utf-8-sig')
    print('==== AUDIT ====',json.dumps(audit,ensure_ascii=False,indent=2))
    if len(out):print('\n==== TR1-4 (EXPLORATORY PROXY, NOT PUBLISHED VALUES) ====\n',out[['Environment','Mechanism','TR','N_species','particles_L','95CI']].to_string(index=False))
    print('\nOutput directory:',args.output.resolve())


if __name__=='__main__':main()
