import re
import json
import glob
import numpy as np
import pandas as pd
from pathlib import Path
import astropy.units as u
from math import factorial
import astropy.constants as const

import xml.etree.ElementTree as ET


def xml_to_dict(element):
    if not list(element):
        result = element.text.strip() if element.text else None
        if element.attrib:
            result = {'@' + k: v for k, v in element.attrib.items()}
            if element.text and element.text.strip():
                result['#text'] = element.text.strip()
        return {element.tag: result}
    
    result = {}
    for child in element:
        child_dict = xml_to_dict(child)
        tag, value = list(child_dict.items())[0]

        if tag in result:
            if not isinstance(result[tag], list):
                result[tag] = [result[tag]]
            result[tag].append(value)
        else:
            result[tag] = value

    if element.attrib:
        result.update(('@' + k, v) for k, v in element.attrib.items())

    return {element.tag: result}


def xml2csv(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    xml_dict = xml_to_dict(root)
    peasoup = xml_dict.get('peasoup_search') or {}
    candidates_dict = peasoup.get('candidates') or {}
    candidates = candidates_dict.get('candidate')
    
    if candidates is None:
        csv_cands = pd.DataFrame(columns=['period', 'dm', 'acc', 'snr'])
    elif not isinstance(candidates, list):
        candidates = [candidates]
        csv_cands = pd.DataFrame(candidates)
        csv_cands = csv_cands.astype(np.float64)[['period', 'dm', 'acc', 'snr']]
    else:
        csv_cands = pd.DataFrame(candidates)
        csv_cands = csv_cands.astype(np.float64)[['period', 'dm', 'acc', 'snr']]

    pepoch = float(xml_dict['peasoup_search'].get('segment_parameters', {'segment_pepoch': None})['segment_pepoch'])
    fftsize = float(xml_dict['peasoup_search']['search_parameters']['size'])
    dt = float(xml_dict['peasoup_search']['header_parameters']['tsamp'])
    
    return csv_cands, pepoch, fftsize, dt


def add_PSR_rv_curve(pulsar, time, rv_seg, pepoch_ref):

    tref = pulsar.obs.obs_len * pepoch_ref
    dt = time - tref

    if pulsar.binary.period:
        rv_seg += pulsar.binary.get_radial_velocity_coord(time + pulsar.orbit_ref)
    else:
        v_deriv = pulsar.AX_list
        rv_seg = 0
        for n, deriv in enumerate(v_deriv):
            rv_seg += deriv * dt**(n + 1) / factorial(n + 1)

    return rv_seg

def get_freq_bounds(rv, pulsar, pepoch):
    F0 = pulsar.FX_list[0]
    obs = pulsar.obs

    if pulsar.pulsar_pars['frame'] == 'bary':
        F0 *= (1 - obs.earth_radial_velocity(pepoch)/const.c.value)

    rv_min = np.min(rv)
    rv_max = np.max(rv)

    F_max = F0 * (1 - rv_min / const.c.value)
    F_min = F0 * (1 - rv_max / const.c.value)

    return F_min, F_max

def DM_curve(pulsar, snr_limit):
    period = pulsar.PX_list[0]
    snr = pulsar.pulsar_pars['SNR']
    W_int = pulsar.pulsar_pars['duty_cycle'] * period

    snr_scale = snr/np.sqrt((period-W_int)/W_int)

    W_eff = period / (1 + (snr_limit/snr_scale)**2)
    freq = pulsar.obs.freq_arr
    freq_diff = np.abs(1/freq[0]**2 - 1/freq[-1]**2)

    DM_range = np.sqrt(W_eff**2 - W_int**2) / (pulsar.prop_effect.DM_const * freq_diff)
    return DM_range


def create_PULSARX_candfile(cands, candfile_path):

    with open(candfile_path, 'w') as file:
        file.write("#id DM accel F0 F1 S/N\n")
        for i, cand in cands.iterrows():
            file.write(f"{i} {cand['dm']} {cand['acc']} {1/cand['period']} 0 {cand['snr']}\n")


def pulsarx_par2csv(injection_report, results_dir):
    psr_candfiles = []
    for psr in injection_report:
        cand_file = glob.glob(f"{results_dir}/{psr['ID']}*.cands")[0]
        cand_df = pd.read_csv(cand_file, skiprows=11, engine='python', sep=r'\s+').iloc[0]
        fold_pars = [psr['ID'], *cand_df[['f0_new', 'f0_err', 'dm_new', 'dm_err', 'acc_new', 'acc_err', 'S/N_new', 'boxcar_width']].values]
        psr_candfiles.append(fold_pars)
    
    df_cands = pd.DataFrame(psr_candfiles, columns=['PSR_ID', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width'])
    return df_cands


def pulsarx_cand2csv(cand_file):
    candidates = []
    cand_df = pd.read_csv(cand_file, skiprows=11, engine='python', sep=r'\s+')
    for _, row in cand_df.iterrows():
        fold_pars = row[['#id', 'f0_new', 'f0_err', 'dm_new', 'dm_err', 'acc_new', 'acc_err', 'S/N_new', 'boxcar_width']].values
        candidates.append(fold_pars)
    
    df_cands = pd.DataFrame(candidates, columns=['cand_ID', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width'])
    return df_cands


def correct_fftsize_offset(period, acc, fftsize, nsamples, dt):
    pdot = acc * period / const.c.value
    return period - pdot * (fftsize - nsamples) * dt / 2


def presto_bestprof2csv(injection_report, results_dir):
    NUM = r'[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?'
    psr_folds = []
    for psr in injection_report:
        bestprof_file = glob.glob(f"{results_dir}/{psr['ID']}*.bestprof")[0]
        with open(bestprof_file) as file:
            text = file.read()

        p_ms, p_err_ms = re.search(rf'P_topo \(ms\)\s*=\s*({NUM})\s*\+/-\s*({NUM})', text).groups()
        pdot, pdot_err = re.search(rf"P'_topo \(s/s\)\s*=\s*({NUM})\s*\+/-\s*({NUM})", text).groups()
        dm = re.search(rf'Best DM\s*=\s*({NUM})', text).group(1)
        sigma = re.search(rf'\(~({NUM})\s*sigma\)', text)

        P0 = float(p_ms) / 1000
        P0_err = float(p_err_ms) / 1000

        F0 = 1 / P0
        F0_err = P0_err / P0**2
        acc = const.c.value * float(pdot) / P0
        acc_err = const.c.value * float(pdot_err) / P0
        snr = float(sigma.group(1)) if sigma else 0.0

        psr_folds.append([psr['ID'], F0, F0_err, float(dm), 0.0, acc, acc_err, snr, 0.0])

    df_folds = pd.DataFrame(psr_folds, columns=['PSR_ID', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width'])
    return df_folds


def dspsr_best2csv(injection_report, results_dir):
    psr_folds = []
    for psr in injection_report:
        best_file = glob.glob(f"{results_dir}/{psr['ID']}*_dspsr.best")[0]
        with open(best_file) as file:
            rows = [line.split() for line in file if not line.startswith('#')]

        # rows: 0=BC_prd, 1=TC_prd, 2=garbage (pdmp.C writes a stale value here), 3=DM_val, 4=BC_freq, 5=width S/N
        DM, _, DM_err = (float(v) for v in rows[3])
        F0, F0_err = (float(v) for v in rows[4])
        width, snr = (float(v) for v in rows[5])

        psr_folds.append([psr['ID'], F0, F0_err, DM, DM_err, 0.0, 0.0, snr, width])

    df_folds = pd.DataFrame(psr_folds, columns=['PSR_ID', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width'])
    return df_folds
