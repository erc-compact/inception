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

from pipeline_tools import glob_psr


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
        csv_cands = pd.DataFrame(columns=['period', 'dm', 'acc', 'snr', '@id'])
    elif not isinstance(candidates, list):
        candidates = [candidates]
        csv_cands = pd.DataFrame(candidates)
        csv_cands = csv_cands.astype(np.float64)[['period', 'dm', 'acc', 'snr', '@id']]
    else:
        csv_cands = pd.DataFrame(candidates)
        csv_cands = csv_cands.astype(np.float64)[['period', 'dm', 'acc', 'snr', '@id']]
    csv_cands = csv_cands.rename(columns={'@id': 'xml_id'}).astype({'xml_id': int})

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
        for n, deriv in enumerate(v_deriv):
            rv_seg += deriv * dt**(n + 1) / factorial(n + 1)

    return rv_seg

def get_freq_bounds(rv, pulsar):
    F0 = pulsar.FX_list[0]

    rv_min = np.min(rv)
    rv_max = np.max(rv)

    F_max = F0 * (1 - rv_min / const.c.value)
    F_min = F0 * (1 - rv_max / const.c.value)

    return F_min, F_max

def observed_channel_profiles(pulsar_model):
    pm = pulsar_model
    observed_profile = pm.get_observed_profile()

    phase = np.linspace(0, 1, pm.emission.profile_length, endpoint=False)

    return np.array([observed_profile(phase, chan) for chan in range(pm.obs.n_chan)])


def dm_response_curve(pulsar_model, dm_offsets, profiles=None):
    pm = pulsar_model
    obs = pm.obs

    if profiles is None:
        profiles = observed_channel_profiles(pm)

    nbins = profiles.shape[1]
    phase = np.linspace(0, 1, nbins, endpoint=False)

    inv_f2 = 1.0 / obs.freq_arr**2
    inv_f2 = inv_f2 - inv_f2[0]

    response = np.empty(len(dm_offsets))
    for i, ddm in enumerate(np.asarray(dm_offsets)):
        turns = pm.prop_effect.DM_const * ddm * inv_f2 / pm.period

        stacked = np.zeros(nbins)
        for chan in range(obs.n_chan):
            stacked += np.interp((phase - turns[chan]) % 1.0, phase,
                                 profiles[chan], period=1.0)

        response[i] = np.sqrt(np.sum((stacked - stacked.mean())**2))

    return response / response.max()


def _dm_crossing(offsets, response, level, i_above, i_below):
    frac = (response[i_above] - level) / (response[i_above] - response[i_below])
    return offsets[i_above] + frac * (offsets[i_below] - offsets[i_above])


def dm_match_bounds(pulsar_model, level=0.3):
    pm = pulsar_model
    obs = pm.obs

    band = abs(1.0 / obs.low_f**2 - 1.0 / obs.high_f**2)
    span = pm.period / (pm.prop_effect.DM_const * band)

    profiles = observed_channel_profiles(pm)

    for _ in range(10):
        offsets = np.linspace(-span, span, 201)
        response = dm_response_curve(pm, offsets, profiles)

        peak = int(np.argmax(response))
        below = response < level
        left = np.flatnonzero(below[:peak])
        right = np.flatnonzero(below[peak:])

        if left.size and right.size:
            low = _dm_crossing(offsets, response, level, left[-1] + 1, left[-1])
            high = _dm_crossing(offsets, response, level, peak + right[0] - 1, peak + right[0])
            return -low, high

        span *= 3

    return span, span


def create_PULSARX_candfile(cands, candfile_path):
    with open(candfile_path, 'w') as file:
        file.write("#id DM accel F0 F1 S/N\n")
        for i, cand in cands.iterrows():
            file.write(f"{i} {cand['dm']} {cand['acc']} {1/cand['period']} 0 {cand['snr']}\n")


def pulsarx_par2csv(injection_report, results_dir):
    psr_candfiles = []
    for psr in injection_report:
        cand_file = glob_psr(results_dir, psr['ID'], '.cands')
        if not cand_file:
            continue
        cand_df = pd.read_csv(cand_file[0], skiprows=11, engine='python', sep=r'\s+').iloc[0]
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


def get_freq_deriv(r, z, w, time_length):
    z0 = z - w / 2
    r0 = r - z0 / 2 - w / 6

    f0 = r0 / time_length
    f1 = z0 / time_length ** 2
    f2 = w / time_length ** 3

    return f0, f1, f2


def presto_sift2csv(sift_csv_path):
    df = pd.read_csv(sift_csv_path)
    if 'w' not in df.columns:
        df['w'] = 0.0

    if 'T' not in df.columns:
        df['T'] = df['r'] * df['P(ms)'] / 1000.0
    df['F0'], df['F1'], df['F2'] = get_freq_deriv(df['r'], df['z'], df['w'], df['T'])

    df['period'] = 1.0 / df['F0']

    df['F_match'] = df['r'] / df['T']

    parsed = df['file'].str.extract(r'_SEG_(?P<seg_i>\d+)_(?P<seg_n>\d+)_DS(?P<downsample>\d+)_DM[\d.]+')
    df['seg_i'] = parsed['seg_i'].astype(int)
    df['seg_n'] = parsed['seg_n'].astype(int)
    df['downsample'] = parsed['downsample'].astype(int)
    df['segment'] = df['seg_i'].astype(str) + '_' + df['seg_n'].astype(str)

    df = df.rename(columns={'DM': 'dm', 'SNR': 'snr'})
    return df


def _parse_bestprof(bestprof_file):
    NUM = r'[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?'
    with open(bestprof_file) as file:
        text = file.read()

    for frame in ('topo', 'bary'):
        period = re.search(rf'P_{frame} \(ms\)\s*=\s*({NUM})\s*\+/-\s*({NUM})', text)
        period_dot = re.search(rf"P'_{frame} \(s/s\)\s*=\s*({NUM})\s*\+/-\s*({NUM})", text)
        if period and period_dot:
            break
    else:
        raise ValueError(f'No topocentric or barycentric period in {bestprof_file}.')

    p_ms, p_err_ms = period.groups()
    pdot, pdot_err = period_dot.groups()
    dm = re.search(rf'Best DM\s*=\s*({NUM})', text).group(1)
    sigma = re.search(rf'\(~({NUM})\s*sigma\)', text)

    P0 = float(p_ms) / 1000
    P0_err = float(p_err_ms) / 1000

    F0 = 1 / P0
    F0_err = P0_err / P0**2
    acc = const.c.value * float(pdot) / P0
    acc_err = const.c.value * float(pdot_err) / P0
    snr = float(sigma.group(1)) if sigma else 0.0

    return {'F0': F0, 'F0_err': F0_err, 'DM': float(dm), 'DM_err': 0.0,
            'acc': acc, 'acc_err': acc_err, 'SNR': snr, 'width': 0.0}


def presto_bestprof2csv(injection_report, results_dir):
    psr_folds = []
    for psr in injection_report:
        bestprof_file = glob_psr(results_dir, psr['ID'], '.bestprof')
        if not bestprof_file:
            continue
        p = _parse_bestprof(bestprof_file[0])
        psr_folds.append([psr['ID'], p['F0'], p['F0_err'], p['DM'], p['DM_err'], p['acc'], p['acc_err'], p['SNR'], p['width']])

    df_folds = pd.DataFrame(psr_folds, columns=['PSR_ID', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width'])
    return df_folds


def presto_cand_bestprof2csv(psr_ids, results_dir):
    psr_folds = []
    for psr_id in psr_ids:
        bestprof_files = glob.glob(f"{results_dir}/{glob.escape(str(psr_id))}_CAND*.bestprof")
        if not bestprof_files:
            continue

        parsed = [(_parse_bestprof(f), f) for f in bestprof_files]
        p, best_file = max(parsed, key=lambda x: abs(x[0]['SNR']))
        psr_folds.append([psr_id, p['F0'], p['F0_err'], p['DM'], p['DM_err'], p['acc'], p['acc_err'], p['SNR'], p['width'],
                          Path(best_file).stem])

    df_folds = pd.DataFrame(psr_folds, columns=['PSR_ID', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width', 'fold'])
    return df_folds


def dspsr_best2csv(injection_report, results_dir):
    psr_folds = []
    for psr in injection_report:
        best_file = glob_psr(results_dir, psr['ID'], '_dspsr.best')
        if not best_file:
            continue
        best_file = best_file[0]
        with open(best_file) as file:
            rows = [line.split() for line in file if not line.startswith('#')]

        DM, _, DM_err = (float(v) for v in rows[3])
        F0, F0_err = (float(v) for v in rows[4])
        width, snr = (float(v) for v in rows[5])

        psr_folds.append([psr['ID'], F0, F0_err, DM, DM_err, 0.0, 0.0, snr, width])

    df_folds = pd.DataFrame(psr_folds, columns=['PSR_ID', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width'])
    return df_folds
