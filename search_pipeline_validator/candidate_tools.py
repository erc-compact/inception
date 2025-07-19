import json
import glob
import numpy as np
import pandas as pd
from pathlib import Path
import astropy.units as u
import astropy.constants as const

import xml.etree.ElementTree as ET

import pipeline_tools as inj_tools

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))




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


def xml2csv(xml_file, output, save=True):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    xml_dict = xml_to_dict(root)
    csv_cands = pd.DataFrame(xml_dict['peasoup_search']['candidates']['candidate'])
    csv_cands = csv_cands.astype(np.float64)[['period', 'dm', 'acc', 'snr']]
    if save:
        csv_cands.to_csv(output)
    return csv_cands


def par_cand2csv(injection_report, work_dir, output, match=[]):
    psr_candfiles = []
    for psr in injection_report['pulsars']:
        cand_file = glob.glob(f"{work_dir}/{psr['ID']}*.cands")[0]
        cand_df = pd.read_csv(cand_file, skiprows=11, engine='python', sep=r'\s+').iloc[0]
        fold_pars = [psr['ID'], *cand_df[['f0_new', 'dm_new', 'acc_new', 'S/N_new']].values]
        psr_candfiles.append(fold_pars)
    
    df_cands = pd.DataFrame(psr_candfiles, columns=['ID', 'f0', 'dm', 'acc', 'SNR'])
    if match:
        matched = []
        report, tol = match
        for i, row in df_cands.iterrows():
            ID, f0, dm, acc, SNR = row
            psr = report['pulsars'][i]

            p_cond = np.abs(psr['PX'][0]-1/f0)/psr['PX'][0] <= tol['p0']
            dm_cond = np.abs(psr['DM'] - dm) <= tol['dm']
            snr_cond= (SNR/psr['SNR'] <= tol['snr_upp']) and (SNR/psr['SNR'] >= tol['snr_low'])

            if p_cond and dm_cond and snr_cond:
                matched.append(True)
            else:
                matched.append(False)
        
        df_cands['match'] = matched

    df_cands.to_csv(output)


def fold_cand2csv(cand_file, output):
    candidates = []
    cand_df = pd.read_csv(cand_file, skiprows=11, engine='python', sep=r'\s+')
    for _, row in cand_df.iterrows():
        fold_pars = row[['#id', 'f0_new', 'dm_new', 'acc_new', 'S/N_new']].values
        candidates.append(fold_pars)
    
    df_cands = pd.DataFrame(candidates, columns=['ID', 'f0', 'dm', 'acc', 'SNR'])
    df_cands.to_csv(output)


def freq_correction(pulsar, pepoch_ref=0.5):
    obs = pulsar.obs
    topo_mjd_mid = obs.obs_start +  obs.obs_len*pepoch_ref * u.s.to(u.day)

    radial_velocity = obs.earth_radial_velocity(topo_mjd_mid)

    if pulsar.binary.period:
        bary_sec = obs.topo2bary(np.array([topo_mjd_mid]), mjd=False, interp=False) + pulsar.orbit_ref
        binary_radial_velocity = pulsar.binary.get_radial_velocity_coord(bary_sec)

        radial_velocity += binary_radial_velocity

    doppler_shift =  1 - radial_velocity/const.c.value
    return doppler_shift[0]


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


def correct_fftsize_offset(period, acc, fftsize, nsamples, dt):
    pdot = acc * period / const.c.value
    return period - pdot * (fftsize - nsamples) * dt / 2


def fit_orbit(pulsar, pepoch_ref=0.5, mode='accel', coord_frame=True):
    from scipy.optimize import curve_fit

    topo_sec = np.linspace(0, pulsar.obs.obs_len, 1000)
    topo_mjd = pulsar.obs.sec2mjd(topo_sec)

    if mode == 'accel':
        def fit_func(t, vel, acc):
            return vel + acc*t
        
    elif mode == 'jerk':
        def fit_func(t, vel, acc, jerk):
            return vel + acc*t + 0.5*jerk*t**2

    r_v = np.zeros_like(topo_sec)
    if pulsar.binary.period:
        bary_sec = pulsar.obs.topo2bary(topo_mjd, mjd=False, interp=False) + pulsar.orbit_ref
        r_v += pulsar.binary.get_radial_velocity_coord(bary_sec)
    if coord_frame:
        r_v += pulsar.obs.earth_radial_velocity(topo_mjd)

    fit_pars = curve_fit(fit_func, topo_sec - pulsar.obs.obs_len*pepoch_ref, r_v)

    vel_max_ind = topo_sec[np.argmax(r_v)]/pulsar.obs.obs_len
    vel_min_ind = topo_sec[np.argmin(r_v)]/pulsar.obs.obs_len
    return fit_pars, vel_max_ind, vel_min_ind


def create_cand_file_acc(cands, cand_file_path):

    with open(cand_file_path, 'w') as file:
        file.write("#id DM accel F0 F1 S/N\n")
        for i, cand in cands.iterrows():
            file.write(f"{i} {cand['dm']} {cand['acc']} {1/cand['period']} 0 {cand['snr']}\n")


def merge_cand_file_acc(cand_list, cand_file_path):

    cand_data = []
    for cand_file in cand_list:
        data = pd.read_csv(cand_file, delimiter=' ')
        for i, row in data.iterrows():
            cand_data.append(row.values[1:].tolist())

    with open(cand_file_path, 'w') as file:
        file.write("#id DM accel F0 F1 S/N\n")
        for i, cand in enumerate(cand_data):
            file.write(f"{i} {cand[0]} {cand[1]} {cand[2]} {cand[3]} {cand[4]}\n")


class CandMatcher:
    def __init__(self, injection_report, candidates, filterbank, fftsize, ephem, corr_period=False, override_length=0):
        from injector.io_tools import FilterbankReader
        from injector.setup_manager import SetupManager

        self.cands = pd.read_csv(candidates, index_col=0) if type(candidates) == str else candidates

        self.fb = FilterbankReader(filterbank, stats_samples=0)
        self.fftsize = fftsize

        self.setup = SetupManager(injection_report, filterbank, ephem, generate=False, override_length=override_length)

        if corr_period:
            self.correct_periods(fftsize, override_length)

    def correct_periods(self, fftsize, override_length):
        n_samples = override_length if override_length else self.fb.n_samples
        new_periods = correct_fftsize_offset(self.cands['period'], self.cands['acc'], fftsize, n_samples, self.fb.dt)

        self.cands = self.cands.rename(columns={'period': 'period_input'})
        self.cands['period'] = new_periods


    def match_candidates(self, pepoch_ref=0.5, snr_limit=3, max_harmonic=4):
        self.cands['cand_id'] = np.arange(len(self.cands))

        pulsar_cands = {}
        for pm in self.setup.pulsar_models:

            pulsar_acc_fit, vel_max_ind, vel_min_ind = fit_orbit(pm, pepoch_ref=pepoch_ref, mode='accel')

            doppler_shift_ref = freq_correction(pm, pepoch_ref=pepoch_ref)
            doppler_shift_min_vel= freq_correction(pm, pepoch_ref=vel_min_ind)
            doppler_shift_max_vel = freq_correction(pm, pepoch_ref=vel_max_ind)

            fft_bin = 1/(self.fftsize*pm.obs.dt)
            F0 = pm.FX_list[0]
            cands_F = (1/self.cands['period'])
            harmonic_div = np.round(((cands_F/F0)))    
            harmonic_div[harmonic_div > max_harmonic] = 1
            cands_F /= harmonic_div

            nbins_offset = (F0*doppler_shift_ref - cands_F) / fft_bin

            F_min, F_max = min(F0*doppler_shift_min_vel, F0*doppler_shift_max_vel), max(F0*doppler_shift_min_vel, F0*doppler_shift_max_vel)
            freq_cond = ((cands_F >= F_min-fft_bin) & (cands_F <= F_max+fft_bin)) 
            
            DM_limit = DM_curve(pm, snr_limit)
            dm_offset =  (pm.prop_effect.DM - self.cands['dm'])
            snr_offset = (pm.SNR - self.cands['snr'])
            dm_cond = np.abs(dm_offset) <= DM_limit

            
            accel_drift = (pulsar_acc_fit[0][1] - self.cands['acc']) * F0 * doppler_shift_ref / (const.c.value * fft_bin**2)


            candidates = self.cands[freq_cond & dm_cond]
            candidates['nbins_offset'] = nbins_offset[freq_cond & dm_cond]
            candidates['accel_bin_drift'] = accel_drift[freq_cond & dm_cond]
            candidates['dm_offset'] = dm_offset[freq_cond & dm_cond]
            candidates['snr_offset'] = snr_offset[freq_cond & dm_cond]
    
            candidates_sorted = candidates.sort_values(by='snr', key=abs, ascending=False)

            pulsar_cands[pm.ID] = candidates_sorted

        return pulsar_cands


    def generate_files(self, candidate_root, max_cand_per_inj=-1, pepoch_ref=0.5, snr_limit=3, max_harmonic=4, create_candfile=True):
        pulsar_cands = self.match_candidates(pepoch_ref=pepoch_ref, snr_limit=snr_limit, max_harmonic=max_harmonic)
        cands_data = []
        for pm in self.setup.pulsar_models:
            candidates = pulsar_cands[pm.ID].reset_index(drop=True)
            for i, row in candidates.iterrows():
                if (i < max_cand_per_inj) or (max_cand_per_inj == -1):
                    cands_data.append([pm.ID, *row.values])
                else:
                    break
        
        cands_df = pd.DataFrame(cands_data, columns=['inj_id', *pulsar_cands[next(iter(pulsar_cands))].columns])
        cands_df.to_csv(f'{candidate_root}.csv')
        if create_candfile:
            create_cand_file_acc(cands_df, f'{candidate_root}.candfile')



