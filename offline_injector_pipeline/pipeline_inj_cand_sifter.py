import os
import re
import json
import numpy as np
import pandas as pd
from pathlib import Path
import astropy.units as u
import astropy.constants as const
from scipy.optimize import fsolve
import xml.etree.ElementTree as ET

from inception.injector.setup_manager import SetupManager
from inception.injector.io_tools import print_exe


class CandFinder:
    def __init__(self, fb, injection_plan, injection_report):
        self.work_dir = os.getcwd()

        self.inj_report, self.beam_id = self.parse_report(injection_report)
        self.inj_plan = self.resolve_inj_plan(injection_plan)

        print_exe('Setting up pulsar models ...')
        self.setup = SetupManager(self.inj_plan, fb, generate=True)
        print_exe('Done. Sifting ...')

    def filter_df(self, df, snr_limit=5, pfact=1, adjust=0, period_key='period'):
        pulsar_cands = []
        for pm in self.setup.pulsar_models:
            period = pm.PX_list[0]*pfact
            dm = pm.prop_effect.DM
            vc = self.get_doppler_velocity(pm)
            d_dm = self.get_DM_smear(pm, snr_limit=snr_limit, pfact=pfact)[0]

            doppler_min, doppler_max = min(1-vc, 1+vc)*(1-adjust), max(1-vc, 1+vc)*(1+adjust)
            p_cond = (df[period_key] > period * doppler_min) & (df[period_key] < period * doppler_max)
            dm_cond = (df['dm'] > dm - d_dm) & (df['dm'] < dm + d_dm)

            psr_cands = df[p_cond & dm_cond]
            psr_cands['pulsar_id'] = pm.ID
            psr_cands['delta_dm'] = d_dm
            psr_cands['doppler_max'] = period * doppler_max
            psr_cands['doppler_min'] = period * doppler_min
            pulsar_cands.append(psr_cands)
        return pd.concat(pulsar_cands)
        
    @staticmethod
    def parse_report(json_file):
        with open(json_file, 'r') as file:
            pars = json.load(file)
        
        fb_name = pars['injection_report']['fb']
        match = re.search(r'cfbf\d{5}', fb_name)
        return pars, match.group()
    
    def resolve_inj_plan(self, json_file):
        with open(json_file, 'r') as file:
            pars = json.load(file)

        seeded_inject_file = f'{self.work_dir}/{Path(json_file).name}'
        pars['psr_global']['global_seed'] = self.inj_report['injection_report']['global_seed']
        with open(seeded_inject_file, 'w') as file:
            json.dump(pars, file, indent=4)

        return seeded_inject_file
    
    @staticmethod
    def get_doppler_velocity(pulsar_model):
        obs_arr = np.linspace(0, pulsar_model.obs.obs_len, 10**4)
        dt = obs_arr[1] - obs_arr[0]
        timeseries = obs_arr*u.s.to(u.day) + pulsar_model.obs.obs_start

        earth_delays = pulsar_model.obs.topo2bary_calc(timeseries, mjd=True)*const.c.value
        earth_vel = np.gradient(earth_delays, dt)

        binary_delays = pulsar_model.binary.orbital_delay(obs_arr)*const.c.value
        binary_vel = np.gradient(binary_delays, dt)

        vc = (np.max(np.abs(binary_vel))+np.max(np.abs(earth_vel)))/const.c.value
        return vc
    
    @staticmethod
    def get_DM_smear(pulsar_model, snr_limit, pfact=1):

        def relative_SNR(d_DM):
            duty_cycle = pulsar_model.pulsar_pars['duty_cycle']
            period = pulsar_model.PX_list[0] * pfact
            DM_const = pulsar_model.prop_effect.DM_const

            W_int = period * duty_cycle
            freq = pulsar_model.obs.freq_arr
            W_eff = np.sqrt(W_int**2 + (2*DM_const * d_DM * (1/freq[0]**2 - 1/freq[-1]**2))**2)
            return np.sqrt(np.clip((period-W_eff)/W_eff, 0, None))
        
        def get_DM_step(d_DM):
            return relative_SNR(d_DM)/relative_SNR(0) * pulsar_model.SNR - snr_limit
        
        d_DM = fsolve(get_DM_step, 1e-3)
        return np.abs(d_DM)

    def parse_csv_file(self, cand_csv):
        inj_rows = cand_csv[cand_csv['beam_id'] == self.beam_id]
        inj_rows.sort_values('snr')
        inj_rows['beam_index'] = np.arange(len(inj_rows))
        return inj_rows
    
    def parse_xml_file(self, xml_file):
        tree = ET.parse(xml_file)
        root = tree.getroot()
        xml_dict = self.xml_to_dict(root)
        xml_df = pd.DataFrame(xml_dict['peasoup_search']['candidates']['candidate'])
        return xml_df.astype(np.float64)
    
    def xml_to_dict(self, element):
        if not list(element):
            result = element.text.strip() if element.text else None
            if element.attrib:
                result = {'@' + k: v for k, v in element.attrib.items()}
                if element.text and element.text.strip():
                    result['#text'] = element.text.strip()
            return {element.tag: result}
        
        result = {}
        for child in element:
            child_dict = self.xml_to_dict(child)
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
    