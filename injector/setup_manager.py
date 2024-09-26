import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
import astropy.units as u

from .io_tools import FilterbankReader, print_exe
from .binary_model import PulsarBinaryModel
from .signal_generator import PulsarSignal
from .observatory import Observation



class SetupManager:
    def __init__(self, pulsar_data_path, filterbank_path, ephem_path, output_path):
        self.fb = self.get_filterbank(filterbank_path) 
        self.pulsars = self.get_pulsars(pulsar_data_path)
        self.ephem = self.get_ephem(ephem_path)
        self.output_path = output_path

        self.pulsar_models = self.construct_models()
        self.parfile_paths = self.create_par_files()
        
    @staticmethod
    def get_filterbank(filterbank_path):
        try:
            fb = FilterbankReader(filterbank_path)
        except FileNotFoundError:
            sys.exit(f'Unable to open filterbank file: {filterbank_path}')
        else:
            return fb
        
    def construct_models(self):
        pulsar_models = []
        for pulsar_data in self.pulsars:
            obs = Observation(self.fb, self.ephem, pulsar_data, validate=True)
            binary = PulsarBinaryModel(pulsar_data, validate=True)
            pulsar_models.append(PulsarSignal(obs, binary, pulsar_data, validate=True))

        return pulsar_models

    @staticmethod
    def get_pulsars(pulsar_data_path):
        try:
            with open(pulsar_data_path, 'r') as file:
                read_inject_file = json.load(file)
        except FileNotFoundError:
            sys.exit(f'Unable to find {pulsar_data_path}.')
        except json.JSONDecodeError:
            sys.exit(f'Unable to parse {pulsar_data_path} using JSON.')

        pulsar_list = read_inject_file.get('pulsars', None)
        help=r'See "example.inject" file in https://github.com/erc-compact/inception/tree/main/injector.'
        if pulsar_list:
            for i, pulsar in enumerate(pulsar_list):
                if not pulsar.get('ID', None):
                    sys.exit(f'No "ID" key word found in pulsar {i+1} in {pulsar_data_path}. {help}')
            return pulsar_list
        else:
            sys.exit(f'No "pulsars" key word found in {pulsar_data_path}. {help}')

    @staticmethod
    def get_ephem(ephem):
        if ephem != 'builtin':
            try:
                import jplephem
            except ImportError:
                print_exe('jplephem package not installed. Using built-in ephemeris... ')
                ephem = 'builtin'
        return ephem

    def create_par_files(self):
        parfile_paths = []
        for i in range(len(self.pulsars)):
            parfile_paths.append(self.parfile_creator(i))
        return parfile_paths

    @staticmethod
    def source2str(source):
        ra_hms = np.array([*source.ra.hms])
        dec_dms = np.array([*source.dec.dms])
        ra_sign = '-' if np.sign(ra_hms[0]) == -1 else ''
        dec_sign = '-' if np.sign(dec_dms[0]) == -1 else ''

        ra_str = '{}{:02.0f}:{:02.0f}:{:07.4f}'.format(ra_sign, *np.abs(ra_hms))
        dec_str = '{}{:02.0f}:{:02.0f}:{:07.4f}'.format(dec_sign, *np.abs(dec_dms))
        return ra_str, dec_str

    def parfile_creator(self, i):
        pulsar_model = self.pulsar_models[i]
        parfile_params = {'PSR': f'0000+{i+1:04}i'}

        parfile_params['RAJ'], parfile_params['DECJ'] = self.source2str(pulsar_model.obs.source)
        parfile_params['POSEPOCH'] = pulsar_model.posepoch

        parfile_params['DM'] = pulsar_model.obs.DM

        parfile_params['PEPOCH'] = pulsar_model.pepoch
        for i, freq_deriv in enumerate(pulsar_model.FX_list):
            if freq_deriv != 0:
                parfile_params[f'F{i}'] = str(freq_deriv).replace('e', 'D')
            
        ephem = Path(pulsar_model.obs.ephem).stem.upper()
        parfile_params['EPHEM'] = ephem if (ephem != 'BUILTIN') else 'DE440'

        parfile_params['TZRMJD'] = pulsar_model.obs.obs_start_bary
        parfile_params['TZRFRQ'] = 0

        parfile_params['CLK'] = 'TT(BIPM)'
        parfile_params['UNITS'] = 'TDB'
        parfile_params['TIMEEPH'] = 'FB90'
        parfile_params['T2CMETHOD'] = 'TEMPO'
        parfile_params['CORRECT_TROPOSPHERE'] = 'N'
        parfile_params['PLANET_SHAPIRO'] = 'N'
        parfile_params['DILATEFREQ'] = 'N'

        if pulsar_model.binary.period != 0:
            parfile_params['BINARY'] = 'BT'
            parfile_params['T0'] = pulsar_model.binary.T0
            parfile_params['A1'] = pulsar_model.binary.A1
            parfile_params['PB'] = pulsar_model.binary.period * u.s.to(u.day)
            parfile_params['ECC'] = pulsar_model.binary.e
            parfile_params['OM'] =  np.rad2deg(pulsar_model.binary.AoP)

        par_file = pd.Series(parfile_params)
        par_file_path = self.output_path+f'/{pulsar_model.ID}.par'
        par_file.to_csv(par_file_path, sep='\t', header=False)

        return par_file_path
