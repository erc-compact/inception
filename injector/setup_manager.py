import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path 
import astropy.units as u

from .io_tools import FilterbankReader, print_exe, str2func
from .binary_model import BinaryModel
from .pulsar_model import PulsarModel
from .observation import Observation



class SetupManager:
    def __init__(self, pulsar_data_path, filterbank_path, ephem_path, output_path):
        self.fb = self.get_filterbank(filterbank_path) 
        self.pulsars = self.get_pulsars(pulsar_data_path)
        self.ephem = self.get_ephem(ephem_path)
        self.output_path = output_path

        self.pulsar_models = self.construct_models()
        self.parfile_paths = self.create_parfiles()
        self.mode_resolver()
        
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
            obs = Observation(self.fb, self.ephem, pulsar_data, generate=False)
            binary = BinaryModel(pulsar_data, generate=False)
            pulsar_models.append(PulsarModel(obs, binary, pulsar_data, generate=False))

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

    @staticmethod
    def source2str(source):
        ra_hms = np.array([*source.ra.hms])
        dec_dms = np.array([*source.dec.dms])
        ra_sign = '-' if np.sign(ra_hms[0]) == -1 else ''
        dec_sign = '-' if np.sign(dec_dms[0]) == -1 else ''

        ra_str = '{}{:02.0f}:{:02.0f}:{:07.4f}'.format(ra_sign, *np.abs(ra_hms))
        dec_str = '{}{:02.0f}:{:02.0f}:{:07.4f}'.format(dec_sign, *np.abs(dec_dms))
        return ra_str, dec_str

    def create_parfiles(self):
        parfile_paths = []
        for i in range(len(self.pulsars)):
            pulsar_model = self.pulsar_models[i]
            parfile_params = {'PSR': f'0000+{i+1:04}i'}

            parfile_params['RAJ'], parfile_params['DECJ'] = self.source2str(pulsar_model.obs.source)
            parfile_params['POSEPOCH'] = pulsar_model.posepoch

            parfile_params['DM'] = pulsar_model.prop_effect.DM

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

            parfile_paths.append(par_file_path)

        return parfile_paths
    
    def mode_resolver(self):
        for i, pulsar_pars in enumerate(self.pulsars):
            mode = pulsar_pars.get('mode', 'python')
            polycos_path = pulsar_pars.get('polycos', None)
            
            if mode not in ['pint', 'python']:
                sys.exit(f"Invalid mode for pulsar {pulsar_pars['ID']}. Must be either 'python' or 'pint'")
            if polycos_path:
                mode = 'pint'
            self.pulsars[i]['mode'] = mode

            if mode == 'pint':
                try:
                    import pint.logging as logging      # type: ignore
                    _ = logging.setup('ERROR')  
                    import pint.models as models        # type: ignore
                    from pint.polycos import Polycos    # type: ignore
                except ImportError:
                    sys.exit('pint-pulsar package not installed, cannot use polycos.')
                else:
                    
                    if not polycos_path:
                        polyco_path = self.polycos_creator(self.parfile_paths[i], pulsar_pars, self.pulsar_models[i].obs,  pint_func=[models, Polycos])
                        self.pulsars[i]['polycos'] = polyco_path

    def polycos_creator(self, par_file, pulsar_pars, obs, pint_func): 
        models, Polycos = pint_func
        timing_model = models.get_model(par_file, EPHEM=obs.ephem)

        t_mid = obs.obs_start + obs.obs_len/2 * u.s.to(u.day)
        polco_range = obs.obs_len/2 + 10*u.min.to(u.s)
        start, end = t_mid - polco_range*u.s.to(u.day), t_mid + polco_range*u.s.to(u.day)
        
        polycos_coeff = str2func(pulsar_pars.get('pint_N', 12), 'pint_N', self.ID, int)
        polycos_tspan = str2func(pulsar_pars.get('pint_T', 5), 'pint_T', self.ID, float)
        polycos_coeff = max(1, polycos_coeff)
        polycos_tspan = max(1, polycos_tspan) # minutes
        gen_poly = Polycos.generate_polycos(timing_model, start, end, obs.tempo_id, 
                                            polycos_tspan, polycos_coeff, 
                                            obs.f0, progress=False)
        
        polycos_path = Path(par_file).with_suffix('.polycos')
        gen_poly.write_polyco_file(polycos_path)
        return polycos_path