import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path 
import astropy.units as u

from .io_tools import FilterbankReader, print_exe
from .pulsar_par_parser import PulsarParParser
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

    def get_pulsars(self, pulsar_data_path):
        try:
            with open(pulsar_data_path, 'r') as file:
                read_inject_file = json.load(file)
        except FileNotFoundError:
            sys.exit(f'Unable to find {pulsar_data_path}.')
        except json.JSONDecodeError:
            sys.exit(f'Unable to parse {pulsar_data_path} using JSON.')

        pulsar_list = read_inject_file.get('pulsars', None)
        global_list = read_inject_file.get('psr_global', None)

        pulsar_list, ID_list = self.resolve_ID(pulsar_list, pulsar_data_path)

        pulsar_clean_list = []
        for pulsar in pulsar_list:
            parser = PulsarParParser(pulsar, global_list)
            pulsar_clean_list.append(parser.psr_pars)

        pulsar_clean_list = self.double_pulsar(pulsar_clean_list, ID_list)
        
        return pulsar_clean_list
    
    def resolve_ID(self, pulsar_list, pulsar_data_path):
        ID_list = []
        help="See 'example.inject' file in https://github.com/erc-compact/inception/tree/main/injector."
        if pulsar_list:
            for i, pulsar in enumerate(pulsar_list):
                if not pulsar.get('ID', None):
                    sys.exit(f'No "ID" key word found in pulsar {i+1} in {pulsar_data_path}. {help}')
                else:
                    ID_list.append(pulsar['ID'])
        else:
            sys.exit(f'No "pulsars" key word found in {pulsar_data_path}. {help}')


        non_dict_elements = [item for item in ID_list if not isinstance(item, dict)]
        if len(non_dict_elements) != len(set(non_dict_elements)):
            sys.exit('Pulsar IDs must be unique.')

        for i, ID in enumerate(ID_list):
            if type(ID) == dict:
                rng = np.random.default_rng(ID['seed'])
                seeds = rng.integers(0, 1e12, size=int(ID['replicate']))
                rng_pars = pulsar_list.pop(i)
                for j, seed in enumerate(seeds):
                    psr_pars = rng_pars.copy()
                    psr_pars['ID'] = f'pulsar_{j}'
                    psr_pars['seed'] = seed
                    psr_pars = self.resolve_random(psr_pars)
                    pulsar_list.append(psr_pars)

        return pulsar_list, ID_list

    def resolve_random(self, pulsar_pars):
        seed = pulsar_pars['seed']
        for key, value in pulsar_pars.items():
            rng = np.random.default_rng(seed)
            if type(value) == dict:
                units = value.get('units', 1)
                if units == 'T_obs':
                    units = self.fb.header['tsamp'] * self.fb.n_samples
                elif units == 'dt':
                    units = self.fb.header['tsamp']

                if value['rng'] == 'choice':
                    p = value.get('weights', np.ones_like(value['samples']))
                    pulsar_pars[key] = rng.choice(a=value['samples'], p=p/np.sum(p))
                elif value['rng'] == 'uniform':
                    pulsar_pars[key] = rng.uniform(low=value['low']*units, high=value['high']*units)
                elif value['rng'] == 'loguniform':
                    pulsar_pars[key] = np.exp(rng.uniform(low=np.log(value['low']*units), high=np.log(value['high']*units)))
                elif value['rng'] == 'normal':
                    pulsar_pars[key] = rng.normal(loc=value['mean']*units, scale=value['sigma']*units)

                elif value['rng'] == 'split_uniform':
                    p = value.get('weights', [1, 1])
                    which_range = rng.choice(a=[0, 1], p=p/np.sum(p))
                    val_range = value['lower_range'] if which_range == 0 else value['upper_range']
                    pulsar_pars[key] = rng.uniform(low=val_range[0]*units, high=val_range[1]*units)

                    binary = value.get('binary', [1, 1])
                    if not binary[which_range]:
                        pulsar_pars['binary_period'] = 0

        return pulsar_pars

    def double_pulsar(self, pulsar_list, ID_list):
        dble_psr_pars = ['RAJ', 'DECJ', 'separation', 'position_angle', 'beam_fwhm', 'DM', 'scattering_time', 'scattering_index', 'DM_smear',
                      'binary_period', 'T0', 'inc', 'ecc', 'LoAN']
        
        for i in range(len((pulsar_list))):
            double_psr = pulsar_list[i]['double_pulsar']
            if double_psr:
                if double_psr in ID_list:
                    c_psr = pulsar_list[ID_list.index(double_psr)]
                    if (not c_psr['double_pulsar']):
                        for par in dble_psr_pars:
                            pulsar_list[i][par] = c_psr[par]
                        pulsar_list[i]['M1'] = c_psr['M2']
                        pulsar_list[i]['M2'] = c_psr['M1']
                        pulsar_list[i]['AoP'] = c_psr['AoP'] + 180
                        pulsar_list[i]['A1'] = PulsarParParser.orbit_par_converter(c_psr['binary_period'], find='A1', 
                                                                                     M1=c_psr['M2'], M2=c_psr['M1'], inc=c_psr['inc']) 
                    else:
                         sys.exit(f"Only one 'double_pulsar' parameter allowed per binary pulsar pair.")                    
                else:
                    sys.exit(f"Invalid double_pulsar ID parameter for pulsar {pulsar_list[i]['ID']}.")
        return pulsar_list

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
            if self.pulsars[i]['create_parfile']:
                parfile_params = {'PSR': f'0000+{i+1:04}i'}

                parfile_params['RAJ'], parfile_params['DECJ'] = self.source2str(pulsar_model.obs.source)
                # parfile_params['POSEPOCH'] = pulsar_model.posepoch

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
                    parfile_params['A1'] = pulsar_model.binary.a1_sini_c
                    parfile_params['PB'] = pulsar_model.binary.period * u.s.to(u.day)
                    parfile_params['ECC'] = pulsar_model.binary.e
                    parfile_params['OM'] =  np.rad2deg(pulsar_model.binary.AoP+pulsar_model.binary.LoAN)

                par_file = pd.Series(parfile_params)
                par_file_path = self.output_path+f'/{pulsar_model.ID}.par'
                par_file.to_csv(par_file_path, sep='\t', header=False)

                parfile_paths.append(par_file_path)

        return parfile_paths
    
    def mode_resolver(self):
        for i, pulsar_pars in enumerate(self.pulsars):
            mode = pulsar_pars.get('mode', 'python')
            polycos_path = pulsar_pars['polycos']
            
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

        polycos_coeff = max(1, pulsar_pars['pint_N'])
        polycos_tspan = max(1, pulsar_pars['pint_T']) # minutes
        gen_poly = Polycos.generate_polycos(timing_model, start, end, obs.tempo_id, 
                                            polycos_tspan, polycos_coeff, 
                                            obs.f0, progress=False)
        
        polycos_path = Path(par_file).with_suffix('.polycos')
        gen_poly.write_polyco_file(polycos_path)
        return polycos_path
    