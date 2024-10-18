import re
import sys
import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path 
import astropy.units as u
from sympy import lambdify, symbols, Function

from .io_tools import FilterbankReader, print_exe
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
        if len(ID_list) != len(set(ID_list)):
            sys.exit('Pulsar IDs must be unique.')

        pulsar_clean_list = []
        for i, pulsar in enumerate(pulsar_list):
            parser = PulsarParamParser(pulsar, global_list)
            pulsar_clean_list.append(parser.psr_pars)

        pulsar_clean_list = self.double_pulsar(pulsar_clean_list)
        
        return pulsar_clean_list

    def double_pulsar(self, pulsar_list, ID_list):
        for i in range(len((pulsar_list))):
            double_psr = pulsar_list[i]['double_pulsar']
            if double_psr:
                if double_psr in ID_list:
                    companion_psr = pulsar_list[ID_list.index(double_psr)]
                    # check if double pulsar keyword exists in companion_psr 
                    # pulsar_list[i]['key'] = 'value'
                    
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
    


class PulsarParamParser:
    def __init__(self, pulsar_pars, global_pars):
        self.get_argument_parser()
        self.parse_inputs(pulsar_pars, global_pars)
    
    def dict_to_args(self, params):
        args_list = []
        for key, value in params.items():
            args_list.append(f'--{key}')
            args_list.append(str(value))
        return args_list
    
    def str2func(value, par, id, func):
        try:
            converted_value = func(value)
        except ValueError:
            sys.exit(f'Error: Invalid {par} for pulsar {id}')
        return converted_value
    
    def parse_inputs(self, pulsar_pars, global_pars):
        psr_pars = global_pars.copy()
        psr_pars.update(pulsar_pars.copy())

        FX_list, PX_list = self.get_spin_params(psr_pars)

        psr_args_list = self.dict_to_args(psr_pars)
        psr_args, _ = self.parser.parse_known_args(psr_args_list)

        clean_psr_pars = vars(psr_args)
        clean_psr_pars['FX'] = FX_list
        clean_psr_pars['PX'] = PX_list
        
        self.psr_pars = clean_psr_pars
    
    def get_argument_parser(self):
        class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
            def add_usage(self, usage, actions, groups, prefix=None): pass

        parser = argparse.ArgumentParser(description='Pulsar parameters', formatter_class=CustomFormatter)
        parser.add_argument('--ID', metavar='(str)', required=True, type=float, help='Identifier for injected pulsar')
        
        parser.add_argument('--RAJ', metavar='(hh:mm:ss)', required=False, type=str, help='Right Ascension (J2000) (default: beam centre)')
        parser.add_argument('--DECJ', metavar='(dd:mm:ss)', required=False, type=str, help='Declination (J2000) (default: beam centre)')
        parser.add_argument('--separation', metavar='(arcmin)', required=False, default=0, type=float, help='Pulsar offset from beam centre')
        parser.add_argument('--position_angle', metavar='(deg)', required=False, default=0, type=float, help='Position angle of pulsar')
        parser.add_argument('--beam_fwhm', metavar='(arcmin)', required=False, default=0, type=float, help='FWHM of injected beam, required for (separation, position_angle)')

        parser.add_argument('--PEPOCH', metavar='(MJD)', required=False, type=float, help='Reference epoch for pulsar spin values (default: obseravtion start, barycentre)')
        parser.add_argument('--FX', metavar='(Hz)', required=False, type=float, help='Xth frequency derivative of pulsar spin')
        parser.add_argument('--PX', metavar='(sec)', required=False, type=float, help='Xth period derivative of pulsar spin')
        parser.add_argument('--phase_offset', metavar='(phase)', required=False, default=0, type=float, help='Phase offset from PEPOCH')
        parser.add_argument('--DM', metavar='(pc/cm^3)', required=True, type=float, help='Dispersion measure')
        parser.add_argument('--SNR', required=True, type=float, help='Injected signal-to-noise')
        parser.add_argument('--PSD', metavar='(file)', required=False, type=str, help='NumPy .npy file containing a pulsar power spectrum (1D)')
        parser.add_argument('--spectral_index', required=False, default=0, type=float, help='Spectral index of pulsar')
        parser.add_argument('--duty_cycle', metavar='(phase)', required=False, default=0.1, type=float, help='Duty cycle of default gaussian pulse profile')
        parser.add_argument('--profile', metavar='(file)', required=False, type=str, default='default', help='NumPy .npy or EPN .txt file containing a custom pulsar pulse profile (1D or 2D)')
        parser.add_argument('--micro_structure', metavar='(microsec)', required=False, default=0, type=float, help='Mean timescale of pulse microstructure')
        parser.add_argument('--scattering_time', metavar='(millisec)', required=False, default=np.inf, type=float, help='Scattering timescale due to ISM')
        parser.add_argument('--scattering_index', required=False, default=-4, type=float, help='Scattering index to describe frequency evolution')
        parser.add_argument('--DM_smear', metavar='(0 or 1)', required=False, default=0, type=int, help='Smear the pulse profile due to intra-channel DM smearing')

        parser.add_argument('--binary_period', metavar='(hour)', required=False, type=float, help='Period of binary oribit')
        parser.add_argument('--T0', metavar='(MJD)', required=False,  type=float, help='Reference epoch of pulsar periapsis (default: obseravtion start, barycentre)')
        parser.add_argument('--A1', metavar='(light-sec)', required=False, type=float, help='Projected semi-major orbital axis')
        parser.add_argument('--M1', metavar='(M_sun)', required=False,  type=float, help='Mass of pulsar (default: 1.4)')
        parser.add_argument('--M2', metavar='(M_sun)', required=False, type=float, help='Companion mass')
        parser.add_argument('--inc', metavar='(deg)', required=False, type=float, help='Orbital inclination (default: 90)')
        parser.add_argument('--ecc', required=False, default=0, type=float, help='Orbital eccentricity')
        parser.add_argument('--AoP', metavar='(deg)', required=False,  default=0, type=float, help='Argument of periapsis')
        parser.add_argument('--LoAN', metavar='(deg)', required=False,  default=0, type=float, help='Longitude of the ascending node')
        parser.add_argument('--double_pulsar', metavar='(ID)', required=False, type=str, help='Make pulsar the companion of pulsar {double_pulsar:ID}')

        parser.add_argument('--mode', metavar='(str)', required=False, default='python', type=str, help="Inject using analytical 'python' code or polycos from 'pint'")
        parser.add_argument('--pint_N', required=False, default=12, type=int, help='Number of coefficients per timestep for polycos generation')
        parser.add_argument('--pint_T', metavar='(min)', required=False, default=5, type=float, help='Timestep for polycos generation')
        parser.add_argument('--polycos', metavar='(file)', required=False, type=str, help='Polycos file defining Pulsar phase. If supplied, Pint will use this file instead of making one.')

        self.parser = parser
        
    @staticmethod
    def spin_pars_converter(spin_values, spin_types):
        FX_values = []
        PX_values = []
        t = symbols('t')

        def converter(spin_type):
            SX = Function(f'{spin_type}0')(t)
            SX_inv = 1/SX

            S_diff_funcs = [SX.diff(t, n) for n in range(len(spin_types))]
            S_diff_sym = symbols([f'S{x}' for x in range(len(spin_types))])
            S_diff_dict = dict(zip(S_diff_funcs, S_diff_sym))

            def deriv_calc(x, values):
                converter_func = lambdify([*S_diff_sym[:x+1]], SX_inv.diff(t, x).subs(S_diff_dict))
                return converter_func(*values)
            
            return deriv_calc
        
        P_converter = converter('P')
        F_converter = converter('F')

        for x in range(len(spin_types)):
            if spin_types[x][0] == 'F':
                FX_values.append(spin_values[x])
                PX_values.append(P_converter(x, FX_values))
            elif spin_types[x][0] == 'P':
                PX_values.append(spin_values[x])
                FX_values.append(F_converter(x, PX_values))

        return FX_values, PX_values
    
    def get_spin_params(self, pulsar_pars):
        ID = pulsar_pars['ID']
        p0_find = pulsar_pars.get('P0', 0) 
        f0_find = pulsar_pars.get('F0', 0)
        if (not p0_find) and (not f0_find):
            sys.exit(f'No P0 or F0 found for pulsar {ID}. Pulsar must be spinning!')

        pattern = re.compile(r'^[PF]\d+$')
        matching_keys = [key for key in pulsar_pars.keys() if pattern.match(key)]
        matching_keys.sort(key=lambda x: int(x[1:]))

        spin_types = []
        for i in range(int(matching_keys[-1][1:])+1):
            if (f'P{i}' in matching_keys) and (f'F{i}' not in matching_keys):
                spin_types.append(f'P{i}')
            else:
                spin_types.append(f'F{i}')
        spin_values = [self.str2func(pulsar_pars.get(key, 0), key, ID, float) for key in spin_types]

        return self.spin_pars_converter(spin_values, spin_types)