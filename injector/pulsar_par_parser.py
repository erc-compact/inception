import re
import sys
import argparse
import numpy as np
import astropy.constants as const
from scipy.optimize import fsolve
from sympy import lambdify, symbols, Function

from .binary_model import BinaryModel


class PulsarParParser:
    def __init__(self, pulsar_pars, global_pars):
        self.get_argument_parser()
        self.parse_inputs(pulsar_pars, global_pars)
    
    def get_argument_parser(self):
        class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
            def add_usage(self, usage, actions, groups, prefix=None): pass

        parser = argparse.ArgumentParser(description='Pulsar parameters in --signal', formatter_class=CustomFormatter)
        parser.add_argument('--ID', metavar='(str)', required=True, type=str, help='Identifier for injected pulsar')
        parser.add_argument('--seed', metavar='(positive int)', required=False, default=0, type=int, help='Random number generator seed for pulsar')
        parser.add_argument('--create_parfile', metavar='(0 or 1)', required=False, default='1', type=str, help="If '1', then create a TEMPO style parfile for injected pulsar")
        parser.add_argument('--frame', metavar='(topo or bary)', required=False, default='bary', type=str, help="barycentric/topocentric injection")

        parser.add_argument('--RAJ', metavar='(hh:mm:ss)', required=False, type=str, help='Right Ascension (J2000) (default: beam centre)')
        parser.add_argument('--DECJ', metavar='(dd:mm:ss)', required=False, type=str, help='Declination (J2000) (default: beam centre)')
        parser.add_argument('--separation', metavar='(arcmin)', required=False, default=0, type=float, help='Pulsar offset from beam centre')
        parser.add_argument('--position_angle', metavar='(deg)', required=False, default=0, type=float, help='Position angle of pulsar')
        parser.add_argument('--beam_fwhm', metavar='(arcmin)', required=False, default=0, type=float, help='FWHM of injected beam, required for (separation, position_angle)')

        parser.add_argument('--PEPOCH', metavar='(MJD)', required=False, type=float, help='Reference epoch for pulsar spin values (default: obseravtion start, barycentre)')
        parser.add_argument('--FX', metavar='(Hz)', required=False, type=list, help='Xth frequency derivative of pulsar spin')
        parser.add_argument('--PX', metavar='(sec)', required=False, type=list, help='Xth period derivative of pulsar spin')
        parser.add_argument('--AX', metavar='(m/s^(2+X))', required=False, type=list, help='Xth acceleration derivative')
        parser.add_argument('--ACCEPOCH', metavar='(obs phase)', required=False,  default=0.5, type=float, help='Reference epoch for AX, (0 - 1) T_obs')
        parser.add_argument('--presto_z', metavar='(z)', required=False, type=float, help='presto acceleration')
        parser.add_argument('--presto_w', metavar='(w)', required=False, type=float, help='presto jerk')
        parser.add_argument('--phase_offset', metavar='(phase)', required=False, default=0, type=float, help='Phase offset from PEPOCH')

        parser.add_argument('--DM', metavar='(pc/cm^3)', required=False, default=0, type=float, help='Dispersion measure')
        parser.add_argument('--SNR', metavar='(-)', required=True, type=float, help='Injected signal-to-noise')
        parser.add_argument('--PSD', metavar='(file)', required=False, type=str, help='NumPy .npy file containing a pulsar power spectrum (1D)')
        parser.add_argument('--spectral_index', metavar='(-)', required=False, default=0, type=float, help='Spectral index of pulsar')
        parser.add_argument('--duty_cycle', metavar='(phase)', required=False, default=0.1, type=float, help='Duty cycle of default gaussian pulse profile')
        parser.add_argument('--profile', metavar='(file/dict)', required=False, default='default', help='NumPy .npy or EPN .txt file containing a custom pulsar pulse profile (1D or 2D), or multi-component dictionary')
        parser.add_argument('--micro_structure', metavar='(microsec)', required=False, default=0, type=float, help='Mean timescale of pulse microstructure')
        parser.add_argument('--scattering_time', metavar='(millisec)', required=False, default=0, type=float, help='Scattering timescale due to ISM')
        parser.add_argument('--scattering_index', metavar='(-)', required=False, default=-4, type=float, help='Scattering index to describe frequency evolution')
        parser.add_argument('--DM_smear', metavar='(0 or 1)', required=False, default=0, type=int, help='Smear the pulse profile due to intra-channel DM smearing')

        parser.add_argument('--binary_period', metavar='(hour)', required=False, type=float, help='Period of binary oribit')
        parser.add_argument('--T0', metavar='(MJD)', required=False,  type=float, help='Reference epoch of pulsar periapsis (default: obseravtion start, barycentre)')
        parser.add_argument('--x', metavar='(light-sec)', required=False, type=float, help='Projected semi-major orbital axis')
        parser.add_argument('--M1', metavar='(M_sun)', required=False,  type=float, help='Mass of pulsar (default: 1.4)')
        parser.add_argument('--M2', metavar='(M_sun)', required=False, type=float, help='Companion mass')
        parser.add_argument('--inc', metavar='(deg)', required=False, type=float, help='Orbital inclination (default: 90)')
        parser.add_argument('--ecc', metavar='(-)', required=False, default=0, type=float, help='Orbital eccentricity')
        parser.add_argument('--AoP', metavar='(deg)', required=False,  default=0, type=float, help='Argument of periapsis')
        parser.add_argument('--LoAN', metavar='(deg)', required=False,  default=0, type=float, help='Longitude of the ascending node')
        parser.add_argument('--double_pulsar', metavar='(ID)', required=False, type=str, help='Make pulsar the companion of pulsar {double_pulsar:ID}')

        parser.add_argument('--mode', metavar='(str)', required=False, default='python', type=str, help="Inject using analytical 'python' code or polycos from 'pint'")
        parser.add_argument('--pint_N', metavar='(-)', required=False, default=12, type=int, help='Number of coefficients per timestep for polycos generation')
        parser.add_argument('--pint_T', metavar='(min)', required=False, default=5, type=float, help='Timestep for polycos generation')
        parser.add_argument('--polycos', metavar='(file)', required=False, type=str, help='.polycos: Pint will use this file, .par: Pint will make polycos from par file, (None) Pint will make from injection params.')

        self.parser = parser
    
    def dict_to_args(self, params):
        args_list = []
        for key, value in params.items():
            args_list.append(f'--{key}')
            args_list.append(str(value))
        return args_list
    
    @staticmethod
    def str2func(value, par, id, func):
        try:
            converted_value = func(value)
        except ValueError:
            sys.exit(f'Error: Invalid {par} for pulsar {id}')
        return converted_value
    
    def parse_inputs(self, pulsar_pars, global_pars):
        if global_pars:
            psr_pars = global_pars.copy()
            psr_pars.update(pulsar_pars)
        else:
            psr_pars = pulsar_pars

        FX_list, PX_list = self.get_spin_params(psr_pars)
        AX_list = self.get_accel(psr_pars)
        RAJ = psr_pars.pop('RAJ', None)
        DECJ = psr_pars.pop('DECJ', None)
        psr_args_list = self.dict_to_args(psr_pars)
        psr_args, _ = self.parser.parse_known_args(psr_args_list)

        clean_psr_pars = vars(psr_args)
        clean_psr_pars['RAJ'] = RAJ
        clean_psr_pars['DECJ'] = DECJ
        clean_psr_pars['FX'] = FX_list
        clean_psr_pars['PX'] = PX_list
        clean_psr_pars['AX'] = AX_list

        if type(clean_psr_pars['profile']) == str:
            if clean_psr_pars['profile'][0] == "{":
                clean_psr_pars['profile'] = eval(clean_psr_pars['profile'])

        clean_psr_pars = self.calc_binary_pars(clean_psr_pars)
        self.psr_pars = clean_psr_pars
    
    def get_accel(self, pulsar_pars):
        ID = pulsar_pars['ID']
        presto_z = self.str2func(pulsar_pars.get('presto_z', 0), 'presto_z', ID, float)
        presto_w = self.str2func(pulsar_pars.get('presto_w', 0), 'presto_w', ID, float)

        if presto_z or presto_w:
            accel_vals = ['presto', presto_z, presto_w]
        else:
            pattern = re.compile(r'^[A]\d+$')
            accel_keys = [key for key in pulsar_pars.keys() if pattern.match(key)]
            if accel_keys:
                accel_keys.sort(key=lambda x: int(x[1:]))
                accel_vals = [self.str2func(pulsar_pars.get(f'A{i}', 0), accel_keys[i], ID, float) for i in range(int(accel_keys[-1][1:])+1)]
            else:
                accel_vals = []
        return accel_vals
        
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
        if not pulsar_pars.get('SNR', 0):
            sys.exit(f'SNR value is required for pulsar {ID}.')

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
    
    @staticmethod
    def orbit_par_converter(period, find='A1', A1=None, M1=None, M2=None, inc=None):
        T = const.G.value * period**2 / (4*np.pi**2)

        def get_M1(M2_):
            sini = np.sin(np.deg2rad(inc))
            term1 = np.sqrt(T) * ((sini * M2_*const.M_sun.value)/(const.c.value * A1))**(3/2)
            return term1/const.M_sun.value - M2_

        if find == 'inc':
            sini3 = (M1 + M2)**2/M2**3 *1/const.M_sun.value * 1/T * (const.c.value * A1) ** 3
            return np.rad2deg(np.arcsin(sini3**(1/3)))
        
        elif find == 'M1':
            return get_M1(M2)
        
        elif find == 'M2':
            return fsolve(lambda M2_: M1 - get_M1(M2_), 1)
        
        elif find == 'A1':
            sini = np.sin(np.deg2rad(inc))
            return BinaryModel.get_semi_major(period, M1+M2) * M2/(M1+M2) * sini / const.c.value
        
    def calc_binary_pars(self, pulsar_pars):
        if not pulsar_pars['binary_period']:
            return pulsar_pars
        period = abs(pulsar_pars['binary_period']) * 3600

        A1 = pulsar_pars['x']
        M1 = pulsar_pars['M1']
        M2 = pulsar_pars['M2']
        inc = pulsar_pars['inc']
        M1_default, inc_default = 1.4, 90
        if A1:
            if not (M1 and M2 and inc):
                if (not M1) and (not M2) and (not inc):
                    M1 = M1_default
                    inc = inc_default
                    M2 = self.orbit_par_converter(period, find='M2', A1=A1, M1=M1_default, inc=inc_default) 
                elif (M1) and (not M2) and (not inc):
                    inc = inc_default
                    M2 = self.orbit_par_converter(period, find='M2', A1=A1, M1=M1, inc=inc_default)                
                elif (not M1) and (M2) and (not inc):
                    inc = inc_default
                    M1 = self.orbit_par_converter(period, find='M1', A1=A1, M2=M2, inc=inc_default)                
                elif (not M1) and (not M2) and (inc):
                    M1 = M1_default
                    M2 = self.orbit_par_converter(period, find='M2', A1=A1, M1=M1_default, inc=inc)                 
                elif (M1) and (M2) and (not inc):
                    inc = self.orbit_par_converter(period, find='inc', A1=A1, M1=M1, M2=M2)                 
                elif (M1) and (not M2) and (inc):
                    M2 = self.orbit_par_converter(period, find='M2', A1=A1, M1=M1, inc=inc)                 
                elif (not M1) and (M2) and (inc):
                    M1 = self.orbit_par_converter(period, find='M1', A1=A1, M2=M2, inc=inc)                
                else:
                    A1_calc = self.orbit_par_converter(period, find='A1', M1=M1, M2=M2, inc=inc) 
                    if A1_calc != A1:
                        sys.exit(f"Error: Inputed A1 for pulsar {pulsar_pars['ID']} is not compatible with inputed M1, M2 and inc. Only three out of these four parameters are required.")
        elif M1 and M2 and inc:
            A1 = self.orbit_par_converter(period, find='A1', M1=M1, M2=M2, inc=inc) 
        elif M1 and M2 and (not inc):
            inc = inc_default
            A1 = self.orbit_par_converter(period, find='A1', M1=M1, M2=M2, inc=inc) 
        elif (not M1) and M2 and inc:
            M1 = M1_default
            A1 = self.orbit_par_converter(period, find='A1', M1=M1, M2=M2, inc=inc) 
        elif (not M1) and M2 and (not inc):
            M1, inc = M1_default, inc_default
            A1 = self.orbit_par_converter(period, find='A1', M1=M1, M2=M2, inc=inc) 
        else:
            sys.exit(f"Error: Pulsar {pulsar_pars['ID']} requires either an A1 or M1, M2 and inc.")

        pulsar_pars['binary_period'] = period
        pulsar_pars['x'] = abs(A1)
        pulsar_pars['M1'] = abs(M1)
        if type(M2) == np.ndarray:
            M2 = M2[0]
        if type(inc) == np.ndarray:
            inc = inc[0]
        pulsar_pars['M2'] = abs(M2)
        pulsar_pars['inc'] = inc
        pulsar_pars['ecc'] = abs(pulsar_pars['ecc'])

        return pulsar_pars