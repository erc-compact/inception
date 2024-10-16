import re
import sys
import warnings
import numpy as np 
import pandas as pd
from pathlib import Path
import astropy.units as u
from math import factorial
from sympy import lambdify, symbols, Function
from scipy.integrate import quad, IntegrationWarning
from scipy.interpolate import interp1d, RegularGridInterpolator

from .propagation_effects import PropagationEffects
from .micro_structure import MicroStructure
from .io_tools import str2func


class PulsarModel:
    def __init__(self, obs, binary, pulsar_pars, generate=True):
        self.get_mode(pulsar_pars)
        self.ID = pulsar_pars['ID']
        self.obs = obs
        self.binary = binary
        
        self.get_epochs(pulsar_pars)
        self.get_spin_params(pulsar_pars)
        self.get_spin_functions(pulsar_pars)
        self.get_period_start()

        self.spectra = self.get_spectra(pulsar_pars)
        self.intrinsic_profile_chan = self.get_intrinsic_profile(pulsar_pars)
        self.micro_structure = str2func(pulsar_pars.get('micro_structure', 0), 'micro_structure', self.ID, float)
        self.prop_effect = PropagationEffects(self.obs, pulsar_pars, self.profile_length, self.period, self.spectra)

        if generate:
            self.observed_profile_chan = self.get_observed_profile()
            self.observed_profile = self.vectorise_observed_profile()

        self.calculate_SNR(pulsar_pars, generate)
    
    def get_mode(self, pulsar_pars):
        self.mode = pulsar_pars.get('mode', 'python')
        self.polycos_path = pulsar_pars.get('polycos', None)
        if self.mode == 'python':
            self.generate_signal = self.generate_signal_python
        elif self.mode == 'pint':
            self.generate_signal = self.generate_signal_polcos

    def get_observed_profile(self):
        scatterd_profile = self.prop_effect.ISM_scattering(self.intrinsic_profile_chan)
        smeared_profile = self.prop_effect.intra_channel_DM_smearing(scatterd_profile)
        return smeared_profile

    def get_epochs(self, pulsar_pars):
        pepoch = pulsar_pars.get('PEPOCH', self.obs.obs_start_bary)
        posepoch = pulsar_pars.get('POSEPOCH', pepoch)
        T0 = pulsar_pars.get('T0', self.obs.obs_start_bary)

        pepoch_float = str2func(pepoch, 'PEPOCH', self.ID, float) 
        posepoch_float = str2func(posepoch, 'POSPOCH', self.ID, float) 
        T0_float = str2func(T0, 'T0', self.ID, float) 

        spin_ref = (self.obs.obs_start_bary - pepoch_float) * u.day.to(u.s)
        pos_ref = (self.obs.obs_start_bary - posepoch_float) * u.day.to(u.s)
        orbit_ref = (self.obs.obs_start_bary - T0_float) * u.day.to(u.s)

        self.pepoch = pepoch_float
        self.posepoch = posepoch_float
        self.binary.T0 = T0_float
        self.spin_ref = spin_ref
        self.pos_ref = pos_ref
        self.orbit_ref = orbit_ref
    
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
        p0_find = pulsar_pars.get('P0', 0) 
        f0_find = pulsar_pars.get('F0', 0)
        if (not p0_find) and (not f0_find):
            sys.exit(f'No P0 or F0 found for pulsar {self.ID}. Pulsar must be spinning!')

        pattern = re.compile(r'^[PF]\d+$')
        matching_keys = [key for key in pulsar_pars.keys() if pattern.match(key)]
        matching_keys.sort(key=lambda x: int(x[1:]))

        spin_types = []
        for i in range(int(matching_keys[-1][1:])+1):
            if (f'P{i}' in matching_keys) and (f'F{i}' not in matching_keys):
                spin_types.append(f'P{i}')
            else:
                spin_types.append(f'F{i}')
        spin_values = [str2func(pulsar_pars.get(key, 0), key, self.ID, float) for key in spin_types]

        self.FX_list, self.PX_list  = self.spin_pars_converter(spin_values, spin_types)

    def get_spin_functions(self, pulsar_pars):
        t = symbols('t')
        FX = symbols([f'F{x}' for x in range(len(self.FX_list))])
        freq_derivs = dict(zip(FX, self.FX_list))

        phase_symbolic = sum([FX[n]*t**(n+1)/factorial(n+1) for n in range(len(self.FX_list))])
        phase_offset = str2func(pulsar_pars.get('phase_offset', 0), 'phase_offset', self.ID, float)
        phase_func_abs = lambdify(t, phase_symbolic.subs(freq_derivs))
        self.phase_func = lambda t: phase_func_abs(t) + phase_offset

        spin_symbolic = sum([FX[n]*t**(n)/factorial(n) for n in range(len(self.FX_list))])
        self.spin_func = lambdify(t, spin_symbolic.subs(freq_derivs))
    
    def coord2proper_time(self, bary_times):
        return bary_times+self.spin_ref - self.binary.orbital_delay(bary_times+self.orbit_ref)
    
    def get_period_start(self):
        T_mid_proper = self.obs.obs_start_bary + self.spin_ref 
        self.period = 1/self.spin_func(T_mid_proper)

    def get_spectra(self, pulsar_pars):
        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                return False
            
        spectra_input = pulsar_pars.get('spectra', 0)
        if is_number(spectra_input):
            return lambda freq: (freq/self.obs.f0)**float(spectra_input)
        else:
            suffix = Path(spectra_input).suffix
            if suffix == '.npy':
                try:
                    spectra_arr = np.load(spectra_input)
                except FileNotFoundError:
                    sys.exit(f'Unable to load {spectra_input} numpy spectra for pulsar {self.ID}.')
            else:
                sys.exit(f'Pulsar {self.ID} has an invalid spectra file extension: {spectra_input}. Must be a numpy .npy file.')

            freq_min = np.min(self.obs.freq_arr) - abs(self.obs.df)/2
            freq_max = np.max(self.obs.freq_arr) + abs(self.obs.df)/2
            freq_range = np.linspace(freq_min, freq_max, len(spectra_arr))
            spectra_interp = interp1d(freq_range, spectra_arr)
            spectra_norm = spectra_interp(self.obs.f0)
            return lambda freq: spectra_interp(freq)/spectra_norm
    
    def get_intrinsic_profile(self, pulsar_pars):
        profile = pulsar_pars.get('profile', 'default')
        if profile == 'default':
            duty_cycle = str2func(pulsar_pars.get('duty_cycle', 0.1), 'duty_cycle', self.ID, float)
            pulse_sigma = (duty_cycle)/(2*np.sqrt(2*np.log(2)))
            self.profile_length = 1000
            def intrinsic_pulse(phase, chan_num=0): 
                return np.exp(-(phase-0.5)**2/(2*(pulse_sigma)**2)) * self.spectra(self.obs.freq_arr[chan_num])
            
        else:
            suffix = Path(profile).suffix
            if suffix == '.npy':
                try:
                    profile_arr = np.load(profile)
                except FileNotFoundError:
                    sys.exit(f'Unable to load {profile} numpy pulse profile for pulsar {self.ID}.')
            elif suffix == '.txt':
                try:
                    epn_profile = pd.read_csv(profile, delimiter=' ', 
                                            names=['col0', 'col1', 'col2', 'intensity'])
                    profile_arr = epn_profile['intensity'].values
                except (FileNotFoundError, KeyError, ValueError):
                    sys.exit(f'Unable to load {profile} EPN pulse profile for pulsar {self.ID}.')
            else:
                sys.exit(f'Pulsar {self.ID} has an invalid profile file extension: {profile}. Must be a numpy .npy or EPN .txt file.')

            if profile_arr.ndim == 1:
                phase_range = np.linspace(0, 1, len(profile_arr))
                self.interp_func = interp1d(phase_range, profile_arr/np.max(profile_arr))
                self.profile_length = len(profile_arr)
                def intrinsic_pulse(phase, chan_num=0): 
                    return self.interp_func(phase) * self.spectra(self.obs.freq_arr[chan_num]) 
                
            elif (profile_arr.dim == 2) and (len(profile_arr) == self.nchans):
                phase_range = np.linspace(0, 1, len(profile_arr.T))
                self.interp_funcs = [interp1d(phase_range, profile_arr[i]/np.max(profile_arr)) for i in range(self.nchans)]
                self.profile_length = len(profile_arr.T)
                def intrinsic_pulse(phase, chan_num): 
                    return self.interp_funcs[chan_num](phase) * self.spectra(self.obs.freq_arr[chan_num])
                
            else:
                sys.exit(f'Unable to interpret {profile} numpy pulse profile for pulsar {self.ID}.')

        return intrinsic_pulse
    
    def calculate_SNR(self, pulsar_pars, generate):
        SNR_obs = str2func(pulsar_pars.get('SNR', 0), 'SNR', self.ID, float)
        SNR_intrinsic = str2func(pulsar_pars.get('SNRi', 0), 'SNRi', self.ID, float)
        if (not SNR_obs) and (not SNR_intrinsic):
            sys.exit(f'SNR value is required for pulsar {self.ID}.')

        if generate:
            if SNR_obs:
                SNR = SNR_obs
                SNR_function = self.observed_profile_chan
            elif SNR_intrinsic:
                SNR = SNR_intrinsic
                SNR_function = self.intrinsic_profile_chan
        
            integrated_profile = 0
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", IntegrationWarning)
                for chan in range(self.obs.n_chan):
                    integrated_profile += quad(SNR_function, args=(chan,), a=0, b=1, epsabs=1e-5)[0]

            self.SNR_scale = SNR * self.obs.fb_std * np.sqrt(self.period * self.obs.n_chan / self.obs.obs_len) / integrated_profile
    
    def vectorise_observed_profile(self):
        phases = self.prop_effect.phase
        freqs = self.obs.freq_arr
        pulse_arr  = np.vstack([self.observed_profile_chan(phases, chan) for chan in range(self.obs.n_chan)]).T
        self.grid_interp = RegularGridInterpolator((phases, freqs), pulse_arr)

        def observed_profile_function(phase, freq):
            grid = np.array([phase.ravel(), freq.ravel()]).T 
            return self.grid_interp(grid).reshape(phase.shape) * self.SNR_scale
        
        return observed_profile_function
    
    def get_polyco_interp(self):
        from pint.polycos import Polycos # type: ignore
        polycos_model = Polycos.read(self.polycos_path)
        interp_topo_mjd = self.obs.observation_span(n_samples=10**6)
        abs_phase_interp = polycos_model.eval_abs_phase(interp_topo_mjd).value

        self.polycos = interp1d(interp_topo_mjd.astype(np.float64), abs_phase_interp.astype(np.float64))
    
    def get_pulse(self, phase_abs, freq):
        if self.micro_structure:
            pulse_generator = MicroStructure(phase_abs, freq, self.micro_structure, self.period, self.observed_profile)
            return pulse_generator.pulse_profile()
        else:
            return self.observed_profile(phase_abs % 1, freq)
        
    def get_phase(self, bary_times):
        T_proper = self.coord2proper_time(bary_times)
        phase_abs = self.phase_func(T_proper)
        return phase_abs 
    
    def generate_signal_polcos(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        DM_array = np.tile(self.prop_effect.DM_delays, (len(timeseries),1))

        topo_times = self.obs.sec2mjd(timeseries)
        phase_array = np.tile(topo_times, (len(self.obs.freq_arr),1)).T
        phase_time = (phase_array + DM_array*u.s.to(u.day))
        
        phase = self.polycos(phase_time)
        return self.get_pulse(phase, freq_array)
        
    def generate_signal_python(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        DM_array = np.tile(self.prop_effect.DM_delays, (len(timeseries),1))
        obs_freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))

        topo_times = self.obs.sec2mjd(timeseries)
        bary_times = self.obs.topo2bary(topo_times, mjd=False, interp=True)
        bary_array = np.tile(bary_times, (len(self.obs.freq_arr),1)).T

        phase_array = self.get_phase(bary_array + DM_array)
        return self.get_pulse(phase_array, obs_freq_array)

   