import re
import sys
import numpy as np 
from pathlib import Path
import astropy.units as u
from math import factorial
import astropy.constants as const
from scipy.interpolate import interp1d
from sympy import lambdify, symbols, Function

from injector.io_tools import str2func



class PulsarSignal:
    def __init__(self, obs, binary, pulsar_pars, validate=False, par_file=None):

        self.ID = pulsar_pars['ID']
        self.obs = obs
        self.binary = binary
        
        pepoch, posepoch, T0, spin_ref, pos_ref, orbit_ref = self.get_epochs(pulsar_pars)
        self.pepoch = pepoch
        self.posepoch = posepoch
        self.binary.T0 = T0
        self.spin_ref = spin_ref
        self.pos_ref = pos_ref
        self.orbit_ref = orbit_ref

        self.scattering_time = pulsar_pars.get('scattering_time', 'inf')
        self.duty_cycle = str2func(pulsar_pars.get('duty_cycle', 0.1), 'duty_cycle', self.ID, float)
        self.spectral_index = str2func(pulsar_pars.get('spectral_index', 0), 'spectral_index', self.ID, float)
        self.SNR = str2func(pulsar_pars.get('SNR', 100), 'SNR', self.ID, float)
        self.power_constant = self.SNR2power()/(self.obs.f0**self.spectral_index)

        self.FX_list, self.PX_list = self.get_spin_pars(pulsar_pars)
        self.phase_offset = str2func(pulsar_pars.get('phase_offset', 0), 'phase_offset', self.ID, float)

        self.pulse_profile = self.get_profile_func(pulsar_pars)
        self.spin_function = self.get_phase_function()

        self.mode = pulsar_pars.get('mode', 'python')
        self.polycos_path = pulsar_pars.get('polycos', None)
        self.polycos_coeff = pulsar_pars.get('pint_N', 12)
        self.polycos_tspan = pulsar_pars.get('pint_T', 5)
        
        self.polycos_reader = None
        self.polycos = None
        self.generate_signal = None

        if not validate:
            self.mode_resolver(par_file)


    def get_epochs(self, pulsar_pars):
        pepoch = pulsar_pars.get('PEPOCH', self.obs.obs_start_bary)
        posepoch = pulsar_pars.get('POSEPOCH', pepoch)
        T0 = pulsar_pars.get('T0', pepoch)

        pepoch_float = str2func(pepoch, 'PEPOCH', self.ID, float) 
        posepoch_float = str2func(posepoch, 'POSPOCH', self.ID, float) 
        T0_float = str2func(T0, 'T0', self.ID, float) 

        spin_ref = (self.obs.obs_start_bary - pepoch_float) * u.day.to(u.s)
        pos_ref = (self.obs.obs_start_bary - posepoch_float) * u.day.to(u.s)
        orbit_ref = (self.obs.obs_start_bary - T0_float) * u.day.to(u.s)
        return pepoch_float, posepoch_float, T0_float, spin_ref, pos_ref, orbit_ref
    
    
    def get_spin_pars(self, pulsar_pars):
        p0_find = pulsar_pars.get('P0', 0) 
        f0_find = pulsar_pars.get('F0', 0)
        if (not p0_find) and (not f0_find):
            raise Exception("InputError: No P0 or F0 found. Pulsar must be spinning!")

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

        FX_list, PX_list = spin_pars_converter(spin_values, spin_types)
        return FX_list, PX_list
    
    def get_profile_func(self, pulsar_pars):
        profile = pulsar_pars.get('profile', None)
        if profile:
            try:
                profile_arr = np.load(profile)
            except FileNotFoundError:
                sys.exit(f'Unable to load {profile} pulse profile for pulsar {self.ID}.')

            phase_range = np.linspace(0, 1, len(profile_arr))
            interp_func = interp1d(phase_range, profile_arr/np.max(profile_arr))

            def intrinsic_pulse(phase): return interp_func(phase)   
        else:
            pulse_sigma = (self.duty_cycle)/(2*np.sqrt(2*np.log(2)))
            def intrinsic_pulse(phase): return  np.exp(-(phase-0.5)**2/(2*(pulse_sigma)**2))
        
        if (self.scattering_time == 'inf') or (self.scattering_time == '0'):
            return lambda phase, freq: self.spectral_power(freq) * intrinsic_pulse(phase)   
        else:
            scattering_time = str2func(self.scattering_time, 'scattering_time', self.ID, float)
            profile_dt = 1000
            
            def scattering_kernal():
                kernal = np.exp(-scattering_time * np.linspace(0, 1, profile_dt*5))
                return kernal / np.sum(kernal)
            
            phase_range = np.linspace(0, 1, profile_dt)
            sample_profile = intrinsic_pulse(phase_range)
            tile_profile = np.tile(sample_profile, 5)
            
            convolve_profile = np.convolve(scattering_kernal(), tile_profile, mode='same')
            resample_profile = convolve_profile[profile_dt//2: (profile_dt*3)//2]
            interp_convolve = interp1d(phase_range.astype(np.float64), resample_profile.astype(np.float64))

            def observed_pulse(phase, freq): 
                return self.spectral_power(freq) * interp_convolve(phase)   
            
            return observed_pulse
        
    def polycos_creator(self, par_file, pint_func): # file is made xcpu 
        models, Polycos = pint_func
        timing_model = models.get_model(par_file, EPHEM=self.obs.ephem)

        t_mid = self.obs.obs_start + self.obs.obs_len/2 * u.s.to(u.day)
        polco_range = self.obs.obs_len/2 + 10*u.min.to(u.s)
        start, end = t_mid - polco_range*u.s.to(u.day), t_mid + polco_range*u.s.to(u.day)
        
        polycos_coeff = max(1, self.polycos_coeff)
        polycos_tspan = max(1, self.polycos_tspan) # minutes
        gen_poly = Polycos.generate_polycos(timing_model, start, end, self.obs.tempo_id, 
                                            polycos_tspan, polycos_coeff, 
                                            self.obs.f0, progress=False)
        
        polycos_path = Path(par_file).with_suffix('.polycos')
        self.polycos_path = polycos_path
        gen_poly.write_polyco_file(polycos_path)

    def mode_resolver(self, par_file):
        if self.mode not in ['pint', 'python']:
            self.mode = 'python' 

        polycos_path = self.polycos_path
        if polycos_path:
            self.mode == 'pint'

        if (self.mode == 'pint'):
            try:
                import pint.logging as logging
                _ = logging.setup('ERROR')  
                import pint.models as models
                from pint.polycos import Polycos
            except ImportError:
                sys.exit('pint-pulsar package not installed, cannot use polycos.')

        if (self.mode == 'pint') and not polycos_path:
            self.polycos_creator(par_file, pint_func=[models, Polycos])

        if (self.mode == 'pint'):
            self.polycos_reader = Polycos
            self.get_polyco_interp()

        self.get_generator()
        
    def get_polyco_interp(self):
        polycos_model = self.polycos_reader.read(self.polycos_path)
        interp_topo_mjd = self.obs.observation_span(n_samples=10**6)
        abs_phase_interp = polycos_model.eval_abs_phase(interp_topo_mjd).value

        self.polycos = interp1d(interp_topo_mjd.astype(np.float64), abs_phase_interp.astype(np.float64))

    def SNR2power(self):
        sigma_pt = self.obs.fb_std
        n_sample = self.obs.n_samples
        Weq_t = self.duty_cycle/(2*np.sqrt(2*np.log(2)))*np.sqrt(2*np.pi)
        Amp = self.SNR * sigma_pt / (np.sqrt(Weq_t) * np.sqrt(self.obs.n_chan) * np.sqrt(n_sample))
        
        return Amp * self.obs.get_beam_snr()

    def spectral_power(self, freq):
        return self.power_constant * (freq)**self.spectral_index
    
    def period_doppler_shift(self, vr, period, rel=True):
        beta =  vr/const.c.value
        if rel:
            return period * np.sqrt((1 + beta) / (1 - beta))
        else:
            return period * (1 + beta)
    
    def get_phase_function(self):
        t = symbols('t')
        FX = symbols([f'F{x}' for x in range(len(self.FX_list))])

        phase_symbolic = sum([FX[n]*t**(n+1)/factorial(n+1) for n in range(len(self.FX_list))])
        freq_derivs = dict(zip(FX, self.FX_list))

        return lambdify(t, phase_symbolic.subs(freq_derivs))

    def get_phase(self, bary_times):
        T_proper = bary_times+self.spin_ref - self.binary.orbital_delay(bary_times+self.orbit_ref)

        phase_abs = self.phase_offset + self.spin_function(T_proper)
        return phase_abs % 1 
    
    def generate_signal_polcos(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        DM_array = np.tile(self.obs.DM_delays, (len(timeseries),1))

        topo_times = self.obs.sec2mjd(timeseries)
        phase_array = np.tile(topo_times, (len(self.obs.freq_arr),1)).T
        phase_time = (phase_array + DM_array*u.s.to(u.day))
        
        phase = self.polycos(phase_time) % 1
        return self.pulse_profile(phase, freq_array)
        
    def generate_signal_python(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        DM_array = np.tile(self.obs.DM_delays, (len(timeseries),1))
        obs_freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        
        topo_times = self.obs.sec2mjd(timeseries)
        bary_times = self.obs.topo2bary(topo_times, mjd=False, interp=True)

        bary_array = np.tile(bary_times, (len(self.obs.freq_arr),1)).T
        phase_array = self.get_phase(bary_array + DM_array)

        return self.pulse_profile(phase_array, obs_freq_array)

    def get_generator(self):        
        if self.mode == 'pint':
            self.generate_signal = self.generate_signal_polcos
        elif self.mode == 'python':
            self.generate_signal = self.generate_signal_python
        

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