import sys
import numpy as np 
import pandas as pd
from pathlib import Path
import astropy.units as u
from math import factorial
from scipy.stats import norm
import astropy.constants as const
from sympy import lambdify, symbols
from scipy.interpolate import interp1d, RegularGridInterpolator

from .propagation_effects import PropagationEffects
from .micro_structure import MicroStructure


class PulsarModel:
    def __init__(self, obs, binary, pulsar_pars, generate=True):
        self.get_mode(pulsar_pars)
        self.ID = pulsar_pars['ID']
        self.seed = pulsar_pars['seed']
        self.SNR = pulsar_pars['SNR']
        self.pulsar_pars = pulsar_pars
        self.obs = obs
        self.binary = binary
        
        self.get_epochs(pulsar_pars)
        self.PX_list = pulsar_pars['PX']
        self.FX_list = pulsar_pars['FX']
        self.AX_list = pulsar_pars['AX']
        self.get_spin_functions(pulsar_pars)

        self.spectra = self.get_spectra(pulsar_pars)
        self.intrinsic_profile_chan = self.get_intrinsic_profile(pulsar_pars)
        self.micro_structure = pulsar_pars['micro_structure']
        self.prop_effect = PropagationEffects(self.obs, pulsar_pars, self.profile_length, self.period, self.spectra)

        if generate:
            self.observed_profile_chan = self.get_observed_profile()
            self.calculate_SNR()
            self.observed_profile = self.vectorise_observed_profile()


    def get_mode(self, pulsar_pars):
        self.mode = pulsar_pars['mode'] if pulsar_pars['mode'] else 'python'
        self.polycos_path = pulsar_pars['polycos']
        if self.mode == 'python':
            self.generate_signal = self.generate_signal_python
        elif self.mode == 'pint':
            self.generate_signal = self.generate_signal_polcos

    def get_observed_profile(self):
        scatterd_profile = self.prop_effect.ISM_scattering(self.intrinsic_profile_chan)
        smeared_profile = self.prop_effect.intra_channel_DM_smearing(scatterd_profile)
        return smeared_profile

    def get_epochs(self, pulsar_pars):
        pepoch = pulsar_pars['PEPOCH'] if pulsar_pars['PEPOCH'] else self.obs.obs_start_bary
        # posepoch = pulsar_pars['POSEPOCH'] if pulsar_pars['POSEPOCH'] else pepoch
        T0 = pulsar_pars['T0'] if pulsar_pars['T0'] else self.obs.obs_start_bary

        spin_ref = (self.obs.obs_start_bary - pepoch) * u.day.to(u.s)
        # pos_ref = (self.obs.obs_start_bary - posepoch) * u.day.to(u.s)
        orbit_ref = (self.obs.obs_start_bary - T0) * u.day.to(u.s)

        self.pepoch = pepoch
        # self.posepoch = posepoch
        self.binary.T0 = T0
        self.spin_ref = spin_ref
        # self.pos_ref = pos_ref
        self.orbit_ref = orbit_ref
    
    def get_spin_functions(self, pulsar_pars):
        t, c = symbols('t, c')
        phase_offset = pulsar_pars['phase_offset']
        n_freq, n_accel = len(self.FX_list), len(self.AX_list)

        FX = symbols([f'F{x}' for x in range(n_freq)])
        freq_derivs = dict(zip(FX, self.FX_list))

        spin_symbolic = sum([FX[n]*t**n/factorial(n) for n in range(n_freq)])
        self.spin_func = lambdify(t, spin_symbolic.subs(freq_derivs))
        self.period = 1/self.spin_func(self.spin_ref)

        if n_accel:
            AX = symbols([f'A{x}' for x in range(n_accel)])
            accel_derivs = dict(zip(AX, self.AX_list))

            Vel_symbolic = sum([AX[n]*t**(n+1)/factorial(n+1) for n in range(n_accel)])
            spin_doppler = spin_symbolic * (1 - Vel_symbolic/c)
            phase_symbolic = spin_doppler.integrate(t)  
            phase_func_abs = lambdify([t, c], phase_symbolic.subs({**freq_derivs, **accel_derivs}))
            self.phase_func = lambda t: phase_func_abs(t, const.c.value) + phase_offset

        else:
            phase_symbolic = sum([FX[n]*t**(n+1)/factorial(n+1) for n in range(n_freq)])
            phase_func_abs = lambdify(t, phase_symbolic.subs(freq_derivs))
            self.phase_func = lambda t: phase_func_abs(t) + phase_offset

    def get_spectra(self, pulsar_pars):            
        PSD_file = pulsar_pars['PSD']
        if PSD_file:
            suffix = Path(PSD_file).suffix
            if suffix == '.npy':
                try:
                    spectra_arr = np.load(PSD_file)
                except FileNotFoundError:
                    sys.exit(f'Unable to load {PSD_file} numpy spectra for pulsar {self.ID}.')
            else:
                sys.exit(f'Pulsar {self.ID} has an invalid spectra file extension: {PSD_file}. Must be a numpy .npy file.')

            freq_min = np.min(self.obs.freq_arr) - abs(self.obs.df)/2
            freq_max = np.max(self.obs.freq_arr) + abs(self.obs.df)/2
            freq_range = np.linspace(freq_min, freq_max, len(spectra_arr))
            spectra_interp = interp1d(freq_range, spectra_arr)
            spectra_norm = spectra_interp(self.obs.f0)
            return lambda freq: spectra_interp(freq)/spectra_norm
        else:
            spectral_index = pulsar_pars['spectral_index']
            return lambda freq: (freq/self.obs.f0)**float(spectral_index)
    
    @staticmethod
    def parse_profile(profile, pulse_i):
        dc = profile['duty_cycle'][pulse_i]
        phase = profile['phase'][pulse_i]
        amp =  profile['amp'][pulse_i]

        phase_range = np.linspace(0, 1, 1000)
        pulse_sigma = (dc)/(2*np.sqrt(2*np.log(2)))
        pulse = norm(phase, pulse_sigma).pdf(phase_range)
        pulse /= np.max(pulse)
        return pulse * amp
    
    def get_intrinsic_profile(self, pulsar_pars):
        profile = pulsar_pars['profile']
        if profile == 'default':
            duty_cycle = pulsar_pars['duty_cycle']
            pulse_sigma = (duty_cycle)/(2*np.sqrt(2*np.log(2)))
            self.profile_length = 1000
            def intrinsic_pulse(phase, chan_num=0): 
                return np.exp(-(phase-0.5)**2/(2*(pulse_sigma)**2)) * self.spectra(self.obs.freq_arr[chan_num])
            
        else:
            if type(profile) == str:
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

            elif type(profile) == dict:
                profile_arr = np.zeros(1000)
                for pulse_i in range(len(profile['phase'])):
                    profile_arr += self.parse_profile(profile, pulse_i)

            else:
                sys.exit(f'Invalid pulse profile for pulsar {self.ID}.')

            if profile_arr.ndim == 1:
                phase_range = np.linspace(0, 1, len(profile_arr))
                self.interp_func = interp1d(phase_range, profile_arr/np.max(profile_arr))
                self.profile_length = len(profile_arr)
                def intrinsic_pulse(phase, chan_num=0): 
                    return self.interp_func(phase) * self.spectra(self.obs.freq_arr[chan_num]) 
                
            elif (profile_arr.dim == 2) and (len(profile_arr) == self.obs.n_chan):
                phase_range = np.linspace(0, 1, len(profile_arr.T))
                self.interp_funcs = [interp1d(phase_range, profile_arr[i]/np.max(profile_arr)) for i in range(self.obs.n_chan)]
                self.profile_length = len(profile_arr.T)
                def intrinsic_pulse(phase, chan_num): 
                    return self.interp_funcs[chan_num](phase) * self.spectra(self.obs.freq_arr[chan_num])
                
            else:
                sys.exit(f'Unable to interpret {profile} numpy pulse profile for pulsar {self.ID}.')

        return intrinsic_pulse

    def calculate_SNR(self):

        beam_scale = self.obs.get_beam_snr() 

        n_chan = self.obs.n_chan
        p0 = self.PX_list[0]
        n_pulse = self.obs.obs_len/p0

        nbins = int(np.round(p0/self.obs.dt))
        phase = np.linspace(0, 1, nbins)
        intrinsic_profile_sum = np.sum([self.intrinsic_profile_chan(phase, chan) for chan in range(n_chan)], axis=0) 
        profile_energy_scale = np.sum((intrinsic_profile_sum*n_pulse)**2)
        noise_energy = self.obs.fb_std ** 2 * (n_pulse * n_chan)
        snr = profile_energy_scale/noise_energy

        self.SNR_scale = self.SNR / np.sqrt(snr) * beam_scale
     
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
            pulse_generator = MicroStructure(phase_abs, freq, self.micro_structure, self.period, self.observed_profile, self.seed)
            return pulse_generator.pulse_profile()
        else:
            return self.observed_profile(phase_abs % 1, freq)
    
    def coord2proper_time(self, bary_times):
        return bary_times+self.spin_ref - self.binary.orbital_delay(bary_times+self.orbit_ref)
    
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
    
    def generate_signal_python_topo(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        DM_array = np.tile(self.prop_effect.DM_delays, (len(timeseries),1))
        obs_freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        bary_array = np.tile(timeseries, (len(self.obs.freq_arr),1)).T

        phase_array = self.get_phase(bary_array + DM_array)
        return self.get_pulse(phase_array, obs_freq_array)

   