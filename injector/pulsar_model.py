import numpy as np 
import astropy.units as u
from math import factorial
import astropy.constants as const
from sympy import lambdify, symbols
from scipy.interpolate import interp1d, RegularGridInterpolator

from .propagation_effects import PropagationEffects
from .pulsar_emission import PulsarEmission
from .micro_structure import MicroStructure


class PulsarModel:
    def __init__(self, obs, binary, pulsar_pars, generate=True):
        self.mode = pulsar_pars['mode'] if pulsar_pars['mode'] else 'python'

        self.ID = pulsar_pars['ID']
        self.seed = pulsar_pars['seed']
        self.SNR = pulsar_pars['SNR']
        self.pulsar_pars = pulsar_pars
        self.obs = obs
        self.binary = binary

        self.PX_list = pulsar_pars['PX']
        self.FX_list = pulsar_pars['FX']
        self.AX_list = pulsar_pars['AX']
        
        self.micro_structure = pulsar_pars['micro_structure']
        
        self.get_epochs()
        self.get_spin_functions(pulsar_pars)
        
        self.emission = PulsarEmission(obs, pulsar_pars)
        self.intrinsic_profile_chan = self.emission.get_intrinsic_profile()
        self.prop_effect = PropagationEffects(self.obs, pulsar_pars, self.emission.profile_length, self.period, self.emission.spectra)

        if generate:
            self.get_mode_generators(pulsar_pars, generate)

            self.observed_profile_chan = self.get_observed_profile()
            self.calculate_SNR()

            self.observed_profile = self.vectorise_observed_profile()

    def get_mode_generators(self, pulsar_pars, generate):
        if self.mode == 'python':
            if  pulsar_pars['frame'] == 'topo':
                self.generate_signal = self.generate_signal_python_topo
            else: 
                self.generate_signal = self.generate_signal_python_bary
        elif self.mode == 'pint':
            self.polycos_path = pulsar_pars['polycos']
            self.get_polyco_interp(generate)
            if  pulsar_pars['frame'] == 'topo':
                self.generate_signal = self.generate_signal_polcos_topo
            else:
                self.generate_signal = self.generate_signal_polcos_bary

    def get_observed_profile(self):
        smeared_profile = self.prop_effect.intra_channel_DM_smearing(self.intrinsic_profile_chan)
        scatterd_profile = self.prop_effect.ISM_scattering(smeared_profile)
        return scatterd_profile
    
    def get_epochs(self):
        self.binary.T0 = self.obs.T0
        self.posepoch = self.obs.posepoch
        self.pepoch = self.obs.pepoch
        self.spin_ref = self.obs.spin_ref
        self.orbit_ref = self.obs.orbit_ref
        self.accepoch = self.obs.accepoch

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
            self.phase_func = lambda t: phase_func_abs(t - self.accepoch, const.c.value) + phase_offset

        else:
            phase_symbolic = sum([FX[n]*t**(n+1)/factorial(n+1) for n in range(n_freq)])
            phase_func_abs = lambdify(t, phase_symbolic.subs(freq_derivs))
            self.phase_func = lambda t: phase_func_abs(t) + phase_offset

    def calculate_SNR(self):
        beam_scale = self.obs.get_beam_snr() 

        n_chan = self.obs.n_chan
        p0 = self.pulsar_pars.get('P0_SNR', self.period)
        n_pulse = self.obs.obs_len/p0

        nbins = self.emission.profile_length
        phase = np.linspace(0, 1, nbins)

        intrinsic_profile_sum = np.sum([self.intrinsic_profile_chan(phase, chan) for chan in range(n_chan)], axis=0) 
        profile_energy_scale = np.sum((intrinsic_profile_sum*n_pulse)**2)

        N_subints = 2**10
        dt_sub = self.obs.obs_len / N_subints
        sub_times = (np.arange(N_subints) + 0.5) * dt_sub
        LC = self.emission.light_curve(sub_times)
        LC /= np.mean(LC)
        profile_energy_scale *= np.mean(LC**2)

        samples_per_bin =  nbins / (p0 / self.obs.dt)
        noise_energy = self.obs.fb_std ** 2 * (n_pulse * n_chan) * samples_per_bin
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
    
    def get_polyco_interp(self, generate_range):
        from pint.polycos import Polycos # type: ignore
        polycos_model = Polycos.read(self.polycos_path)
        interp_topo_mjd = self.obs.observation_span(generate_range, n_samples=10**5, return_mjd=True)
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
    
    def generate_signal_polcos_bary(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        DM_array = np.tile(self.prop_effect.DM_delays, (len(timeseries),1))

        topo_times = self.obs.sec2mjd(timeseries)
        phase_array = np.tile(topo_times, (len(self.obs.freq_arr),1)).T
        phase_time = (phase_array + DM_array*u.s.to(u.day))

        bary_times = self.obs.topo2bary(timeseries, return_mjd=False, interp=True)
        bary_array = np.tile(bary_times, (len(self.obs.freq_arr),1)).T

        phase = self.polycos(phase_time) + self.get_phase(bary_array + DM_array)
        LC = self.emission.light_curve(timeseries)[:, None]
        return self.get_pulse(phase, freq_array) * LC
    
    def generate_signal_polcos_topo(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        DM_array = np.tile(self.prop_effect.DM_delays, (len(timeseries),1))

        topo_times = self.obs.sec2mjd(timeseries)
        phase_array = np.tile(topo_times, (len(self.obs.freq_arr),1)).T
        phase_time = (phase_array + DM_array*u.s.to(u.day))
        topo_array = np.tile(timeseries, (len(self.obs.freq_arr),1)).T

        phase = self.polycos(phase_time) + self.get_phase(topo_array + DM_array)
        LC = self.emission.light_curve(timeseries)[:, None]
        return self.get_pulse(phase, freq_array) * LC
        
    def generate_signal_python_bary(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        DM_array = np.tile(self.prop_effect.DM_delays, (len(timeseries),1))
        obs_freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))

        bary_times = self.obs.topo2bary(timeseries, return_mjd=False, interp=True)
        bary_array = np.tile(bary_times, (len(self.obs.freq_arr),1)).T

        phase_array = self.get_phase(bary_array + DM_array)
        LC = self.emission.light_curve(timeseries)[:, None]
        return self.get_pulse(phase_array, obs_freq_array) * LC
    
    def generate_signal_python_topo(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        DM_array = np.tile(self.prop_effect.DM_delays, (len(timeseries),1))
        obs_freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        topo_array = np.tile(timeseries, (len(self.obs.freq_arr),1)).T

        phase_array = self.get_phase(topo_array + DM_array)
        LC = self.emission.light_curve(timeseries)[:, None]
        return self.get_pulse(phase_array, obs_freq_array) * LC

   