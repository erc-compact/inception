import numpy as np 
import astropy.units as u
from math import factorial
import astropy.constants as const
from sympy import lambdify, symbols, diff
from scipy.interpolate import interp1d

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
            self.get_mode_generators(pulsar_pars)

            self.observed_profile_chan = self.get_observed_profile()
            self.calculate_SNR()

            self.observed_profile = self.vectorise_observed_profile()

    def get_mode_generators(self, pulsar_pars):
        if self.mode == 'python':
            if  pulsar_pars['frame'] == 'topo':
                self.generate_signal = self.generate_signal_python_topo
            else: 
                self.generate_signal = self.generate_signal_python_bary
        elif self.mode == 'pint':
            self.polycos_path = pulsar_pars['polycos']
            self.get_polyco_interp()
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

            doppler = spin_doppler.subs({**freq_derivs, **accel_derivs, c: const.c.value})
            self.FX_doppler = [float(diff(doppler, t, k).subs(t, 0)) for k in range(n_freq + n_accel)]

        else:
            phase_symbolic = sum([FX[n]*t**(n+1)/factorial(n+1) for n in range(n_freq)])
            phase_func_abs = lambdify(t, phase_symbolic.subs(freq_derivs))
            self.phase_func = lambda t: phase_func_abs(t) + phase_offset
            self.FX_doppler = list(self.FX_list)

    def calculate_SNR(self):
        n_chan = self.obs.n_chan
        p0 = self.pulsar_pars.get('P0_SNR') or self.period
        n_pulse = self.obs.obs_len/p0

        nbins = self.emission.profile_length
        phase = np.linspace(0, 1, nbins)

        intrinsic_profile_sum = np.sum([self.intrinsic_profile_chan(phase, chan) for chan in range(n_chan)], axis=0) 
        if self.pulsar_pars.get('SNR_def', 'total') == 'pulsed':
            intrinsic_profile_sum -= np.mean(intrinsic_profile_sum)
        profile_energy_scale = np.sum((intrinsic_profile_sum*n_pulse)**2)

        samples_per_bin =  nbins / (p0 / self.obs.dt)
        noise_energy = self.obs.fb_std ** 2 * (n_pulse * n_chan) * samples_per_bin
        snr = profile_energy_scale/noise_energy

        self.SNR_scale = self.SNR / np.sqrt(snr) * self.obs.beam_scale
     
    def vectorise_observed_profile(self):
        phases = self.prop_effect.phase
        n_phase = len(phases)
        table = np.vstack([self.observed_profile_chan(phases, chan) for chan in range(self.obs.n_chan)]).ravel() * self.SNR_scale

        def observed_profile_function(phase, chan):
            x = phase * (n_phase - 1)
            i0 = np.minimum(x.astype(np.int64), n_phase - 2)
            frac = x - i0
            base = chan * n_phase + i0
            return table[base] * (1 - frac) + table[base + 1] * frac

        return observed_profile_function
    
    def get_polyco_interp(self):
        from pint.polycos import Polycos # type: ignore
        polycos_model = Polycos.read(self.polycos_path)
        interp_topo_mjd = self.obs.observation_span(n_samples=10**5, return_mjd=True)
        abs_phase_interp = polycos_model.eval_abs_phase(interp_topo_mjd).value

        self.polycos = interp1d(interp_topo_mjd.astype(np.float64), abs_phase_interp.astype(np.float64))
    
    def get_pulse(self, phase_abs, chan):
        if self.micro_structure:
            pulse_generator = MicroStructure(phase_abs, chan, self.micro_structure, self.period, self.observed_profile, self.seed)
            return pulse_generator.pulse_profile()
        else:
            return self.observed_profile(phase_abs % 1, chan)
    
    def coord2proper_time(self, bary_times):
        return bary_times+self.spin_ref - self.binary.orbital_delay(bary_times+self.orbit_ref)
    
    def get_phase(self, bary_times):
        T_proper = self.coord2proper_time(bary_times)
        phase_abs = self.phase_func(T_proper)
        return phase_abs 
    
    def generate_signal_polcos_bary(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        DM_delays = self.prop_effect.DM_delays[None, :]

        topo_times = self.obs.sec2mjd(timeseries)
        phase_time = topo_times[:, None] + DM_delays*u.s.to(u.day)

        bary_times = self.obs.topo2bary(timeseries, return_mjd=False, interp=True)

        phase = self.polycos(phase_time) + self.get_phase(bary_times[:, None] + DM_delays)
        gain_map = self.emission.gain(bary_times)
        return self.get_pulse(phase, np.arange(self.obs.n_chan)[None, :]) * gain_map
    
    def generate_signal_polcos_topo(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        DM_delays = self.prop_effect.DM_delays[None, :]

        topo_times = self.obs.sec2mjd(timeseries)
        phase_time = topo_times[:, None] + DM_delays*u.s.to(u.day)

        phase = self.polycos(phase_time) + self.get_phase(timeseries[:, None] + DM_delays)
        gain_map = self.emission.gain(timeseries)
        return self.get_pulse(phase, np.arange(self.obs.n_chan)[None, :]) * gain_map
        
    def generate_signal_python_bary(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        bary_times = self.obs.topo2bary(timeseries, return_mjd=False, interp=True)

        phase_array = self.get_phase(bary_times[:, None] + self.prop_effect.DM_delays[None, :])
        gain_map = self.emission.gain(bary_times)
        return self.get_pulse(phase_array, np.arange(self.obs.n_chan)[None, :]) * gain_map
    
    def generate_signal_python_topo(self, n_samples, sample_start=0):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)

        phase_array = self.get_phase(timeseries[:, None] + self.prop_effect.DM_delays[None, :])
        gain_map = self.emission.gain(timeseries)
        return self.get_pulse(phase_array, np.arange(self.obs.n_chan)[None, :]) * gain_map

   