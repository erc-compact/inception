import sys
import numpy as np 
import pandas as pd
from pathlib import Path
from scipy.interpolate import interp1d, PchipInterpolator


class PulsarEmission:
    def __init__(self, obs, pulsar_pars):
        self.obs = obs
        self.pulsar_pars = pulsar_pars
        self.ID = pulsar_pars['ID']

        self.spectrum = pulsar_pars['spectrum']
        self.gain_path = pulsar_pars['gain_map']
        self.profile = pulsar_pars['profile']

        self.build_emission_model()

    def build_emission_model(self):

        self.spectra = self.get_spectra()
        self.gain = self.get_gain_map()
        self.intrinsic_profile_chan = self.get_intrinsic_profile()

    def get_spectra(self):            
        if type(self.spectrum) == str:
            suffix = Path(self.spectrum).suffix
            if suffix == '.npy':
                try:
                    spectra_arr = np.load(self.spectrum)
                except FileNotFoundError:
                    sys.exit(f'Unable to load {self.spectrum} numpy spectra for pulsar {self.ID}.')
            else:
                sys.exit(f'Pulsar {self.ID} has an invalid spectra file extension: {self.spectrum}. Must be a numpy .npy file.')

            freq_min = np.min(self.obs.freq_arr) - abs(self.obs.df)/2
            freq_max = np.max(self.obs.freq_arr) + abs(self.obs.df)/2
            freq_range = np.linspace(freq_min, freq_max, len(spectra_arr))

            norm = np.mean(np.abs(spectra_arr))
            return lambda freq: np.interp(freq, freq_range, spectra_arr) / norm
        else:
            return lambda freq: (freq/self.obs.f0)**float(self.spectrum)
        
    def get_gain_map(self):
        if self.gain_path:
            try:
                gain = np.load(self.gain_path)
            except FileNotFoundError:
                sys.exit(f'Unable to load {self.gain_path} numpy gain map for pulsar {self.ID}.')
        else:
            return lambda t: np.ones((len(np.atleast_1d(t)), self.obs.n_chan))

        n_chan = self.obs.n_chan
        gain_axis = self.pulsar_pars.get('gain_axis', 'time')
        if gain.ndim == 1:
            if gain_axis == 'time':
                n_subints = len(gain)
                gain_times = (np.arange(n_subints) + 0.5) * self.obs.obs_len / n_subints
                gain = np.tile(gain[:, None], (1, n_chan))

            elif gain_axis == 'freq':
                if len(gain) != n_chan:
                    raise ValueError(f'Frequency gain must have {n_chan} values, got {len(gain)}.')

                gain_times = np.array([0., self.obs.obs_len])
                gain = np.tile(gain, (2, 1))

            else:
                raise ValueError("gain_axis must be 'time' or 'freq'.")

        elif gain.ndim == 2:
            if gain.shape[0] != n_chan:
                raise ValueError(f'2D gain must have shape ({n_chan}, n_subints), got {gain.shape}.')

            n_subints = gain.shape[1]

            gain_times = (np.arange(n_subints) + 0.5) * self.obs.obs_len / n_subints
            gain = gain.T
        else:
            raise ValueError('Gain map must be a 1D or 2D numpy array.')

        return interp1d(gain_times, gain, axis=0,
                        kind='linear',bounds_error=False, fill_value=(gain[0], gain[-1]))
        
    def get_intrinsic_profile(self):

        if self.profile == "default":
            phase_arr, profile_arr = self.build_default_profile()

        elif isinstance(self.profile, str):
            phase_arr, profile_arr = self.load_profile()

        elif isinstance(self.profile, dict):
            phase_arr, profile_arr = self.build_parametric_profile()

        else:
            sys.exit(f"Invalid pulse profile for pulsar {self.ID}.")

        return self.make_intrinsic_profile(phase_arr, profile_arr)
    
    @staticmethod
    def single_pulse(phase, centre, sigma):
        return np.exp(-(phase - centre)**2 / (2 * sigma**2))
    
    @staticmethod
    def DC2sigma(DC):
        return DC / (2 * np.sqrt(2 * np.log(2)))
    
    def build_default_profile(self):
        duty_cycle = self.pulsar_pars["duty_cycle"]
        pulse_sigma = self.DC2sigma(duty_cycle)

        profile_length = 1000
        phase = np.linspace(-1, 2, profile_length * 3)

        profile = (
            self.single_pulse(phase, 0.5, pulse_sigma)
            + self.single_pulse(phase, -0.5, pulse_sigma)
            + self.single_pulse(phase, 1.5, pulse_sigma)
        )

        mask = (phase > -0.5) & (phase < 1.5)
        return phase[mask], profile[mask]
    
    def make_intrinsic_profile(self, phase_arr, profile_arr):
        spectral_weights = self.spectra(self.obs.freq_arr)

        if profile_arr.ndim == 1:

            interp = interp1d(phase_arr, profile_arr/np.max(profile_arr))
            self.profile_length = profile_arr.size

            return lambda phase, chan=0: interp(phase) * spectral_weights[chan]

        elif (profile_arr.ndim == 2) and (profile_arr.shape[0] == self.obs.n_chan):
            phase_arr = np.linspace(0, 1, len(profile_arr.T))
            interps = [interp1d(phase_arr, row / np.max(profile_arr)) for row in profile_arr]
            self.profile_length = profile_arr.shape[1]

            return lambda phase, chan: interps[chan](phase) * spectral_weights[chan]
        
        # elif (profile_arr.ndim == 3) and (profile_arr.shape[0] == self.obs.n_chan):
        #     self.profile_length = profile_arr.shape[2]
        #     n_subints = profile_arr.shape[1]
        #     dt_sub = self.obs.obs_len / n_subints
        #     time_arr = (np.arange(n_subints) + 0.5) * dt_sub

        #     interps = [interp1d(time_arr, profile_arr[chan] / np.max(profile_arr), axis=0, bounds_error=False, fill_value="extrapolate")
        #                 for chan in range(self.obs.n_chan)]
            
        #     return lambda time, chan: interps[chan](time) * spectral_weights[chan]

        else:
            sys.exit(f'Unable to interpret {self.profile} numpy pulse profile for pulsar {self.ID}.')
        
    def profile_component(self, phase_range, pulse_i):
        DC = self.profile["duty_cycle"][pulse_i]
        centre = self.profile["phase"][pulse_i]
        amp = self.profile["amp"][pulse_i]

        pulse_sigma = self.DC2sigma(DC)
        pulse = self.single_pulse(phase_range, centre, pulse_sigma)

        pulse /= pulse.max()
        return pulse * amp

    def build_parametric_profile(self):
        profile_length = 1000
        phase_range = np.linspace(0, 1, profile_length)
        multi_profile = np.zeros(profile_length)

        for i in range(len(self.profile["phase"])):
            multi_profile += self.profile_component(phase_range, i)

        return phase_range, multi_profile

    def load_profile(self):
        suffix = Path(self.profile).suffix
        if suffix == ".npy":
            try:
                profile_arr = np.load(self.profile)
                phase_range = np.linspace(0, 1, len(profile_arr))
                return phase_range, profile_arr
            
            except FileNotFoundError:
                sys.exit(f"Unable to load {self.profile} numpy pulse profile for pulsar {self.ID}.")

        elif suffix == ".txt":
            try:
                epn_profile = pd.read_csv(self.profile, delimiter=' ', 
                                        names=['col0', 'col1', 'col2', 'intensity'])
                profile_arr = epn_profile['intensity'].values
                phase_range = np.linspace(0, 1, len(profile_arr))
                return phase_range, profile_arr

            except (FileNotFoundError, KeyError, ValueError):
                sys.exit(f"Unable to load {self.profile} EPN pulse profile for pulsar {self.ID}.")

        sys.exit(f"Pulsar {self.ID} has an invalid profile file extension: {self.profile}. Must be a numpy .npy or EPN .txt file.")
    