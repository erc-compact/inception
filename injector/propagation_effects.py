import sys
import numpy as np 
import astropy.units as u
from astropy.time import Time
import astropy.constants as const
from scipy.integrate import quad_vec
from scipy.ndimage import convolve1d
from scipy.interpolate import interp1d


class PropagationEffects:
    def __init__(self, obs, pulsar_pars, profile_length=1, period=0, spectra=0):
        self.ID = pulsar_pars['ID']
        self.obs = obs
        self.pulsar_pars = pulsar_pars
        self.get_DM_delays()

        self.period = period
        self.spectra = spectra
        self.phase = np.linspace(0, 1, profile_length)
          
    def get_DM_delays(self):
        self.DM = self.pulsar_pars['DM']
        self.cDM = self.pulsar_pars.get('cDM', 0)
        self.DM_const = (const.e.si**2/(8*np.pi**2*const.m_e*const.c) /(const.eps0) * u.pc.to(u.m)*u.m).value*1e-6   # MHz^2 pc^-1 cm^3 s

        ref_freq = self.pulsar_pars['DM_ref']
        if ref_freq == 'inf':
            self.DM_delays = -self.DM * self.DM_const / self.obs.freq_arr**2  
        else:
            self.DM_delays = -self.DM * self.DM_const * (1/self.obs.freq_arr**2  - 1/self.obs.high_f**2)

    def scattering_relation_TPA(self, freq):
        DM_term = 3.6e-6 * self.DM**2.2 * (1 + 0.00194*self.DM**2)
        f_term = (freq/327)**-self.pulsar_pars.get('scattering_index', 4)
        Tau_s = DM_term * f_term / 1000

        return Tau_s
    
    def scattering_relation_CORDES(self, freq):
        f_term = (freq/1000)**-self.pulsar_pars.get('scattering_index', 4)
        A_ms = 2.98e-7
        a = 1.4
        B = 3.55e-5
        b = 3.1

        sigma_T = 0.76 * self.pulsar_pars.get('scattering_sigma', 0)

        A = A_ms/1000
        Tau_s = A * self.DM ** a * (1 + B * self.DM ** b) * f_term

        return 10 ** (np.log10(Tau_s) + sigma_T)

    def ISM_scattering(self, intrinsic_pulse):
        ref_scattering_time = self.pulsar_pars['scattering_time']

        if (ref_scattering_time == '') or (ref_scattering_time == '0') or (ref_scattering_time == 0):
            return intrinsic_pulse
        else:
            if ref_scattering_time == 'TPA':
                scattering_relation = self.scattering_relation_TPA
            elif ref_scattering_time == 'CORDES':
                scattering_relation = self.scattering_relation_CORDES
            else:
                scattering_index = self.pulsar_pars['scattering_index']
                ref_scattering_time = float(ref_scattering_time)
            
                def scattering_relation(freq):
                    scattering_phase = ref_scattering_time * u.ms.to(u.s) 
                    return scattering_phase * (freq/self.obs.f0) ** -scattering_index 
                        
            def scattering_kernel(nchan):
                scattering_time = scattering_relation(self.obs.freq_arr[nchan]) / self.period
                kernal = np.exp(-self.phase/scattering_time)
                return kernal / np.sum(kernal)
            
            conv_profiles = [convolve1d(intrinsic_pulse(self.phase, nchan), scattering_kernel(nchan), mode='wrap') for nchan in range(self.obs.n_chan)]
            phase_shift = [np.roll(profile, int(round(len(self.phase)/2))) for profile in conv_profiles]
            self.scattered_profiles = [interp1d(self.phase, profile) for profile in phase_shift]

            def scattered_pulse(phase, chan_num):
                return self.scattered_profiles[chan_num](phase)
            
            return scattered_pulse
        
    def intra_channel_DM_smearing(self, intrinsic_pulse):

        if self.pulsar_pars["DM_smear"] == 'off':
            return intrinsic_pulse

        elif (self.pulsar_pars["profile"] == "default") and (self.pulsar_pars["DM_smear"] == 'approx'):
            return self.gaussian_DM_smearing()

        elif (self.pulsar_pars["DM_smear"] == 'exact'):
            return self.numerical_DM_smearing(intrinsic_pulse)
        
        else:
            sys.exit(f"Invalid DM_smear flag for pulsar {self.ID}. 'approx' can only be used for profile == 'defualt'")

    def gaussian_DM_smearing(self):
        duty_cycle = self.pulsar_pars["duty_cycle"]
        W_int = self.period * duty_cycle
        
        def single_pulse(phase, DC):
            sigma = DC / (2 * np.sqrt(2 * np.log(2)))
            return np.exp(-(phase - 0.5) ** 2 / (2 * sigma**2))

        self.smeared_profiles = []
        for chan, freq in enumerate(self.obs.freq_arr):

            chan_top = freq + self.obs.df / 2
            chan_bottom = freq - self.obs.df / 2

            smear = self.DM_const * (self.DM - self.cDM) * (chan_top**-2 - chan_bottom**-2)
            W_eff = np.sqrt(W_int**2 + smear**2)

            eff_DC = W_eff / self.period
            phase = np.linspace(-1, 2, len(self.phase) * 3)

            profile = (
                single_pulse(phase, eff_DC)
                + single_pulse(phase - 1, eff_DC)
                + single_pulse(phase + 1, eff_DC)
            )

            profile /= profile.max()
            profile *= W_int / W_eff

            mask = (phase > -0.5) & (phase < 1.5)
            interp = interp1d(phase[mask], profile[mask] * self.spectra(freq))
            self.smeared_profiles.append(interp)

        def smeared_pulse(phase, chan_num):
            return self.smeared_profiles[chan_num](phase)

        return smeared_pulse
    
    def numerical_DM_smearing(self, intrinsic_pulse):

        def dm_smearing_kernel(nchan):
            freq = self.obs.freq_arr[nchan]
            chan_top = freq + abs(self.obs.df) / 2
            chan_bottom = freq - abs(self.obs.df) / 2

            smear = abs(self.DM_const * (self.DM - self.cDM) * (chan_top**-2 - chan_bottom**-2))

            smear_phase = smear / self.period
            width = max(1, int(round(smear_phase * len(self.phase))))
            width = min(width, len(self.phase))
            if width % 2 == 0:
                width += 1

            kernel = np.ones(width)
            kernel /= kernel.sum()

            return kernel

        conv_profiles = [convolve1d(intrinsic_pulse(self.phase, nchan), dm_smearing_kernel(nchan), mode="wrap") for nchan in range(self.obs.n_chan)]
        self.smeared_profiles = [ interp1d( self.phase, profile, bounds_error=False, fill_value="extrapolate") for profile in conv_profiles ]

        def smeared_pulse(phase, chan_num):
            return self.smeared_profiles[chan_num](phase)

        return smeared_pulse

    

    