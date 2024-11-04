import sys
import numpy as np 
import astropy.units as u
from astropy.time import Time
import astropy.constants as const
from scipy.integrate import quad_vec
from scipy.ndimage import convolve1d
from scipy.interpolate import interp1d


class PropagationEffects:
    def __init__(self, obs, pulsar_pars, profile_length, period, spectra):
        self.ID = pulsar_pars['ID']
        self.obs = obs
        self.pulsar_pars = pulsar_pars
        self.get_DM_delays()

        self.period = period
        self.spectra = spectra
        self.phase = np.linspace(0, 1, profile_length)
          
    def get_DM_delays(self):
        self.DM = self.pulsar_pars['DM']
        self.DM_const = (const.e.si**2/(8*np.pi**2*const.m_e*const.c) /(const.eps0) * u.pc.to(u.m)*u.m).value*1e-6   # Mhz^2 pc^-1 cm^3 s
        self.DM_delays = -self.DM * self.DM_const / self.obs.freq_arr**2  

    def scattering_relation_1(self, freq):
        DM_term = 0.154*np.log10(self.DM) + 1.07*np.log10(self.DM)**2
        f_term = -3.86*np.log10(freq * (u.MHz.to(u.GHz)))
        return 10**(-6.46 + DM_term + f_term) * u.ms.to(u.s) / self.period

    def ISM_scattering(self, intrinsic_pulse):
        ref_scattering_time = self.pulsar_pars['scattering_time']

        if (ref_scattering_time == 0):
            return intrinsic_pulse
        else:
            if ref_scattering_time == 12321:
                scattering_relation = self.scattering_relation_1
            else:
                scattering_index = self.pulsar_pars['scattering_index']
            
                def scattering_relation(freq):
                    scattering_phase = ref_scattering_time * u.ms.to(u.s) / self.period
                    return scattering_phase * (freq/self.obs.f0) ** scattering_index 
                        
            def scattering_kernal(nchan):
                scattering_time = scattering_relation(self.obs.freq_arr[nchan])
                kernal = np.exp(-self.phase/scattering_time)
                return kernal / np.sum(kernal)
            
            conv_profiles = [convolve1d(intrinsic_pulse(self.phase, nchan), scattering_kernal(nchan), mode='wrap') for nchan in range(self.obs.n_chan)]
            phase_shift = [np.roll(profile, int(round(len(self.phase)/2))) for profile in conv_profiles]
            self.scattered_profiles = [interp1d(self.phase, profile) for profile in phase_shift]

            def scattered_pulse(phase, chan_num):
                return self.scattered_profiles[chan_num](phase)
            
            return scattered_pulse
    

    def intra_channel_DM_smearing(self, intrinsic_pulse):
        
        def channel_profile(phase, freq, channel_freq, chan_num):
            dt = self.DM_const * self.DM * (channel_freq**-2 - freq**-2)
            phase_freq = phase + dt/self.period
            PSD_corr = self.spectra(freq)/self.spectra(channel_freq)
            return intrinsic_pulse(phase_freq % 1, chan_num) * PSD_corr
        
        if not self.pulsar_pars['DM_smear']:
            return intrinsic_pulse
        else:
            self.smeared_profiles = []
            for chan in range(self.obs.n_chan):
                channel_freq = self.obs.freq_arr[chan]
                channel_top = abs(self.obs.df)/2 + channel_freq
                channel_bottom = -abs(self.obs.df)/2 + channel_freq

                integrate_freq = quad_vec(lambda freq: channel_profile(self.phase, freq, channel_freq, chan),
                                        a=channel_bottom, b=channel_top,
                                        epsabs=1e-3, epsrel=1e-3)
                
                self.smeared_profiles.append(interp1d(self.phase, integrate_freq[0]/np.abs(self.obs.df)))
                # if integrate_freq[1] < 1:
                #     sys.exit(f'Unable to DM smear pulse profile for pulsar {self.pulsar_pars['ID']}.')

            def smeared_pulse(phase, chan_num):
                return self.smeared_profiles[chan_num](phase)
                
            return smeared_pulse
        

