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
        self.DM_const = (const.e.si**2/(8*np.pi**2*const.m_e*const.c) /(const.eps0) * u.pc.to(u.m)*u.m).value*1e-6   # Mhz^2 pc^-1 cm^3 s
        self.DM_delays = -self.DM * self.DM_const / self.obs.freq_arr**2  

    def scattering_relation_TPA(self, freq):
        DM_term = 3.6e-6 * self.DM**2.2 * (1 + 0.00194*self.DM**2)
        f_term = (freq/327)**-self.pulsar_pars.get('scattering_index', 4)
        Tau_s = DM_term * f_term

        return Tau_s/self.period

    def ISM_scattering(self, intrinsic_pulse):
        ref_scattering_time = self.pulsar_pars['scattering_time']

        if (ref_scattering_time == '') or (ref_scattering_time == 0):
            return intrinsic_pulse
        else:
            if ref_scattering_time == 'TPA':
                scattering_relation = self.scattering_relation_TPA
            else:
                scattering_index = self.pulsar_pars['scattering_index']
                ref_scattering_time = float(ref_scattering_time)
            
                def scattering_relation(freq):
                    scattering_phase = ref_scattering_time * u.ms.to(u.s) / self.period
                    return scattering_phase * (freq/self.obs.f0) ** -scattering_index 
                        
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
            if self.pulsar_pars['profile'] == 'default': # untested
                duty_cycle = self.pulsar_pars['duty_cycle']
                W_int = self.period * duty_cycle
                snr_int = np.sqrt((self.period-W_int) / W_int)

                def intrinsic_pulse(phase, duty_cycle=duty_cycle, snr_scale=1, chan_num=0): 
                    pulse_sigma = (duty_cycle)/(2*np.sqrt(2*np.log(2)))
                    return np.exp(-(phase-0.5)**2/(2*(pulse_sigma)**2)) * self.spectra(self.obs.freq_arr[chan_num]) * snr_scale
                
                self.smeared_values = []
                for chan in range(self.obs.n_chan):
                    channel_freq = self.obs.freq_arr[chan]
                    chan_top = channel_freq + self.obs.df/2
                    chan_bottom = channel_freq - self.obs.df/2
                    
                    W_eff = np.sqrt(W_int**2 + (self.DM_const * self.DM * (chan_top**-2 - chan_bottom**-2 ))**2)
                    snr_scale = np.sqrt((self.period-W_eff) / W_eff) / snr_int

                    self.smeared_values.append([W_eff, snr_scale])

                def smeared_pulse(phase, chan_num):
                    W_eff, snr_scale = self.smeared_values[chan_num]
                    return intrinsic_pulse(phase, W_eff/self.period, snr_scale, chan_num)

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
                    #     sys.exit(f"Unable to DM smear pulse profile for pulsar {self.pulsar_pars['ID']}.")

                def smeared_pulse(phase, chan_num):
                    return self.smeared_profiles[chan_num](phase)
                
            return smeared_pulse
        

