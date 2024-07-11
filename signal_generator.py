import numpy as np 
import astropy.units as u
import astropy.constants as const


class PulsarSignal:
    def __init__(self, pulse_period, duty_cycle, spectral_index, SNR, DM, obs, binary):

        self.obs = obs
        self.binary = binary

        self.p0 = pulse_period / 1000
        self.duty_cycle = duty_cycle
        self.phase_offset = 0.5
        self.pulse_sigma = (duty_cycle*self.p0)/(2*np.sqrt(2*np.log(2)))
        self.spectral_index = spectral_index
        self.power_constant = self.SNR2power(SNR)/(self.obs.low_f)**spectral_index

        self.DM = DM
        self.DM_const = (const.e.si**2/(8*np.pi**2*const.m_e*const.c) /(const.eps0) * u.pc.to(u.m)*u.m).value*1e-6 
        #^ Mhz^2 pc^-1 cm^3 s
        self.DM_delays = self.DM * self.DM_const * (1/self.obs.high_f**2 - 1/self.obs.freq_arr**2)

    def SNR2power(self, SNR):
        sigma_pt = self.obs.fb_std
        n_sample = self.obs.n_samples
        Weq_t = self.duty_cycle/(2*np.sqrt(2*np.log(2)))*np.sqrt(2*np.pi)
        Amp = SNR * sigma_pt / (np.sqrt(Weq_t) * np.sqrt(self.obs.n_chan) * np.sqrt(n_sample))
        return Amp

    def spectral_power(self, freq):
        return self.power_constant * (freq)**self.spectral_index

    def pulse(self, timeseries, freq):
        phase = timeseries%self.p0 - self.phase_offset*self.p0
        return  self.spectral_power(freq) * np.exp(-(phase)**2/(2*self.pulse_sigma**2))
    
    def period_doppler_shift(self, vr, rel=False):
        beta =  vr/const.c.value
        if rel:
            return self.p0 * np.sqrt((1 + beta) / (1 - beta))
        else:
            return self.p0 * (1 + beta)
    
    def get_pulse_info(self, timeseries):
        total_time = len(timeseries)*self.obs.dt + max(self.DM_delays)
        n_pulses = np.round(total_time/self.p0)
        already_pulsed = np.floor((timeseries[0]-self.obs.dt)/self.p0)
        return int(n_pulses), max(int(already_pulsed), 0)

    def generate_signal(self, n_samples, sample_start=0, bary=False):
        timeseries = np.linspace(self.obs.dt*sample_start, self.obs.dt*(n_samples+sample_start-1), n_samples)
        time_array = np.tile(timeseries, (len(self.obs.freq_arr),1)).T

        freq_array = np.tile(self.obs.freq_arr, (len(timeseries),1))
        DM_array = np.tile(self.DM_delays, (len(timeseries),1))
        pulse_offset = - DM_array

        if bary:
            bary_delay = self.obs.topo_delay(timeseries*u.s.to(u.day) + self.obs.ref_time)
            bary_delay_arr = np.tile(bary_delay, (len(self.obs.freq_arr),1)).T
            pulse_offset += bary_delay_arr

        if self.binary.period != 0:
            orbit_delay = self.binary.get_roemer_delay(timeseries) 
            orbit_delay_chan = np.tile(orbit_delay, (len(self.obs.freq_arr),1)).T
            pulse_offset += orbit_delay_chan

        filterbank = self.pulse(time_array-pulse_offset, freq_array)
        return filterbank
