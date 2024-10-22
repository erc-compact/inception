import sys
import numpy as np 
import astropy.units as u
from astropy.time import Time
import astropy.constants as const
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris, solar_system

from .propagation_effects import PropagationEffects

class Observation:
    telescope_id = {64: ['Meerkat', 'mk'],  1: ['Arecibo', 'ao'], 4: ['Parkes', 'pk'], 5: ['Jodrell', 'jb'], 
                    6: ['GBT', 'gb'], 7: ['GMRT', 'gm'], 8: ['Effelsberg', 'ef']}

    def __init__(self, filterbank, ephem, pulsar_pars, generate=[]):
        solar_system_ephemeris.set(ephem) 
        self.ephem = ephem
        fb_header = filterbank.header
        
        self.get_pointing_data(fb_header, pulsar_pars)
        self.get_beam_data(pulsar_pars)
        self.get_obs_data(fb_header, filterbank)
        self.get_freq_data(fb_header)

        self.obs_start_bary = self.topo2bary([self.obs_start])[0]
        if generate:
            self.prop_effect = PropagationEffects(self, pulsar_pars, 1, '', '')
            self.barycentre_delays_interp = self.generate_interp(generate)

    def get_pointing_data(self, fb_header, pulsar_pars):
        self.telescope_ID, self.tempo_id = Observation.telescope_id[fb_header['telescope_id']]
        self.observatory = EarthLocation.of_site(self.telescope_ID)

        self.obs_pointing = SkyCoord(ra=self.convert_coord(fb_header['src_raj']), 
                                     dec=self.convert_coord(fb_header['src_dej']), 
                                     unit=(u.hourangle, u.deg), frame='icrs')
        self.source = self.get_coords(pulsar_pars)
    
    def get_beam_data(self, pulsar_pars):
        self.beam_fwhm = pulsar_pars['beam_fwhm'] # beam reader, one beam/ multi beams- meta map

    def get_obs_data(self, fb_header, filterbank):
        self.obs_start = fb_header['tstart']
        self.dt = fb_header['tsamp']
        self.n_chan = fb_header['nchans']
        self.nbits = fb_header['nbits']
        self.n_samples = filterbank.n_samples
        self.obs_len = self.n_samples * self.dt
        self.fb_mean = filterbank.fb_mean
        self.fb_std = filterbank.fb_std

    def get_freq_data(self, fb_header):
        self.freq_arr = np.linspace(fb_header['fch1'], 
                                    fb_header['fch1'] + fb_header['nchans']*fb_header['foff'], 
                                    fb_header['nchans'], endpoint=False)
        self.f0 = fb_header['fch1'] + fb_header['foff'] * (fb_header['nchans']-1)/2
        self.low_f = min(self.freq_arr)
        self.high_f = max(self.freq_arr)
        self.df = self.freq_arr[1]-self.freq_arr[0]
        
    @staticmethod
    def convert_coord(coord_str):
        coord_str = str(coord_str)
        is_negative = coord_str.startswith('-')
        if is_negative:
            coord_str = coord_str[1:]
        
        hours = coord_str[:2]
        minutes = coord_str[2:4]
        seconds = coord_str[4:]
        formatted_coord = f'{hours}:{minutes}:{seconds}'
        if is_negative:
            formatted_coord = '-' + formatted_coord

        return formatted_coord   
    
    def get_coords(self, pulsar_pars):   
        pulsar_ra, pulsar_dec = pulsar_pars['RAJ'],  pulsar_pars['DECJ']
        if pulsar_ra and pulsar_dec:
            try: 
                source = SkyCoord(ra=pulsar_ra, dec=pulsar_dec, unit=(u.hourangle, u.deg), frame='icrs')
            except ValueError:
                sys.exit(f"Invalid RA/DEC for pulsar {pulsar_pars['ID']}")
            else:
                return source
        else:
            return self.obs_pointing.directional_offset_by(pulsar_pars['position_angle']*u.deg, pulsar_pars['separation']*u.arcmin)

    def get_beam_snr(self):
        if self.beam_fwhm != 0:
            beam_sigma = self.beam_fwhm/(2*np.sqrt(2*np.log(2))) # arcmin
            beam_offset = self.source.separation(self.obs_pointing).arcmin
            beam_SNR = np.exp(-(beam_offset)**2/(2*(beam_sigma)**2))
            return beam_SNR
        else:
            return 1
    
    def sec2mjd(self, time_sec):
        return time_sec * u.s.to(u.day) + self.obs_start
    
    def observation_span(self, obs_range, n_samples, mjd=True):
        if type(obs_range) == list:
            pad_time = 1
            lower_bound, upper_bound = np.min(self.prop_effect.DM_delays) + obs_range[0] * self.dt - pad_time, obs_range[1] * self.dt + pad_time
        else:
            lower_bound, upper_bound = np.min(self.prop_effect.DM_delays), self.obs_len + self.dt
        obs_sec = np.linspace(lower_bound, upper_bound, n_samples)
        if mjd:
            return obs_sec * u.s.to(u.day) + self.obs_start
        else:
            return obs_sec
    
    def topo2bary_calc(self, topo_time, mjd=True):
        time_scale = Time(topo_time, format='mjd', scale='utc')
        L_hat  = self.source.cartesian.xyz.value.astype(np.float64)

        ep = solar_system.get_body_barycentric('earth', time_scale)
        op, _ = self.observatory.get_gcrs_posvel(time_scale)
        pos = ep.xyz.value.astype(np.float64) + op.xyz.value.astype(np.float64)*u.m.to(u.km)

        re_dot_L = np.sum(pos.T * L_hat, axis=1)
        bary_delays_sec = re_dot_L * u.km.to(u.lightsecond) 
        bary_times = time_scale.tdb.value + bary_delays_sec * u.s.to(u.day)
        if mjd:
            return bary_times
        else:
            return (bary_times - self.obs_start_bary)*u.day.to(u.s)
    
    def topo2bary(self, topo_times, mjd=True, interp=False):
        if interp:
            bary_delays = self.barycentre_delays_interp(topo_times)
        else:
            bary_delays = self.topo2bary_calc(topo_times, mjd=mjd)
        return bary_delays
        
    def generate_interp(self, obs_range):
        time_samples = self.observation_span(obs_range, n_samples=10**4)
       
        delta_time = self.topo2bary_calc(time_samples, mjd=False)
        return interp1d(time_samples, delta_time, kind='cubic')
    



    