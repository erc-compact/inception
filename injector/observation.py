import sys
import numpy as np 
import astropy.units as u
from astropy.time import Time
import astropy.constants as const
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris, solar_system

from .io_tools import str2func


class Observation:
    telescope_id = {64: ['Meerkat', 'mk'], 6: ['GBT', 'gb']}

    def __init__(self, filterbank, ephem, pulsar_pars, generate=True):
        solar_system_ephemeris.set(ephem) 
        self.ephem = ephem

        fb_header = filterbank.header
        self.telescope_ID, self.tempo_id = Observation.telescope_id[fb_header['telescope_id']]
        self.observatory = EarthLocation.of_site(self.telescope_ID)

        self.obs_pointing = SkyCoord(ra=self.convert_coord(fb_header['src_raj']), 
                                     dec=self.convert_coord(fb_header['src_dej']), 
                                     unit=(u.hourangle, u.deg), frame='icrs')
        self.source = self.get_coords(pulsar_pars)
        self.beam_fwhm = str2func(pulsar_pars.get('beam_fwhm', 0), 'beam_fwhm', pulsar_pars['ID'], float) # beam reader, one beam/ multi beams- meta map
        
        self.obs_start = fb_header['tstart']
        self.dt = fb_header['tsamp']
        self.n_chan = fb_header['nchans']
        self.nbits = fb_header['nbits']
        self.n_samples = filterbank.n_samples
        self.obs_len = self.n_samples * self.dt
        self.fb_mean = filterbank.fb_mean
        self.fb_std = filterbank.fb_std
        
        self.freq_arr = np.linspace(fb_header['fch1'], 
                                    fb_header['fch1'] + fb_header['nchans']*fb_header['foff'], 
                                    fb_header['nchans'], endpoint=False)
        self.f0 = fb_header['fch1'] + fb_header['foff'] * (fb_header['nchans']-1)/2
        self.low_f = min(self.freq_arr)
        self.high_f = max(self.freq_arr)
        self.df = self.freq_arr[1]-self.freq_arr[0]

        self.DM = str2func(pulsar_pars.get('DM', None), 'DM', pulsar_pars['ID'], float)
        self.DM_const = (const.e.si**2/(8*np.pi**2*const.m_e*const.c) /(const.eps0) * u.pc.to(u.m)*u.m).value*1e-6   # Mhz^2 pc^-1 cm^3 s
        self.DM_delays = -self.DM * self.DM_const / self.freq_arr**2  

        self.obs_start_bary = self.topo2bary([self.obs_start])[0]
        if generate:
            self.barycentre_delays_interp = self.generate_interp()
        
        
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
        pulsar_ra, pulsar_dec = pulsar_pars.get('RAJ', None),  pulsar_pars.get('DECJ', None)
        pulsar_sep, pulsar_PA = pulsar_pars.get('separation', None),  pulsar_pars.get('position_angle', None)
        if pulsar_ra and pulsar_dec:
            try: 
                source = SkyCoord(ra=pulsar_ra, dec=pulsar_dec, unit=(u.hourangle, u.deg), frame='icrs')
            except ValueError:
                str2func(None, 'RA or DEC', pulsar_pars['ID'], float)
            else:
                return source
        elif pulsar_sep and pulsar_PA:
            sep_float = str2func(pulsar_sep, 'separation', pulsar_pars['ID'], float)
            PA_float = str2func(pulsar_PA, 'position_angle', pulsar_pars['ID'], float)
            return self.obs_pointing.directional_offset_by(PA_float*u.deg, sep_float*u.arcmin)
        else:
            return self.obs_pointing

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
    
    def observation_span(self, n_samples, mjd=True):
        lower_bound, upper_bound = np.min(self.DM_delays), self.obs_len + self.dt
        obs_sec = np.linspace(lower_bound, upper_bound, n_samples)
        if mjd:
            return obs_sec * u.s.to(u.day) + self.obs_start
        else:
            return obs_sec
    
    def barycentre_delays(self, topo_time):
        time_scale = Time(topo_time, format='mjd', scale='utc')
        L_hat  = self.source.cartesian.xyz.value.astype(np.float64)

        ep = solar_system.get_body_barycentric('earth', time_scale)
        op, _ = self.observatory.get_gcrs_posvel(time_scale)
        pos = ep.xyz.value.astype(np.float64) + op.xyz.value.astype(np.float64)*u.m.to(u.km)

        re_dot_L = np.sum(pos.T * L_hat, axis=1)
        return re_dot_L * u.km.to(u.lightsecond)
    
    def topo2bary(self, topo_times, mjd=True, interp=False):
        topo_times_tdb = Time(topo_times, format='mjd', scale='utc', precision=9).tdb.value

        if interp:
            bary_delays = self.barycentre_delays_interp(topo_times)
        else:
            bary_delays = self.barycentre_delays(topo_times)
        bary_times = Time(topo_times_tdb + bary_delays  * u.s.to(u.day), format='mjd', scale='tdb')

        if mjd:
            return bary_times.value
        else:
            return (bary_times.value - self.obs_start_bary)*u.day.to(u.s)
        
    def generate_interp(self):
        time_samples = self.observation_span(n_samples=10**4)
       
        delta_time = self.barycentre_delays(time_samples)
        return interp1d(time_samples, delta_time, kind='cubic')
    



    