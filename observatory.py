import numpy as np 
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, get_body, solar_system_ephemeris


class Observation:
    telescope_id = {64: 'Meerkat', 6: 'GBT'}

    def __init__(self, filterbank):
        fb_header = filterbank.header

        self.name = Observation.telescope_id[fb_header['telescope_id']]
        self.observatory = EarthLocation.of_site(self.name)
        self.source = SkyCoord(ra=self.convert_coord(fb_header['src_raj']), 
                               dec=self.convert_coord(fb_header['src_dej']), 
                               unit=(u.hourangle, u.deg), frame='gcrs',
                               obsgeoloc=np.array([*self.observatory.value])*u.m)
        
        self.ref_time = Time(fb_header['tstart'], format='mjd')
        self.dt = fb_header['tsamp']
        self.n_chan = fb_header['nchans']
        self.nbits = fb_header['nbits']
        self.n_samples = filterbank.n_samples
        self.fb_mean = filterbank.fb_mean
        self.fb_std = filterbank.fb_std
        
        self.freq_arr = np.linspace(fb_header['fch1'], 
                                    fb_header['fch1'] + fb_header['nchans']*fb_header['foff'], 
                                    fb_header['nchans'], endpoint=False)
        self.f0 = fb_header['fch1'] + fb_header['foff'] * (fb_header['nchans']-1)/2
        self.low_f = min(self.freq_arr)
        self.high_f = max(self.freq_arr)

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
    
    @staticmethod
    def set_ephemeris(ephemeris='de432s'):
        solar_system_ephemeris.set(ephemeris) 

    def get_topo_delay(self, dt, mjd): #todo
        times = np.linspace(mjd, mjd+dt*u.day, 1000)
        sun_refs =  get_body('sun', times, self.observatory)
        seps = self.source.separation(sun_refs)