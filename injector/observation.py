import sys
import numpy as np 
import astropy.units as u
from astropy.time import Time, TimeDelta
import astropy.constants as const
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris, solar_system, Distance, CartesianRepresentation

from .propagation_effects import PropagationEffects


class Observation:
    telescope_id = {64: ['Meerkat', 'mk'],  1: ['Arecibo', 'ao'], 4: ['Parkes', 'pk'], 5: ['Jodrell', 'jb'], 
                    6: ['GBT', 'gb'], 7: ['GMRT', 'gm'], 8: ['Effelsberg', 'ef']}

    def __init__(self, filterbank, ephem, pulsar_pars, generate=[], override_length=0):
        solar_system_ephemeris.set(ephem) 
        self.ephem = ephem
        fb_header = filterbank.header
        
        self.get_obs_data(fb_header, filterbank, override_length)
        self.get_pointing_data(fb_header, pulsar_pars)
        self.get_beam_data(pulsar_pars)
        self.get_freq_data(fb_header)
        
        self.obs_start_bary = self.topo2bary([self.obs_start], return_mjd=True)[0]

        self.obs_start_time_TIME = Time(self.obs_start, format="mjd", scale="utc")
        self.obs_start_bary_TIME = Time(self.obs_start_bary, format="mjd", scale="tdb")
        self.set_epochs(pulsar_pars)

        self.prop_effect = PropagationEffects(self, pulsar_pars)
        if generate:
            self.barycentre_delays_interp = self.generate_interp(generate)

    def get_pointing_data(self, fb_header, pulsar_pars):
        self.telescope_ID, self.tempo_id = Observation.telescope_id[fb_header['telescope_id']]
        self.observatory = EarthLocation.of_site(self.telescope_ID)

        self.obs_pointing = SkyCoord(ra=self.convert_coord(fb_header['src_raj']), 
                                     dec=self.convert_coord(fb_header['src_dej']), 
                                     unit=(u.hourangle, u.deg), frame='icrs')
        self.source = self.get_coords(pulsar_pars)
    
    def get_beam_data(self, pulsar_pars):
        self.beam_fwhm = pulsar_pars.get('beam_fwhm', 0)

    def get_obs_data(self, fb_header, filterbank, override_length):
        self.obs_start = fb_header['tstart']
        self.dt = fb_header['tsamp']
        self.n_chan = fb_header['nchans']
        self.nbits = fb_header['nbits']
        if override_length == 0:
            self.n_samples = filterbank.n_samples
        else:
            self.n_samples = override_length
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
        integer, decimal = str(coord_str).split('.')
        coord_str = f"{'0'*(6-len(integer))}{integer}.{decimal}"

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
    
    def set_epochs(self, pulsar_pars):

        if pulsar_pars['frame'] == "topo":
            obs_ref = self.obs_start
        elif pulsar_pars['frame'] == "bary":
            obs_ref = self.obs_start_bary

        pepoch = pulsar_pars['PEPOCH'] if pulsar_pars['PEPOCH'] else obs_ref
        if 0 <= np.abs(pepoch) <= 1:
            pepoch = obs_ref + pulsar_pars['PEPOCH'] * self.obs_len * u.s.to(u.day)

        T0 = pulsar_pars['T0'] if pulsar_pars['T0'] else obs_ref

        spin_ref = (obs_ref - pepoch) * u.day.to(u.s)
        orbit_ref = (obs_ref - T0) * u.day.to(u.s)

        self.T0 = T0
        self.pepoch = pepoch
        self.spin_ref = spin_ref
        self.orbit_ref = orbit_ref
        self.accepoch = pulsar_pars['ACCEPOCH'] * self.obs_len

    
    def get_coords(self, pulsar_pars):   
        self.posepoch = pulsar_pars['POSEPOCH'] if pulsar_pars['POSEPOCH'] else Time(self.obs_start, format="mjd", scale="utc").tdb.mjd
        if 0 <= np.abs(self.posepoch) <= 1:
            self.posepoch = self.obs_start + self.posepoch * self.obs_len * u.s.to(u.day)

        pulsar_ra, pulsar_dec = pulsar_pars.get('RAJ', None),  pulsar_pars.get('DECJ', None)
        pm_ra_cosdec, pm_dec = pulsar_pars.get('PMRA', 0),  pulsar_pars.get('PMDEC', 0)
        self.radial_velocity = pulsar_pars.get('RV', 0.0)

        self.parallax = pulsar_pars.get('parallax', None)
        self.dist = pulsar_pars.get('DIST', 10.0)

        if self.parallax is not None:
            self.dist = Distance(parallax=self.parallax * u.mas).kpc
        else:
            self.parallax = 1/self.dist
        
        if pulsar_ra and pulsar_dec:
            try: 
                source = SkyCoord(ra=pulsar_ra, dec=pulsar_dec, unit=(u.hourangle, u.deg), frame='icrs',
                                  pm_ra_cosdec=pm_ra_cosdec * u.mas/u.yr, pm_dec=pm_dec * u.mas/u.yr,
                                  distance=self.dist * u.kpc, radial_velocity=self.radial_velocity * u.km/u.s,
                                  obstime=Time(self.posepoch, format="mjd", scale="tdb"))

            except ValueError:
                sys.exit(f"Invalid RA/DEC for pulsar {pulsar_pars['ID']}")
            else:
                return source
        else:
            self.pos_ang = pulsar_pars.get('position_angle', 0) * u.deg
            self.sep = pulsar_pars.get('separation', 0) * u.arcmin

            source = self.obs_pointing.directional_offset_by(self.pos_ang, self.sep)
            source = SkyCoord(ra=source.ra, dec=source.dec, unit=(u.hourangle, u.deg), frame='icrs',
                                pm_ra_cosdec=pm_ra_cosdec * u.mas/u.yr, pm_dec=pm_dec * u.mas/u.yr,
                                distance=self.dist * u.kpc, radial_velocity=self.radial_velocity * u.km/u.s,
                                obstime=Time(self.posepoch, format="mjd", scale="tdb"))
        return source
        
    def get_source(self, topo_time):
        dt = (topo_time.tdb - self.source.obstime).to(u.s)

        r0 = self.source.cartesian
        v = self.source.velocity

        r = CartesianRepresentation(
            x=r0.x + v.d_x * dt,
            y=r0.y + v.d_y * dt,
            z=r0.z + v.d_z * dt,
        )
        return SkyCoord(r, frame=self.source.frame, obstime=topo_time.tdb)

    def get_beam_snr(self):
        if self.beam_fwhm != 0:
            beam_sigma = self.beam_fwhm/(2*np.sqrt(2*np.log(2)))
            beam_offset = self.source.separation(self.obs_pointing).arcmin
            beam_SNR = np.exp(-(beam_offset)**2/(2*(beam_sigma)**2))
            return beam_SNR
        else:
            return 1
    
    def sec2mjd(self, time_sec):
        return (self.obs_start_time_TIME + TimeDelta(time_sec, format="sec")).mjd
    
    def observation_span(self, obs_range, n_samples, return_mjd=True):
        if type(obs_range) == list:
            pad_time = 1
            lower_bound = np.min(self.prop_effect.DM_delays) + obs_range[0] * self.dt - pad_time
            upper_bound = obs_range[1] * self.dt + pad_time
        else:
            lower_bound = np.min(self.prop_effect.DM_delays)
            upper_bound = self.obs_len + self.dt

        obs_sec = np.linspace(lower_bound, upper_bound, n_samples)

        if return_mjd:
            return (self.obs_start_time_TIME + TimeDelta(obs_sec, format="sec")).mjd
        else:
            return obs_sec
    
    def topo2bary_calc(self, topo_time, return_mjd=True):
        if isinstance(topo_time, Time):
            time_scale = topo_time
        else:
            time_scale = Time(topo_time, format="mjd", scale="utc")

        source = self.get_source(time_scale)
        
        ep = solar_system.get_body_barycentric('earth', time_scale)
        op, _ = self.observatory.get_gcrs_posvel(time_scale)
        pos = ep.xyz.value.astype(np.float64) + op.xyz.value.astype(np.float64)*u.m.to(u.km)

        r_p = source.cartesian.xyz.to_value(u.km).T
        d = np.linalg.norm(r_p, axis=1)

        L_hat = r_p / d[:, None]

        dot = np.sum(pos.T * L_hat, axis=1)
        r2 = np.sum(pos.T**2, axis=1)

        delay = dot - (r2 - dot**2)/(2*d)
        bary_delays_sec = delay / const.c.to_value(u.km/u.s)
        bary_times = time_scale.tdb + TimeDelta(bary_delays_sec, format="sec")
        if return_mjd:
            return bary_times.mjd
        else:
            return (bary_times - self.obs_start_bary_TIME).sec
        
    def bary2topo_calc(self, bary_time, return_mjd=True, tol=1e-12, max_iter=5):
        if isinstance(bary_time, Time):
            bary_time = bary_time.tdb
        else:
            bary_time = Time(bary_time, format="mjd", scale="tdb")

        topo_time = bary_time.utc
        for _ in range(max_iter):
            bary_est = self.topo2bary_calc(topo_time, return_mjd=True)
            dt = bary_time.mjd - bary_est

            topo_time = topo_time + TimeDelta(dt, format="jd")
            if np.max(np.abs(dt)) < tol:
                break

        if return_mjd:
            return topo_time.mjd
        else:
            return (topo_time - self.obs_start_time_TIME).sec
        
    def earth_radial_velocity(self, topo_time):
        topo_time = np.atleast_1d(topo_time)

        time_scale = Time(topo_time, format='mjd', scale='utc')

        _, ep_vel = solar_system.get_body_barycentric_posvel('earth', time_scale)
        _, obs_vel = self.observatory.get_gcrs_posvel(time_scale)

        source = self.get_source(time_scale)
        r_p = source.cartesian.xyz.to_value(u.km).T
        d = np.linalg.norm(r_p, axis=1)
        L_hat = r_p / d[:, None]

        total_vel = (ep_vel.xyz + obs_vel.xyz).to_value(u.m/u.s).T
        radial_velocity = -np.sum(total_vel * L_hat, axis=1)

        return radial_velocity
    
    def topo2bary(self, topo_times, return_mjd=True, interp=False):
        if interp:
            bary_times = topo_times + self.barycentre_delays_interp(topo_times)
        else:
            bary_times = self.topo2bary_calc(topo_times, return_mjd=return_mjd)
        return bary_times
        
    def generate_interp(self, obs_range):
        time_samples = self.observation_span(obs_range, n_samples=10**4, return_mjd=False)
        topo_times = self.obs_start_time_TIME + TimeDelta(time_samples, format="sec")

        bary_times = self.topo2bary_calc(topo_times, return_mjd=False)
        bary_corr = bary_times - time_samples

        return interp1d(time_samples, bary_corr, kind="cubic")
    



    