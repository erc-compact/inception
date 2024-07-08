import numpy as np 
import astropy.constants as const
from scipy.optimize import fsolve
from scipy.interpolate import interp1d


class PulsarBinaryModel:
    def __init__(self, pulsar_mass, companion_mass, period, 
                 eccentricity, inclination, LoAN, AoP):
        self.mass_p = pulsar_mass # Msun
        self.mass_c = companion_mass # Msun
        self.period = period * 3600 # s
        self.a = self.get_semi_major() # m
        self.proj_a = self.a * np.sin(np.deg2rad(inclination))
        self.x = self.proj_a/const.c.value
        self.e = min(eccentricity, 0.999)
        self.I = np.deg2rad(inclination) # rad
        self.LoAN = np.deg2rad(LoAN) # rad
        self.AoP = np.deg2rad(AoP) # rad

        # orbit interpolation
        timeaxis = np.linspace(0, self.period, 10000)
        orbit_delay = self.get_orbital_delay(timeaxis, interp=False)
        self.interp_orbit = interp1d(timeaxis, orbit_delay, kind='cubic')
        
    def get_semi_major(self):
        return  (const.G.value*(self.mass_p+self.mass_c)*const.M_sun.value  * (self.period)**2 / (4*np.pi**2)) ** (1/3)
    
    def radial(self, theta):
        return (1-self.e**2)/(1+self.e*np.cos(theta)) * self.a
    
    def get_coords(self, theta):
        r = self.radial(theta)
        x = r * (np.cos(self.LoAN)*np.cos(self.AoP+theta) - np.sin(self.LoAN)*np.sin(self.AoP+theta)*np.cos(self.I))
        y = r * (np.sin(self.LoAN)*np.cos(self.AoP+theta) + np.cos(self.LoAN)*np.sin(self.AoP+theta)*np.cos(self.I))
        z = r * np.sin(self.AoP+theta)*np.sin(self.I)
        return x, y, z
    
    def star_coord(self, theta, star='pulsar'):
        rx, ry, rz = self.get_coords(theta)
        numerator_mass = -self.mass_c if star=='pulsar' else self.mass_p
        com = numerator_mass/(self.mass_p+self.mass_c) 
        return rx*com, ry*com, rz*com

    def eccentric_anomaly(self, t, time_scale=1, ref_time=0): 
        period = self.period * time_scale
        ecc = self.e
   
        def mean_anomaly(t):
            return 2*np.pi/period * (t-ref_time)

        def root_eccentric_anomaly(E, t):
            return E - ecc*np.sin(E) - mean_anomaly(t)
        
        @np.vectorize
        def find_eccentric_anomaly(t):
            return fsolve(root_eccentric_anomaly, x0=1, args=(t))
        
        return find_eccentric_anomaly(t)
    
    def true_anomaly(self, t, time_scale=1, ref_time=0): 
        period = self.period * time_scale
        ecc = self.e
        tanE2 = np.tan(self.eccentric_anomaly(t, time_scale, ref_time)/2)
        tan_ecc = np.sqrt((1+ecc)/(1-ecc)) * tanE2
        TA = 2*np.arctan(tan_ecc)

        cond_list = [(t-ref_time < period/2), (t-ref_time >= period/2)]
        func_list = [TA,  2*np.pi + TA]
        return np.select(cond_list, func_list)
    
    def get_radial_velocity(self, t, star='pulsar'):
        numerator_mass = self.mass_c if star=='pulsar' else self.mass_p
        com = numerator_mass/(self.mass_p+self.mass_c) 
        Ob = 2*np.pi/self.period 
        ecc = 1/np.sqrt(1-self.e**2) 
        theta = self.true_anomaly(t)
        return com * Ob * self.proj_a * ecc * (np.cos(self.AoP + theta) + self.e*np.cos(self.AoP))
    
    def get_radial_accel(self, t, star='pulsar'):
        numerator_mass = self.mass_c if star=='pulsar' else self.mass_p
        com = numerator_mass/(self.mass_p+self.mass_c) 
        Ob = 2*np.pi/self.period 
        ecc = 1/(1-self.e**2) 
        theta = self.true_anomaly(t)
        return -com * Ob ** 2 * self.proj_a * ecc * (1 + self.e*np.cos(theta))**2 * np.sin(self.AoP + theta)
    
    def get_roemer_delay(self, t, interp=True, ref_time=0):
        if interp:
            return self.interp_orbit((t-ref_time)%self.period)
        else:
            return -self.star_coord(self.true_anomaly(t-ref_time))[2]/const.c.value