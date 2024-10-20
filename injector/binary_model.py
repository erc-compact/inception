import numpy as np 
import astropy.constants as const
from scipy.optimize import fsolve
from scipy.interpolate import interp1d


class BinaryModel:
    def __init__(self, pulsar_pars, generate=True):
        self.mass_p = pulsar_pars['M1'] 
        self.mass_c = pulsar_pars['M2'] 
        self.period = pulsar_pars['binary_period'] 
        if self.period:
            self.P2pi = self.period/(2*np.pi)

            self.e = min(pulsar_pars['ecc'], 0.99)
            self.AoP = np.deg2rad(pulsar_pars['AoP']) 
            self.LoAN = np.deg2rad(pulsar_pars['LoAN']) 
            self.I = np.deg2rad(pulsar_pars['inc'])
            self.sini = np.sin(self.I)

            self.a1_sini_c = pulsar_pars['A1'] 
            self.a = self.get_semi_major(self.period, self.mass_p+self.mass_c) # m
            self.a1 = self.mass_c/(self.mass_p+self.mass_c) * self.a
            self.gamma = self.mass_c**2 * (self.mass_p + 2*self.mass_c) * self.P2pi * self.e / (self.a1 * (self.mass_p + self.mass_c)**2)

        self.T0 = None
        if generate:
            self.orbital_delay = self.generate_interp()

    @staticmethod
    def get_semi_major(period, mass_func):
        return  (const.G.value*(mass_func)*const.M_sun.value  * (period)**2 / (4*np.pi**2)) ** (1/3)
    
    @staticmethod
    def mass_func(period, semi_major_axis):
        return  semi_major_axis ** 3 * 4 * np.pi**2 /(const.G.value * period**2 * const.M_sun.value) 
    
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

    def eccentric_anomaly(self, t, time_scale=1): 
        period = self.period * time_scale
   
        def mean_anomaly(t):
            return 2*np.pi/period * t

        def root_eccentric_anomaly(E, t):
            return E - self.e*np.sin(E) - mean_anomaly(t)
        
        @np.vectorize
        def find_eccentric_anomaly(t):
            return fsolve(root_eccentric_anomaly, x0=1, args=(t))
        
        return find_eccentric_anomaly(t)
    
    def true_anomaly(self, t, time_scale=1): 
        period = self.period * time_scale
        ecc = self.e
        tanE2 = np.tan(self.eccentric_anomaly(t, time_scale)/2)
        tan_ecc = np.sqrt((1+ecc)/(1-ecc)) * tanE2
        TA = 2*np.arctan(tan_ecc)

        cond_list = [(t < period/2), (t>= period/2)]
        func_list = [TA,  2*np.pi + TA]
        return np.select(cond_list, func_list)
    
    def get_radial_velocity(self, t, star='pulsar'):
        numerator_mass = self.mass_c if star=='pulsar' else self.mass_p
        com = numerator_mass/(self.mass_p+self.mass_c) 
        Ob = 1/self.P2pi
        ecc = 1/np.sqrt(1-self.e**2) 
        theta = self.true_anomaly(t)
        return com * Ob * self.a1 * ecc * (np.cos(self.AoP + theta) + self.e*np.cos(self.AoP))
    
    def get_radial_accel(self, t, star='pulsar'):
        numerator_mass = self.mass_c if star=='pulsar' else self.mass_p
        com = numerator_mass/(self.mass_p+self.mass_c) 
        Ob = 1/self.P2pi
        ecc = 1/(1-self.e**2) 
        theta = self.true_anomaly(t)
        return -com * Ob ** 2 * self.a1 * ecc * (1 + self.e*np.cos(theta))**2 * np.sin(self.AoP + theta)

    def get_roemer_delay_proper(self, t):
        E = self.eccentric_anomaly(t)
        return self.a1_sini_c * (np.cos(E) - self.e)*np.sin(self.AoP) + self.a1_sini_c*np.sin(E)*np.sqrt(1-self.e**2)*np.cos(self.AoP)
    
    def get_roemer_delay_coord(self, t):
        E = self.eccentric_anomaly(t)
        sinE, cosE = np.sin(E), np.cos(E)

        alpha = self.a1_sini_c * np.sin(self.AoP)
        beta = np.sqrt(1-self.e**2) * self.a1_sini_c * np.cos(self.AoP)

        term1 = alpha * (cosE - self.e)
        term2 = (beta + self.gamma) * sinE
        term3 = term1 + term2
        term4 = (alpha*sinE - beta*cosE)
        term5 = self.P2pi * (1 - self.e * cosE) 

        return term3 + term4 * term3 / term5
    
    def generate_interp(self):
        if self.period:
            time_series = np.linspace(0, self.period, int(5e4))
            roemer_delay = self.get_roemer_delay_coord(time_series)
            interp_func = interp1d(time_series, roemer_delay)
            return lambda t: interp_func(t % self.period)
        else:
            return lambda t: np.zeros_like(t)
        
    