import sys
import numpy as np 
import astropy.constants as const
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

from .io_tools import str2func

class PulsarBinaryModel:
    def __init__(self, pulsar_pars, validate=False):
        period, A1, M1, M2, inc, ecc, AoP, LoAN = sanitize_binary_pars(pulsar_pars)
        self.mass_p = M1 # Msun
        self.mass_c = M2 # Msun
        self.period = period # s
        if self.period:
            self.P2pi = period/(2*np.pi)

            self.e = min(ecc, 0.99)
            self.AoP = np.deg2rad(AoP) # rad
            self.LoAN = np.deg2rad(LoAN) # rad
            self.I = np.deg2rad(inc)
            self.sini = np.sin(self.I)

            self.A1 = A1
            self.a = self.get_semi_major(period, M1+M2) # m
            self.a1 = M2/(M1+M2) * self.a
            self.gamma =  M2**2 * (M1 + 2*M2) * self.P2pi * self.e / (self.a1 * (M1 + M2)**2)

        self.T0 = None
        if not validate:
            self.orbital_delay = self.generate_interp()

   
    @staticmethod
    def get_semi_major(period, mass_func):
        return  (const.G.value*(mass_func)*const.M_sun.value  * (period)**2 / (4*np.pi**2)) ** (1/3)
    
    @staticmethod
    def mass_func(period, semi_major_axis):
        return  semi_major_axis ** 3 * 4 * np.pi**2 /(const.G.value * period**2 * const.M_sun.value) 
    
    @staticmethod
    def orbit_par_converter(period, find='A1', A1=None, M1=None, M2=None, inc=None):
        T = const.G.value * period**2 / (4*np.pi**2)

        def get_M1(M2_):
            sini = np.sin(np.deg2rad(inc))
            term1 = np.sqrt(T) * ((sini * M2_*const.M_sun.value)/(const.c.value * A1))**(3/2)
            return term1/const.M_sun.value - M2_

        if find == 'inc':
            sini3 = (M1 + M2)**2/M2**3 *1/const.M_sun.value * 1/T * (const.c.value * A1) ** 3
            return np.rad2deg(np.arcsin(sini3**(1/3)))
        
        elif find == 'M1':
            return get_M1(M2)
        
        elif find == 'M2':
            return fsolve(lambda M2_: M1 - get_M1(M2_), 1)
        
        elif find == 'A1':
            sini = np.sin(np.deg2rad(inc))
            return PulsarBinaryModel.get_semi_major(period, M1+M2) * M2/(M1+M2) * sini / const.c.value
    
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
        return com * Ob * self.proj_a * ecc * (np.cos(self.AoP + theta) + self.e*np.cos(self.AoP))
    
    def get_radial_accel(self, t, star='pulsar'):
        numerator_mass = self.mass_c if star=='pulsar' else self.mass_p
        com = numerator_mass/(self.mass_p+self.mass_c) 
        Ob = 1/self.P2pi
        ecc = 1/(1-self.e**2) 
        theta = self.true_anomaly(t)
        return -com * Ob ** 2 * self.proj_a * ecc * (1 + self.e*np.cos(theta))**2 * np.sin(self.AoP + theta)

    def get_roemer_delay_proper(self, t):
        E = self.eccentric_anomaly(t)
        return self.A1 * (np.cos(E) - self.e)*np.sin(self.AoP) + self.A1*np.sin(E)*np.sqrt(1-self.e**2)*np.cos(self.AoP)
    
    def get_roemer_delay_coord(self, t):
        E = self.eccentric_anomaly(t)
        sinE, cosE = np.sin(E), np.cos(E)

        alpha = self.A1 * np.sin(self.AoP)
        beta = np.sqrt(1-self.e**2) * self.A1 * np.cos(self.AoP)

        term1 = alpha * (cosE - self.e)
        term2 = (beta + self.gamma) * sinE
        term3 = term1 + term2
        term4 = (alpha*sinE - beta*cosE)
        term5 = self.P2pi * (1 - self.e * cosE) 

        return term3 + term4 * term3 / term5
    
    def get_orbital_delay(self, t, interp=False):
        if self.period:
            if interp:
                return self.roemer_delay_interp(t % self.period)
            else:
                return self.get_roemer_delay_coord(t)
        else:
            return lambda t: np.zeros_like(t)
        
    def generate_interp(self):
        if self.period:
            time_series = np.linspace(0, self.period, int(3e4))
            roemer_delay = self.get_roemer_delay_coord(time_series)
            interp_func = interp1d(time_series, roemer_delay)
            return lambda t: interp_func(t % self.period)
        else:
            return lambda t: np.zeros_like(t)
        

def sanitize_binary_pars(pulsar_pars):
    period = pulsar_pars.get('binary_period', None)
    if not period:
        return 0, 0, 0, 0, 0, 0, 0, 0
    period_float = str2func(period, 'Binary period', pulsar_pars['ID'], float) * 3600

    ecc = str2func(pulsar_pars.get('ecc', 0), 'ecc', pulsar_pars['ID'], float)
    aop = str2func(pulsar_pars.get('aop', 0), 'aop', pulsar_pars['ID'], float)
    loan = str2func(pulsar_pars.get('laon', 0), 'loan', pulsar_pars['ID'], float)
    
    A1 = pulsar_pars.get('A1', None)
    M1 = pulsar_pars.get('M1', None)
    M2 = pulsar_pars.get('M2', None)
    inc = pulsar_pars.get('inc', None)
    if A1:
        if not (M1 and M2 and inc):
            A1_float = str2func(A1, 'A1', pulsar_pars['ID'], float) 
            M1_default, inc_default = 1.4, 90

            if (not M1) and (not M2) and (not inc):
                inc_float = inc_default
                M1_float = M1_default
                M2_float = PulsarBinaryModel.orbit_par_converter(period_float, find='M2', A1=A1_float, M1=M1_float, inc=inc_float) 
            
            elif (M1) and (not M2) and (not inc):
                inc_float = inc_default
                M1_float = str2func(M1, 'M1', pulsar_pars['ID'], float) 
                M2_float = PulsarBinaryModel.orbit_par_converter(period_float, find='M2', A1=A1_float, M1=M1_float, inc=inc_float) 
            
            elif (not M1) and (M2) and (not inc):
                inc_float = inc_default
                M2_float = str2func(M2, 'M2', pulsar_pars['ID'], float)
                M1_float = PulsarBinaryModel.orbit_par_converter(period_float, find='M1', A1=A1_float, M2=M2_float, inc=inc_float) 
            
            elif (not M1) and (not M2) and (inc):
                inc_float = str2func(inc, 'inc', pulsar_pars['ID'], float)
                M1_float = M1_default
                M2_float = PulsarBinaryModel.orbit_par_converter(period_float, find='M2', A1=A1_float, M1=M1_float, inc=inc_float) 
            
            elif (M1) and (M2) and (not inc):
                M1_float = str2func(M1, 'M1', pulsar_pars['ID'], float) 
                M2_float = str2func(M2, 'M2', pulsar_pars['ID'], float) 
                inc_float = PulsarBinaryModel.orbit_par_converter(period_float, find='inc', A1=A1_float, M1=M1_float, M2=M2_float) 
            
            elif (M1) and (not M2) and (inc):
                M1_float = str2func(M1, 'M1', pulsar_pars['ID'], float) 
                inc_float = str2func(inc, 'inc', pulsar_pars['ID'], float)
                M2_float = PulsarBinaryModel.orbit_par_converter(period_float, find='M2', A1=A1_float, M1=M1_float, inc=inc_float) 
            
            elif (not M1) and (M2) and (inc):
                M2_float = str2func(M2, 'M2', pulsar_pars['ID'], float) 
                inc_float = str2func(inc, 'inc', pulsar_pars['ID'], float)
                M1_float = PulsarBinaryModel.orbit_par_converter(period_float, find='M1', A1=A1_float, M2=M2_float, inc=inc_float) 
            
            else:
                M1_float = str2func(M1, 'M1', pulsar_pars['ID'], float) 
                M2_float = str2func(M2, 'M2', pulsar_pars['ID'], float) 
                inc_float = str2func(inc, 'inc', pulsar_pars['ID'], float)
                A1_calc = PulsarBinaryModel.orbit_par_converter(period_float, find='A1', M1=M1_float, M2=M2_float, inc=inc_float) 
                if A1_calc != A1_float:
                    sys.exit(f"Error: Inputed A1 for pulsar {pulsar_pars['ID']} is not compatible with inputed M1, M2 and inc. Only three out of these four parameters are required.")
            
            
    elif (M1) and (M2) and (inc):
        M1_float = str2func(M1, 'M1', pulsar_pars['ID'], float) 
        M2_float = str2func(M2, 'M2', pulsar_pars['ID'], float) 
        inc_float = str2func(inc, 'inc', pulsar_pars['ID'], float)
        A1_float = PulsarBinaryModel.orbit_par_converter(period_float, find='A1', M1=M1_float, M2=M2_float, inc=inc_float) 
    else:
        sys.exit(f"Error: Pulsar {pulsar_pars['ID']} requires either an A1 or M1, M2 and inc.")

    return abs(period_float), abs(A1_float), abs(M1_float), abs(M2_float), inc_float, abs(ecc), aop, loan