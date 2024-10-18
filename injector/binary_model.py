import sys
import numpy as np 
import astropy.constants as const
from scipy.optimize import fsolve
from scipy.interpolate import interp1d


class BinaryModel:
    def __init__(self, pulsar_pars, generate=True):
        period, A1, M1, M2, inc, ecc, AoP, LoAN = self.calc_binary_pars(pulsar_pars)
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

            self.a1_sini_c = A1
            self.a = self.get_semi_major(period, M1+M2) # m
            self.a1 = M2/(M1+M2) * self.a
            self.gamma =  M2**2 * (M1 + 2*M2) * self.P2pi * self.e / (self.a1 * (M1 + M2)**2)

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
        
    def orbit_par_converter(self, period, find='A1', A1=None, M1=None, M2=None, inc=None):
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
            return self.get_semi_major(period, M1+M2) * M2/(M1+M2) * sini / const.c.value
        
    def calc_binary_pars(self, pulsar_pars):
        period = pulsar_pars.get('binary_period', None)
        if not period:
            return 0, 0, 0, 0, 0, 0, 0, 0
        period_float = abs(pulsar_pars['Binary period']) * 3600

        ecc = abs(pulsar_pars['ecc'])
        aop = pulsar_pars['AoP']
        loan = pulsar_pars['LoAN']
        
        A1 = pulsar_pars['A1']
        M1 = pulsar_pars['M1']
        M2 = pulsar_pars['M2']
        inc = pulsar_pars['inc']
        if A1:
            if not (M1 and M2 and inc):
                M1_default, inc_default = 1.4, 90
                if (not M1) and (not M2) and (not inc):
                    M2 = self.orbit_par_converter(period_float, find='M2', A1=A1, M1=M1_default, inc=inc_default) 
                elif (M1) and (not M2) and (not inc):
                    M2 = self.orbit_par_converter(period_float, find='M2', A1=A1, M1=M1, inc=inc_default)                
                elif (not M1) and (M2) and (not inc):
                    M1 = self.orbit_par_converter(period_float, find='M1', A1=A1, M2=M2, inc=inc_default)                
                elif (not M1) and (not M2) and (inc):
                    M2 = self.orbit_par_converter(period_float, find='M2', A1=A1, M1=M1_default, inc=inc)                 
                elif (M1) and (M2) and (not inc):
                    inc = self.orbit_par_converter(period_float, find='inc', A1=A1, M1=M1, M2=M2)                 
                elif (M1) and (not M2) and (inc):
                    M2 = self.orbit_par_converter(period_float, find='M2', A1=A1, M1=M1, inc=inc)                 
                elif (not M1) and (M2) and (inc):
                    M1 = self.orbit_par_converter(period_float, find='M1', A1=A1, M2=M2, inc=inc)                
                else:
                    A1_calc = self.orbit_par_converter(period_float, find='A1', M1=M1, M2=M2, inc=inc) 
                    if A1_calc != A1:
                        sys.exit(f"Error: Inputed A1 for pulsar {pulsar_pars['ID']} is not compatible with inputed M1, M2 and inc. Only three out of these four parameters are required.")
        elif M1 and M2 and inc:
            A1 = self.orbit_par_converter(period_float, find='A1', M1=M1, M2=M2, inc=inc) 
        else:
            sys.exit(f"Error: Pulsar {pulsar_pars['ID']} requires either an A1 or M1, M2 and inc.")

        return period_float, abs(A1), abs(M1), abs(M2), inc, ecc, aop, loan