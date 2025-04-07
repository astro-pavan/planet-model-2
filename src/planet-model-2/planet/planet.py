import numpy as np

from planet.core import core
from planet.hydrosphere import hydrosphere
from planet.atmosphere import atmosphere

M_earth = 5.972e24 # kg
R_earth = 6371000 # m

magrathea_path = '/data/pt426/Magrathea'

class planet:

    def __init__(self, M, R, F_star, P_surface, T_initial, spec_type, atm_vmrs):
        
        self.mass = M
        self.radius = R
        
        self.water_mass_fraction = None
        self.hydrogen_mass_fraction = None

        self.instellation = F_star
        self.host_star_spectral_type = spec_type

        self.P_surface = P_surface

        self.atmosphere = atmosphere(self.radius, self.mass, atm_vmrs, self.instellation, self.P_surface, T_initial)  

    def generate_planet(self):
        pass

    def integrate_planet(self, M, R, x_H2O):
        
        residual_R = 0

        return residual_R


    
    # print('Generating hydrosphere...')
    # self.hydrosphere = hydrosphere(mass, mass - water_mass, radius, P_surface, T_surface)

    # print('Generating atmosphere...')
    # self.atmosphere = atmo(mass, radius, P_surface, T_surface, 'H2O', instellation, star_spectral_type)


if __name__ == '__main__':

    test_planet = planet(1, 1.1, 1, 1e5, 300, 'G2', '')

        
