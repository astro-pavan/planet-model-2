import numpy as np

from core import core
from hydrosphere import hydrosphere
from atmosphere import atmo

M_earth = 5.972e24 # kg
R_earth = 6371000 # m 

class planet:

    def __init__(self, M, R, F_star, atm_vmrs):
        
        self.mass = None
        self.radius = None
        
        self.water_mass_fraction = None
        self.hydrogen_mass_fraction = None

        self.instellation = None
        self.host_star_spectral_type = None

        self.P_surface = None
        self.T_surface = None

        self.find_R, self.find_M, self.find_x_H2O, self.find_x_H = False, False, False, False
        self.find_star_properties, self.find_surface_conditions = False, False
        self.hydrogen_envelope = False

    def set_bulk_properties(self, mass=None, radius=None, water_mass_fraction=None, hydrogen_mass_fraction=None):

        self.mass = mass
        self.radius = radius
        self.water_mass_fraction = water_mass_fraction
        self.hydrogen_mass_fraction = hydrogen_mass_fraction

        self.find_R = radius is None
        self.find_M = mass is None
        self.find_x_H2O = water_mass_fraction is None
        self.find_x_H = hydrogen_mass_fraction is None

    def set_surface_conditions(self, P_surface, T_surface):
        
        self.P_surface = P_surface
        self.T_surface = T_surface

        self.find_star_properties = True

    def set_atmosphere_properties(self):
        pass

    def set_host_star(self, instellation, spectral_type):

        self.instellation = instellation
        self.host_star_spectral_type = spectral_type

        self.find_surface_conditions = True

    def set_chemical_composition(self):
        pass

    def generate_planet(self):
        pass

    def integrate_planet(self, M, R, x_H2O, x_H):
        
        residual_R = 0

        return residual_R


    
    # print('Generating hydrosphere...')
    # self.hydrosphere = hydrosphere(mass, mass - water_mass, radius, P_surface, T_surface)

    # print('Generating atmosphere...')
    # self.atmosphere = atmo(mass, radius, P_surface, T_surface, 'H2O', instellation, star_spectral_type)


if __name__ == '__main__':

    test_planet = planet(5 * M_earth, 2 * R_earth, 0.1, 1e5, 300, 'H2O', instellation=0.8)

    test_planet.atmosphere.run_AGNI()

        
