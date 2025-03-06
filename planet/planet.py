import numpy as np

from core import core
from hydrosphere import hydrosphere
from atmosphere import atmo

M_earth = 5.972e24
R_earth = 6371000

class planet:

    def __init__(self, mass, radius, water_mass_fraction, P_surface, T_surface, atmosphere_type):

        water_mass = water_mass_fraction * mass

        print('Generating hydrosphere...')
        self.hydrosphere = hydrosphere(mass, mass - water_mass, radius, P_surface, T_surface)

        print('Generating atmosphere...')
        self.atmosphere = atmo(mass, radius, P_surface, T_surface, 'H2O')


if __name__ == '__main__':

    test_planet = planet(5 * M_earth, 2 * R_earth, 0.1, 1e5, 300, 'H2O')

        