from layers.atmosphere import atmosphere
from layers.hydrosphere import hydrosphere
from constants import ABSOLUTE_ZERO, EARTH_ATM
from utils import modify_file_by_lines

import numpy as np
import os

PHREEQC_path = '/data/pt426/phreeqc/bin'

class surface:

    def __init__(self, hydrosphere : hydrosphere, atmosphere : atmosphere):
        
        self.hydrosphere = hydrosphere
        self.atmosphere = atmosphere

        self.P = self.atmosphere.P[-1]
        self.T = self.atmosphere.T[-1]

    def H2O_evaporation(self):

        relative_humidity = 0.8

        P_H2O = relative_humidity * saturation_vapour_pressure(self.T)

        print(P_H2O / 1e2)
        
        x_H2O_old = self.atmosphere.x_gas['H2O']
        x_H2O_new = np.full_like(x_H2O_old, P_H2O / self.P)

        self.atmosphere.change_gas_species('H2O', x_H2O_new)

        print(self.atmosphere.x_gas)

    def gas_ocean_interaction(self):

        pass


def saturation_vapour_pressure(T):
    
    Tc = T + ABSOLUTE_ZERO # Convert to Celsius

    # from Huang (2018)
    P_sat_vap_water = np.exp(34.494 - (4924.99/(Tc + 237.1))) / ((Tc + 105) ** 1.57) # in Pa
    P_sat_vap_ice = np.exp(43.494 - (6545.8/(Tc + 278))) / ((Tc + 868) ** 2) # in Pa

    return np.where(Tc > 0, P_sat_vap_water, P_sat_vap_ice)



def phreeqc_CO2_equilibrium_phase(P, T, pH, C_dissolved, P_CO2, m_water=1, mol_CO2=10):
    
    input_template_file_path = 'templates/phreeqc_co2_equilibrium_input_template.txt'
    input_file_path_new = f'{PHREEQC_path}/input'
    output_file_path = f'{PHREEQC_path}/input.out'
    
    input_file_modifications = {
            4 : f'    temp        {T + ABSOLUTE_ZERO:.1f}        # Temperature in degrees Celsius',
            5 : f'    pressure    {P / EARTH_ATM:.2f}         # Pressure in atmospheres',
            6 : f'    pH          {pH:.1f}         # Initial pH',
            11 : f'    C           {C_dissolved:.2f}         # Total dissolved carbon',
            12 : f'    water       {m_water:.2f}         # mass of water in kg',
            15 : f'    CO2(g)      {np.log10(P_CO2 / EARTH_ATM)}        {mol_CO2:.2f}    # partial pressure in log(atm) and number of moles of CO2',
        }

    modify_file_by_lines(input_template_file_path, input_file_path_new, input_file_modifications)

    wd = os.getcwd()
    os.chdir(PHREEQC_path)
    os.system('./phreeqc input')
    os.chdir(wd)


    