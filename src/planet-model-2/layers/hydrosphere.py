from layers.layer import layer
from EOS.H2O import eos_water
from utils import modify_file_by_lines
from constants import M_EARTH, R_EARTH, EARTH_ATM, ABSOLUTE_ZERO

import numpy as np
import os

PHREEQC_path = '/data/pt426/phreeqc/bin'

class hydrosphere(layer):
    
    def calculate_CO2(self):

        input_file_modifications = {
            4 : f'temp        60.0        # Temperature in degrees Celsius [1, 2]',
            5 : f'pressure    1.5         # Pressure in atmospheres [2]'
        }

        input_file_path = 'templates/phreeqc_co2_equilibrium_input_template.txt'
        input_file_path_new = f'{PHREEQC_path}/input'

        modify_file_by_lines(input_file_path, input_file_path_new, input_file_modifications)

        wd = os.getcwd()
        os.chdir(PHREEQC_path)
        os.system('./phreeqc input')
        os.chdir(wd)


def equilbriate_CO2(P, T, pH_initial, C_dissolved_initial, P_CO2_initial):
    pass


def phreeqc_CO2_equilibrium_phase(P, T, pH, C_dissolved, P_CO2):
    
    input_template_file_path = 'templates/phreeqc_co2_equilibrium_input_template.txt'
    input_file_path_new = f'{PHREEQC_path}/input'
    output_file_path = f'{PHREEQC_path}/input.out'
    
    input_file_modifications = {
            4 : f'    temp        {T + ABSOLUTE_ZERO:.1f}        # Temperature in degrees Celsius [1, 2]',
            5 : f'    pressure    {P / EARTH_ATM:.2f}         # Pressure in atmospheres [2]',
            6 : f'    pH          {pH:.1f}         # Initial pH',
            8 : f'    C           {C_dissolved:.2f}         # Total dissolved carbon [4]',
            11 : f'    CO2(g)      {np.log10(P_CO2 / EARTH_ATM)}        10.0    # partial pressure in log(atm) and number of moles of CO2'
        }

    modify_file_by_lines(input_template_file_path, input_file_path_new, input_file_modifications)

    wd = os.getcwd()
    os.chdir(PHREEQC_path)
    os.system('./phreeqc input')
    os.chdir(wd)