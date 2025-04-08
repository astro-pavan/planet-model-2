from layers.layer import layer, M_Earth, R_EARTH
from EOS.H2O import eos_water
from utils import modify_file_by_lines

import os

PHREEQC_path = '/data/pt426/phreeqc/bin'

class hydrosphere(layer):
    
    def calculate_CO2(self):

        input_file_modifications = {}

        input_file_path = 'templates/phreeqc_co2_equilibrium_input_template.txt'
        input_file_path_new = f'{PHREEQC_path}/input'

        modify_file_by_lines(input_file_path, input_file_path_new, input_file_modifications)

        wd = os.getcwd()
        os.chdir(PHREEQC_path)
        os.system('./phreeqc input')
        os.chdir(wd)




if __name__ == '__main__':

    test = hydrosphere(M_Earth * 5, M_Earth * 4, R_EARTH * 2, 1e5, 300)
    test.calculate_CO2()