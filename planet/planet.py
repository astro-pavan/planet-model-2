import os
import pandas as pd
import numpy as np

from core import core
from hydrosphere import hydrosphere
from atmosphere import atmosphere

magrathea_path = '/data/pt426/Magrathea'

class planet:

    def __init__(self, name, core_mass, iron_fraction, ocean_mass, atmosphere_mass, surface_temp):

        # rewrites config file for MARGRATHEA and deletes prior output file

        config_file_path = f'{magrathea_path}/run/mode0.cfg'
        wd = os.getcwd()

        try:
            os.remove(f'{magrathea_path}/result/{name}.txt')
        except FileNotFoundError:
            pass

        # NOTE : I've changed the water and atmosphere EOS to the tabulated ones.

        config_file_modifications = {
            12: f'mass_of_core={core_mass * iron_fraction:.2f}	# Earth Masses in core',
            13: f'mass_of_mantle={core_mass * (1 - iron_fraction):.2f}	# Earth Masses in mantle',
            14: f'mass_of_hydro={ocean_mass:.2f} # Earth Masses in hydrosphere',
            15: f'mass_of_atm={atmosphere_mass:.2f}		# Earth Masses in atmosphere',
            16: f'surface_temp={surface_temp:.0f}	# Kelvin, top of planet where enclosed mass equals total mass',
            21: f'output_file="./result/{name}.txt"	# Output file name & location'
        }

        try:
            # Read the file into a list of lines
            with open(config_file_path, 'r') as file:
                lines = file.readlines()
            
            # Apply modifications
            for line_number, new_content in config_file_modifications.items():
                if 1 <= line_number <= len(lines):
                    lines[line_number - 1] = new_content + '\n'
                else:
                    print(f"Line {line_number} is out of range. Skipping.")
            
            # Write the modified lines back to the file
            with open(config_file_path, 'w') as file:
                file.writelines(lines)

            print("File modified successfully.")
        
        except FileNotFoundError:
            print(f"Error: The file at '{config_file_path}' was not found.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
                
        # runs MAGRATHEA

        os.chdir(magrathea_path)
        os.system('./planet run/mode0.cfg')
        os.chdir(wd)

        # reads output file and loads into a dataframe

        output_file_path = f'{magrathea_path}/result/{name}.txt'

        df = pd.read_table(output_file_path, sep='	 ', engine='python')

        mat_id = np.where(self.df['Phase'] == 'Fe hcp (Smith)', 401, 0)
        mat_id = np.where(self.df['Phase'] == 'Si PPv (Sakai)', 400, self.mat_id)
        mat_id = np.where(self.df['Phase'] == 'Pv (Dorogokupets)', 400, self.mat_id)
        mat_id = np.where(self.df['Phase'] == 'Rwd (Dorogokupets)', 400, self.mat_id)
        mat_id = np.where(self.df['Phase'] == 'Wds (Dorogokupets)', 400, self.mat_id)
        mat_id = np.where(self.df['Phase'] == 'Fo/Ol (Dorogokupets)', 400, self.mat_id)
        mat_id = np.where(self.df['Phase'] == 'H2O (AQUA)', 304, self.mat_id)
        mat_id = np.where(self.df['Phase'] == 'H/He (Chabrier)', 307, self.mat_id)
        mat_id = np.where(self.df['Phase'] == 'Isothermal Ideal Gas', 307, self.mat_id)

        