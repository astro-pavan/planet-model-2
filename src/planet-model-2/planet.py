import numpy as np
import os
import pandas as pd

from layers.core import core
from layers.hydrosphere import hydrosphere
from layers.atmosphere import atmosphere

from utils import modify_file_by_lines
from constants import M_EARTH, R_EARTH

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

        print('GENERATING ATMOSPHERE...')

        self.atmosphere = atmosphere(self.radius, self.mass, atm_vmrs, self.instellation, spec_type, self.P_surface, T_initial)

        print('ATMOSPHERE GENERATED')

        self.T_surface = self.atmosphere.T[-1]

        print('GENERATING PLANET...')

        df_core, df_hydro = self.run_Magarathea(self.T_surface)

        m_core = np.array(df_core['M (earth)']) * M_EARTH
        r_core = np.array(df_core['Radius (earth)']) * R_EARTH
        P_core = np.array(df_core['P (GPa)']) / 1e9
        T_core = np.array(df_core['T (K)'])  
        rho_core =  np.array(df_core['Density (g cm^-3)']) * 1e3

        m_hydro = np.array(df_hydro['M (earth)']) * M_EARTH
        r_hydro = np.array(df_hydro['Radius (earth)']) * R_EARTH
        P_hydro = np.array(df_hydro['P (GPa)']) / 1e9
        T_hydro = np.array(df_hydro['T (K)'])
        rho_hydro =  np.array(df_hydro['Density (g cm^-3)']) * 1e3

        self.core = core(m_core, r_core, P_core, T_core, rho_core)
        self.hydrosphere = hydrosphere(m_hydro, r_hydro, P_hydro, T_hydro, rho_hydro)

        print('PLANET GENERATED')


    def run_Magarathea(self, T_surface):
        
        mode0_config_file_path = f'{magrathea_path}/run/mode0.cfg'
        mode4_config_file_path = f'{magrathea_path}/run/mode4.cfg'
        input_file_path = f'{magrathea_path}/input/MR.txt'
        mode4_output_file_path = f'{magrathea_path}/result/outputplanetsols.txt'
        mode0_output_file_path = f'{magrathea_path}/result/pl2.txt'

        wd = os.getcwd()

        with open(input_file_path, 'w') as file:
            file.writelines(['M (Earth-masses) 	R (Earth-radii)\n', f'{self.mass / M_EARTH:.2f} {self.radius / R_EARTH:.2f}'])

        mode4_config_file_modifications = {
            37 : f'P_surface={self.P_surface * 10:.1e}				# The pressure level that the broad band optical transit radius probes (in microbar)',
            29 : f'surface_temp={T_surface:.0f}	# Kelvin, top of planet where enclosed mass equals total mass',
            13 : f'input_file="./input/MR.txt"		# Input file name & location'
        }

        modify_file_by_lines(mode4_config_file_path, mode4_config_file_path, mode4_config_file_modifications)

        os.chdir(magrathea_path)
        os.system('./planet run/mode4.cfg')

        mode4_results_table = pd.read_table(mode4_output_file_path, sep='\s+')

        R_core_ocean_boundary = mode4_results_table.at[0, "RMantle"]

        mode0_config_file_modifications = {
            12: f'mass_of_core={mode4_results_table.at[0, "MCore"]:.2f}	# Earth Masses in core',
            13: f'mass_of_mantle={mode4_results_table.at[0, "MMantle"]:.2f}	# Earth Masses in mantle',
            14: f'mass_of_hydro={mode4_results_table.at[0, "MWater"]:.2f} # Earth Masses in hydrosphere',
            15: f'mass_of_atm=0		# Earth Masses in atmosphere',
            16: f'surface_temp={T_surface:.0f}	# Kelvin, top of planet where enclosed mass equals total mass',
            21: f'output_file="./result/pl2.txt"	# Output file name & location',
            25: f'P_surface={self.P_surface * 10:.1e}			# The pressure level that the broad band optical transit radius probes (in microbar)'
        }

        modify_file_by_lines(mode0_config_file_path, mode0_config_file_path, mode0_config_file_modifications)

        os.system('./planet run/mode0.cfg')
        os.chdir(wd)

        df = pd.read_table(mode0_output_file_path, sep='	 ', engine='python')
        
        for col in df.columns:
            if col != 'Phase':
                df[col] = pd.to_numeric(df[col], errors='coerce')

        core_mask = df['Radius (earth)'] < R_core_ocean_boundary

        df_core = df[core_mask]
        df_hydro = df[~core_mask]

        return df_core, df_hydro


if __name__ == '__main__':

    test_planet = planet(1.1 * M_EARTH, 1.1 * R_EARTH, 1, 1e5, 300, 'G2', '')

        
