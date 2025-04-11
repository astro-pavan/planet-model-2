import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import subprocess
import h5py

from layers.core import core
from layers.hydrosphere import hydrosphere
from layers.atmosphere import atmosphere
from surface import surface

from utils import modify_file_by_lines
from constants import M_EARTH, R_EARTH

magrathea_path = '/data/pt426/Magrathea'

class planet:

    def __init__(self, M, R, F_star, P_surface, T_initial, spec_type, atm_vmrs, tidally_locked=True):
        
        self.mass = M
        self.radius = R

        print('GENERATING ATMOSPHERE...')

        self.atmosphere = atmosphere(self.radius, self.mass, atm_vmrs, F_star, spec_type, P_surface, T_initial, tidally_locked)

        print('ATMOSPHERE GENERATED')

        print('GENERATING INTERNAL STRUCTURE...')

        df_core, df_hydro = self.run_Magarathea(P_surface, self.atmosphere.T[-1])

        m_hydro = np.array(df_hydro['M (earth)'])[::-1] * M_EARTH
        r_hydro = np.array(df_hydro['Radius (earth)'])[::-1] * R_EARTH
        P_hydro = np.array(df_hydro['P (GPa)'])[::-1] * 1e9
        T_hydro = np.array(df_hydro['T (K)'])[::-1]
        rho_hydro =  np.array(df_hydro['Density (g cm^-3)'])[::-1] * 1e3

        m_core = np.array(df_core['M (earth)'])[::-1] * M_EARTH
        r_core = np.array(df_core['Radius (earth)'])[::-1] * R_EARTH
        P_core = np.array(df_core['P (GPa)'])[::-1] * 1e9
        T_core = np.array(df_core['T (K)'])[::-1]  
        rho_core =  np.array(df_core['Density (g cm^-3)'])[::-1] * 1e3

        self.hydrosphere = hydrosphere(m_hydro, r_hydro, P_hydro, T_hydro, rho_hydro)
        self.core = core(m_core, r_core, P_core, T_core, rho_core)

        print('INTERNAL STRUCTURE GENERATED')

        self.surface = surface(self.hydrosphere, self.atmosphere)


    def run_Magarathea(self, P_surface, T_surface):
        
        mode0_config_file_path = f'{magrathea_path}/run/mode0.cfg'
        mode4_config_file_path = f'{magrathea_path}/run/mode4.cfg'
        input_file_path = f'{magrathea_path}/input/MR.txt'
        mode4_output_file_path = f'{magrathea_path}/result/outputplanetsols.txt'
        mode0_output_file_path = f'{magrathea_path}/result/pl2.txt'

        try:
            os.remove(mode0_output_file_path)
        except FileNotFoundError:
            pass

        wd = os.getcwd()

        with open(input_file_path, 'w') as file:
            file.writelines(['M (Earth-masses) 	R (Earth-radii)\n', f'{self.mass / M_EARTH:.2f} {self.radius / R_EARTH:.2f}'])

        mode4_config_file_modifications = {
            37 : f'P_surface={P_surface * 10:.1e}				# The pressure level that the broad band optical transit radius probes (in microbar)',
            29 : f'surface_temp={T_surface:.0f}	# Kelvin, top of planet where enclosed mass equals total mass',
            13 : f'input_file="./input/MR.txt"		# Input file name & location'
        }

        modify_file_by_lines(mode4_config_file_path, mode4_config_file_path, mode4_config_file_modifications)

        os.chdir(magrathea_path)
        subprocess.run(['./planet', 'run/mode4.cfg'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        mode4_results_table = pd.read_table(mode4_output_file_path, sep='\s+')

        R_core_ocean_boundary = mode4_results_table.at[0, "RMantle"]

        mode0_config_file_modifications = {
            12: f'mass_of_core={mode4_results_table.at[0, "MCore"]:.4f}	# Earth Masses in core',
            13: f'mass_of_mantle={mode4_results_table.at[0, "MMantle"]:.4f}	# Earth Masses in mantle',
            14: f'mass_of_hydro={mode4_results_table.at[0, "MWater"]:.4f} # Earth Masses in hydrosphere',
            15: f'mass_of_atm=0		# Earth Masses in atmosphere',
            16: f'surface_temp={T_surface:.0f}	# Kelvin, top of planet where enclosed mass equals total mass',
            21: f'output_file="./result/pl2.txt"	# Output file name & location',
            25: f'P_surface={P_surface * 10:.1e}			# The pressure level that the broad band optical transit radius probes (in microbar)'
        }

        modify_file_by_lines(mode0_config_file_path, mode0_config_file_path, mode0_config_file_modifications)

        subprocess.run(['./planet', 'run/mode0.cfg'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.chdir(wd)

        df = pd.read_table(mode0_output_file_path, sep='	 ', engine='python')
        
        for col in df.columns:
            if col != 'Phase':
                df[col] = pd.to_numeric(df[col], errors='coerce')

        core_mask = df['Radius (earth)'] < R_core_ocean_boundary

        df_core = df[core_mask]
        df_hydro = df[~core_mask]

        return df_core, df_hydro
    
    def plot_PT(self):
        plt.plot(self.atmosphere.T, self.atmosphere.P)
        plt.plot(self.hydrosphere.T, self.hydrosphere.P)
        plt.plot(self.core.T, self.core.P)

        plt.yscale('log')
        plt.savefig('PT.png')
        plt.close()

    def plot_PT_ocean(self):
        plt.plot(self.hydrosphere.T, self.hydrosphere.r - np.max(self.hydrosphere.r))

        # plt.yscale('log')
        plt.savefig('PT_ocean.png')
        plt.close()

    def save_to_hdf5(self, filename):
        
        with h5py.File(filename, 'w') as f:
            
            f.attrs['Mass'] = self.mass
            f.attrs['Radius'] = self.radius

            f.attrs['Instellation'] = self.atmosphere.instellation
            f.attrs['Host Star Spectral Type'] = self.atmosphere.host_star_spectral_type

            core_grp = f.create_group('Core')
            hydro_grp = f.create_group('Hydrosphere')
            atm_grp = f.create_group('Atmosphere')

            core_grp.create_dataset('m', data=self.core.m)
            core_grp.create_dataset('r', data=self.core.r)
            core_grp.create_dataset('P', data=self.core.P)
            core_grp.create_dataset('T', data=self.core.T)
            core_grp.create_dataset('rho', data=self.core.rho)

            hydro_grp.create_dataset('m', data=self.hydrosphere.m)
            hydro_grp.create_dataset('r', data=self.hydrosphere.r)
            hydro_grp.create_dataset('P', data=self.hydrosphere.P)
            hydro_grp.create_dataset('T', data=self.hydrosphere.T)
            hydro_grp.create_dataset('rho', data=self.hydrosphere.rho)

            molarity_grp = hydro_grp.create_group('Molarity')
            for species in self.hydrosphere.molarity.keys():
                molarity_grp.attrs[species] = self.hydrosphere.molarity[species]

            atm_grp.create_dataset('m', data=self.atmosphere.m)
            atm_grp.create_dataset('r', data=self.atmosphere.r)
            atm_grp.create_dataset('P', data=self.atmosphere.P)
            atm_grp.create_dataset('T', data=self.atmosphere.T)
            atm_grp.create_dataset('rho', data=self.atmosphere.rho)
            atm_grp.create_dataset('mmw', data=self.atmosphere.mmw)
            
            x_gas_grp = atm_grp.create_group('x_gas')
            for gas in self.atmosphere.x_gas.keys():
                x_gas_grp.create_dataset(gas, self.atmosphere.x_gas[gas])




if __name__ == '__main__':

    test_planet = planet(1 * M_EARTH, 1 * R_EARTH, 1, 1e5, 270, 'G2', 
                         {'N2' : [0.78], 'O2' : [0.21], 'Ar' : [0.009], 'CO2' : [0.0004], 'H2O' : [0.00]},
                           tidally_locked=False)

        
