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

magrathea_path = 'external/Magrathea'

class planet:

    def __init__(self, M, R, F_star, P_surface, T_initial, spec_type, atm_vmrs, tidally_locked=True):
        
        self.mass = M
        self.radius = R

        print('Generating atmosphere...')

        self.atmosphere = atmosphere(self.radius, self.mass, atm_vmrs, F_star, spec_type, P_surface, T_initial, tidally_locked)

        print('Atmosphere genereated')

        print('Generating internal structure...')

        df_core, df_hydro = self.run_Magarathea(P_surface, self.atmosphere.T[-1])

        m_hydro = np.array(df_hydro['M (earth)']) * M_EARTH
        r_hydro = np.array(df_hydro['Radius (earth)']) * R_EARTH
        P_hydro = np.array(df_hydro['P (GPa)']) * 1e9
        T_hydro = np.array(df_hydro['T (K)'])
        rho_hydro =  np.array(df_hydro['Density (g cm^-3)']) * 1e3

        m_core = np.array(df_core['M (earth)']) * M_EARTH
        r_core = np.array(df_core['Radius (earth)']) * R_EARTH
        P_core = np.array(df_core['P (GPa)']) * 1e9
        T_core = np.array(df_core['T (K)'])  
        rho_core =  np.array(df_core['Density (g cm^-3)']) * 1e3

        self.hydrosphere = hydrosphere(m_hydro, r_hydro, P_hydro, T_hydro, rho_hydro)
        self.core = core(m_core, r_core, P_core, T_core, rho_core)

        print('Internal structure generated')

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

        # os.chdir(magrathea_path)
        subprocess.run(['./planet', 'run/mode4.cfg'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=magrathea_path)

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

        subprocess.run(['./planet', 'run/mode0.cfg'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=magrathea_path)
        # os.chdir(wd)

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

            molality_grp = hydro_grp.create_group('Molarity')
            for species in self.hydrosphere.molality.keys():
                molality_grp.attrs[species] = self.hydrosphere.molality[species]

            atm_grp.create_dataset('m', data=self.atmosphere.m)
            atm_grp.create_dataset('r', data=self.atmosphere.r)
            atm_grp.create_dataset('P', data=self.atmosphere.P)
            atm_grp.create_dataset('T', data=self.atmosphere.T)
            atm_grp.create_dataset('rho', data=self.atmosphere.rho)
            atm_grp.create_dataset('mmw', data=self.atmosphere.mmw)
            
            x_gas_grp = atm_grp.create_group('x_gas')
            for gas in self.atmosphere.x_gas.keys():
                x_gas_grp.create_dataset(gas, self.atmosphere.x_gas[gas])

    @classmethod
    def load_from_hdf5(cls, filename):
        # Create an uninitialized Planet instance
        obj = cls.__new__(cls)
        
        # Create uninitialized subcomponents without calling their constructors
        obj.core = type('core', (), {})()
        obj.hydrosphere = type('hydrosphere', (), {})()
        obj.atmosphere = type('atmosphere', (), {})()

        with h5py.File(filename, 'r') as f:
            # Top-level attributes
            obj.mass = f.attrs['Mass']
            obj.radius = f.attrs['Radius']
            obj.atmosphere.instellation = f.attrs['Instellation']
            obj.atmosphere.host_star_spectral_type = f.attrs['Host Star Spectral Type']

            # Core
            core_grp = f['Core']
            obj.core.m = core_grp['m'][()]
            obj.core.r = core_grp['r'][()]
            obj.core.P = core_grp['P'][()]
            obj.core.T = core_grp['T'][()]
            obj.core.rho = core_grp['rho'][()]

            # Hydrosphere
            hydro_grp = f['Hydrosphere']
            obj.hydrosphere.m = hydro_grp['m'][()]
            obj.hydrosphere.r = hydro_grp['r'][()]
            obj.hydrosphere.P = hydro_grp['P'][()]
            obj.hydrosphere.T = hydro_grp['T'][()]
            obj.hydrosphere.rho = hydro_grp['rho'][()]

            molarity_grp = hydro_grp['Molarity']
            obj.hydrosphere.molality = {
                species: molarity_grp.attrs[species]
                for species in molarity_grp.attrs
            }

            # Atmosphere
            atm_grp = f['Atmosphere']
            obj.atmosphere.m = atm_grp['m'][()]
            obj.atmosphere.r = atm_grp['r'][()]
            obj.atmosphere.P = atm_grp['P'][()]
            obj.atmosphere.T = atm_grp['T'][()]
            obj.atmosphere.rho = atm_grp['rho'][()]
            obj.atmosphere.mmw = atm_grp['mmw'][()]

            x_gas_grp = atm_grp['x_gas']
            obj.atmosphere.x_gas = {
                gas: x_gas_grp[gas][()]
                for gas in x_gas_grp.keys()
            }

        return obj


if __name__ == '__main__':

    test_planet = planet(
        1 * M_EARTH,
        1 * R_EARTH, 
        1.0, 
        1e5, 
        280, 
        'G2', 
        {'N2' : [0.977], 'CO2' : [0.003], 'H2O' : [0.02]},
        tidally_locked=False
        )

    # test_planet.save_to_hdf5('output/planet.hdf5')
    test_planet.surface.calculate_surface_conditions(keep_atmosphere_constant=False)
    test_planet.surface.check_carbon()
    # test_planet.
    # print(f'Surface temperature: {test_planet.surface.T:.0f} K')



        
