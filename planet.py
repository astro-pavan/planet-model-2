# file with planet class

import os
import numpy as np
import pandas as pd

import seagen
import woma

import matplotlib.pyplot as plt

magrathea_path = '/data/pt426/Magrathea'

class planet:

    def __init__(self, name, core_mass, mantle_mass, ocean_mass, atmosphere_mass, surface_temp):

        # rewrites config file for MARGRATHEA and deletes prior output file

        config_file_path = f'{magrathea_path}/run/mode0.cfg'
        wd = os.getcwd()

        try:
            os.remove(f'{magrathea_path}/result/{name}.txt')
        except FileNotFoundError:
            pass

        # NOTE : I've changed the water and atmosphere EOS to the tabulated ones.

        config_file_modifications = {
            12: f'mass_of_core={core_mass:.2f}	# Earth Masses in core',
            13: f'mass_of_mantle={mantle_mass:.2f}	# Earth Masses in mantle',
            14: f'mass_of_hydro={ocean_mass:.2f}# Earth Masses in hydrosphere',
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

        self.df = pd.read_table(output_file_path, sep='	 ', engine='python')

        self.M = np.max(self.df['M (earth)']) * 5.9724e24
        self.R = np.max(self.df['Radius (earth)']) * 6371000

    def make_particle_planet(self, n_particles):
        
        # makes MAGRATHEA data into a form SEAGEN likes

        # gets material IDs for each material

        mat_id = np.where(self.df['Phase'] == 'Fe hcp (Smith)', 401, 0)
        mat_id = np.where(self.df['Phase'] == 'Si PPv (Sakai)', 400, mat_id)
        mat_id = np.where(self.df['Phase'] == 'Pv (Dorogokupets)', 400, mat_id)
        mat_id = np.where(self.df['Phase'] == 'Rwd (Dorogokupets)', 400, mat_id)
        mat_id = np.where(self.df['Phase'] == 'Wds (Dorogokupets)', 400, mat_id)
        mat_id = np.where(self.df['Phase'] == 'Fo/Ol (Dorogokupets)', 400, mat_id)
        mat_id = np.where(self.df['Phase'] == 'H2O (AQUA)', 304, mat_id)
        mat_id = np.where(self.df['Phase'] == 'H/He (Chabrier)', 307, mat_id)
        mat_id = np.where(self.df['Phase'] == 'Isothermal Ideal Gas', 307, mat_id)

        woma.load_eos_tables()

        # calculates values in SI units
        
        A1_r = np.array(self.df['Radius (earth)'] * 6371000)
        A1_rho = np.array(self.df['Density (g cm^-3)'] * 1000)
        A1_mat_id = mat_id
        A1_T = np.array(self.df['T (K)'])
        A1_P = np.array(self.df['P (GPa)'] * 1e9)
        A1_u = woma.A1_u_rho_T(A1_rho, A1_T, A1_mat_id)

        # deletes elements in the array where r doesn't increase

        mask = (A1_r - np.roll(A1_r, 1)) != 0

        A1_r = np.array(A1_r[mask])
        A1_rho = np.array(A1_rho[mask])
        A1_mat_id = np.array(A1_mat_id[mask])
        A1_T = np.array(A1_T[mask])
        A1_P = np.array(A1_P[mask])
        A1_u = np.array(A1_u[mask])

        # uses SEAGEN to make particle planet

        particles = seagen.GenSphere(
            n_particles,
            A1_r[1:],
            A1_rho[1:],
            A1_mat_id[1:],
            A1_u[1:],
            A1_T[1:],
            A1_P[1:],
        )

        # calculates smoothing length

        A1_h = np.cbrt(48 * particles.A1_m / (4 / 3 * np.pi * particles.A1_rho)) / 2

        A2_pos = np.transpose([particles.A1_x, particles.A1_y, particles.A1_z])

        return particles, A1_h, A2_pos
    