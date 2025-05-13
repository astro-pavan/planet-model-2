from layers.atmosphere import atmosphere
from layers.hydrosphere import hydrosphere
from constants import ABSOLUTE_ZERO, EARTH_ATM
from utils import modify_file_by_lines

import numpy as np
import pandas as pd
import os
import subprocess

PHREEQC_path = '/data/pt426/phreeqc/bin'

class surface:

    def __init__(self, hydrosphere : hydrosphere, atmosphere : atmosphere):
        
        self.hydrosphere = hydrosphere
        self.atmosphere = atmosphere

        self.P = self.atmosphere.P[-1]
        self.T = self.atmosphere.T[-1]

        self.calculate_surface_conditions()
        #self.calculate_surface_conditions()
        #self.calculate_surface_conditions()

    def calculate_surface_conditions(self):
        self.gas_ocean_equilbrium(keep_atmosphere_constant=True)
        self.atmosphere.run_AGNI(T_iso=self.T, high_spectral_res=False)
        self.atmosphere.run_HELIOS()
        self.T = self.atmosphere.T[-1]
        print(f'{self.T:0f} K')

    def gas_ocean_equilbrium(self, keep_atmosphere_constant=True):

        x_CO2 = self.atmosphere.x_gas['CO2'][-1]
        
        if keep_atmosphere_constant:

            m_water = 1
            mol_CO2 = 1

        else:

            mmw_atm = self.atmosphere.mmw[-1]
            m_atm = self.atmosphere.column_surface_density() / self.hydrosphere.column_surface_density()
            m_water = 1

            mol_atm = m_atm / mmw_atm
            mol_CO2 = mol_atm * x_CO2

        molality_C, x_CO2, x_H2O = CO2_equilibrium(
            self.P, self.T, 
            x_CO2,
            self.hydrosphere.molality['C'],
            m_water,
            mol_CO2,
            constant_atmosphere=keep_atmosphere_constant
            )
        
        if not keep_atmosphere_constant:
            self.atmosphere.change_gas_species('CO2', x_CO2)
        
        self.atmosphere.change_gas_species('H2O', x_H2O)

        self.hydrosphere.molality['C'] = molality_C
                

    def H2O_evaporation(self):

        relative_humidity = 0.8

        Tc = self.T + ABSOLUTE_ZERO # Convert to Celsius

        # from Huang (2018)
        P_sat_vap_water = np.exp(34.494 - (4924.99/(Tc + 237.1))) / ((Tc + 105) ** 1.57) # in Pa
        P_sat_vap_ice = np.exp(43.494 - (6545.8/(Tc + 278))) / ((Tc + 868) ** 2) # in Pa

        saturation_vapour_pressure = P_sat_vap_water if Tc > 0 else P_sat_vap_ice

        P_H2O = relative_humidity * saturation_vapour_pressure
        
        x_H2O_old = self.atmosphere.x_gas['H2O']
        x_H2O_new = np.full_like(x_H2O_old, P_H2O / self.P)

        self.atmosphere.change_gas_species('H2O', x_H2O_new)

    def phreeqc_equilibrium_phase(self, keep_atmosphere_constant=False):
    
        input_template_file_path = 'templates/phreeqc_co2_equilibrium_input_template.txt'
        input_file_path_new = f'{PHREEQC_path}/input'
        output_file_path = f'{PHREEQC_path}/output.txt'

        atmosphere_column_density = self.atmosphere.column_surface_density()
        hydrosphere_column_density = self.hydrosphere.column_surface_density()

        mmw_atmosphere = self.atmosphere.mmw[0]

        mol_atm = atmosphere_column_density / mmw_atmosphere
        mol_CO2 = mol_atm * self.atmosphere.x_gas['CO2'][0]

        if mol_CO2 == 0:
            scale_factor = 10000
        else:
            scale_factor = 1 / mol_CO2

        m_water = hydrosphere_column_density * scale_factor

        P_CO2 = self.P * self.atmosphere.x_gas['CO2'][0]
        
        input_file_modifications = {
                4 : f'    temp        {self.T + ABSOLUTE_ZERO:.1f}        # Temperature in degrees Celsius',
                5 : f'    pressure    {self.P / EARTH_ATM:.2f}         # Pressure in atmospheres',
                6 : f'    pH          {self.hydrosphere.pH:.1f}         # Initial pH',
                8 : f'    Ca          {self.hydrosphere.molality["Ca"] * 1000:.6f}         # Calcium concentration',
                9 : f'    Na          {self.hydrosphere.molality["Na"] * 1000:.6f}         # Sodium concentration',
                10 : f'    Cl          {self.hydrosphere.molality["Cl"] * 1000:.6f}         # Chloride concentration',
                11 : f'    C           {self.hydrosphere.molality["C"] * 1000:.6f}         # Total dissolved carbon',
                12 : f'    water       {m_water:.2f}         # mass of water in kg',
                15 : f'    CO2(g)      {np.max([np.log10(P_CO2 / EARTH_ATM), -10])}        {mol_CO2 * scale_factor:.2f}    # partial pressure in log(atm) and number of moles of CO2',
            }

        modify_file_by_lines(input_template_file_path, input_file_path_new, input_file_modifications)

        wd = os.getcwd()
        os.chdir(PHREEQC_path)
        subprocess.run(['./phreeqc', 'input'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.chdir(wd)

        solution_df = pd.read_table(output_file_path, sep='\s+')

        self.hydrosphere.molality['C'] = solution_df.at[1, 'C']
        self.hydrosphere.molality['Ca'] = solution_df.at[1, 'Ca']
        self.hydrosphere.molality['Na'] = solution_df.at[1, 'Na']
        self.hydrosphere.molality['Cl'] = solution_df.at[1, 'Cl']

        d_CO2 = solution_df.at[1, 'd_CO2(g)']

        if not keep_atmosphere_constant:

            mol_CO2_new = mol_CO2 * (1 + d_CO2)
            x_CO2_new = mol_CO2_new / mol_atm

            self.atmosphere.change_gas_species('CO2', x_CO2_new)

        return d_CO2
    


def phreeqc_calculation(P, T, m_water, mol_CO2, x_CO2, molality_C):
    
    input_template_file_path = 'templates/phreeqc_co2_equilibrium_input_template.txt'
    input_file_path_new = f'{PHREEQC_path}/input'
    output_file_path = f'{PHREEQC_path}/output.txt'

    P_CO2 = x_CO2 * P

    input_file_modifications = {
                4 : f'    temp        {T + ABSOLUTE_ZERO:.1f}        # Temperature in degrees Celsius',
                5 : f'    pressure    {P / EARTH_ATM:.2f}         # Pressure in atmospheres',
                6 : f'    pH          {7:.1f}         # Initial pH',
                8 : f'    Ca          {0:.8f}         # Calcium concentration',
                9 : f'    Na          {0:.8f}         # Sodium concentration',
                10 : f'    Cl          {0:.8f}         # Chloride concentration',
                11 : f'    C           {molality_C:.8f}         # Total dissolved carbon',
                12 : f'    water       {m_water:.2f}         # mass of water in kg',
                15 : f'    CO2(g)      {np.max([np.log10(P_CO2 / EARTH_ATM), -10])}        {mol_CO2:.2f}    # partial pressure in log(atm) and number of moles of CO2',
            }

    modify_file_by_lines(input_template_file_path, input_file_path_new, input_file_modifications)

    wd = os.getcwd()
    os.chdir(PHREEQC_path)
    subprocess.run(['./phreeqc', 'input'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.chdir(wd)

    solution_df = pd.read_table(output_file_path, sep='\s+')

    final_molality_CO2 = solution_df.at[1, 'C']
    final_x_CO2 = (EARTH_ATM * 10 ** solution_df.at[1, 'si_CO2(g)']) / P
    final_x_H2O = (EARTH_ATM * 10 ** solution_df.at[1, 'si_H2O(g)']) / P

    return final_molality_CO2, final_x_CO2, final_x_H2O


def CO2_equilibrium(P, T, x_CO2, molality_C, m_water=1, mol_CO2=1, constant_atmosphere=False):

    molality_C_old, x_CO2_old = 0, x_CO2
    molality_C_new, x_CO2_new = molality_C, x_CO2
    start = True
    n = 0

    while np.abs(molality_C_old - molality_C_new) > 1e-9 or start:

        start = False
        
        molality_C_old = molality_C_new

        if not constant_atmosphere:
            x_CO2_old = x_CO2_new
            
        # print(f'x_CO2 : {x_CO2_old:.6f}, molality_CO2 : {molality_C_old:.6f}')

        molality_C_new, x_CO2_new, x_H2O_new = phreeqc_calculation(P, T, m_water, mol_CO2, x_CO2_old, molality_C_old)

        # print(f"\rCALCULATING HYDROSPHERE CO2 CONTENT (STEPS: {n})", end='', flush=True)
        n += 1

    # print('')

    return molality_C_new, x_CO2_new, x_H2O_new


if __name__ == '__main__':

    molality, x_CO2 = CO2_equilibrium(1e5, 300, 0.004, 0, constant_atmosphere=True)

    print(molality)
    print(x_CO2)

    # print(phreeqc_calculation(1e5, 280, 1, 100, 0.004, molality)[1])

    molality, x_CO2 = CO2_equilibrium(1e5, 295, 0.004, molality, constant_atmosphere=False)

    print(molality)
    print(x_CO2)


