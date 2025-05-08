from layers.layer import layer
from utils import modify_file_by_lines
from constants import G, IDEAL_GAS_CONSTANT

import netCDF4
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import os
import subprocess
import matplotlib.pyplot as plt

AGNI_path = '/home/pt426/AGNI'

spectral_types = {
    'M8' : 'trappist-1.txt',
    'M4.5': 'gj1214.txt',
    'M4' : 'gj1132.txt',
    'M3.5' : 'gj849.txt',
    'M3' : 'l-98-59.txt',
    'M2' : 'gj179.txt',
    'K6' : 'hd85512.txt',
    'K2' : 'eps-eri.txt',
    'K1' : 'hd97658.txt',
    'G2' : 'sun.txt',
    'G0' : 'HIP67522.txt'
}

solar_constant = 1360 # W / m^2

class atmosphere(layer):

    def __init__(self, r_bottom, m_bottom, atm_vmrs, F_star, spec_type, P_surface, T_surface_initial, tidally_locked):

        self.m_bottom, self.r_bottom = m_bottom, r_bottom

        self.g_surface = (G * self.m_bottom) / (self.r_bottom ** 2)

        self.instellation = F_star
        self.host_star_spectral_type = spec_type
        self.P_surface = P_surface
        self.tidally_locked = tidally_locked

        self.run_AGNI(T_iso=T_surface_initial, x_gas=atm_vmrs, high_spectral_res=True)

    def run_AGNI(self, T_iso=None, x_gas=None, high_spectral_res=False, n_levels=25):

        print('FINDING RADIATIVE CONVECTIVE EQULIBRIUM...')

        wd = os.getcwd()
        config_file_path = 'templates/default.toml'
        config_file_path_new = f'{AGNI_path}/res/config/pl2.toml'

        n_spectral_bands = 256 if high_spectral_res else 48

        if x_gas is None:
            x_gas = self.x_gas

        x_N2 = x_gas['N2'][0]
        x_CO2 = x_gas['CO2'][0]
        x_H2O = x_gas['H2O'][0]

        vmr_dict = '{' + f'N2={x_N2:.6f}, CO2={x_CO2:.6f}, H2O={x_H2O:.6f}' + '}'

        config_file_modifications = {
            30 : f'    p_surf          = {self.P_surface / 1e5:.2f}                     # Total surface pressure [bar].',
            15 : f'    radius          = {self.r_bottom:.3e}            # Planet radius at the surface [m].',
            16 : f'    gravity         = {self.g_surface:.2e}              # Gravitational acceleration at the surface [m s-2]',
            9 : f'    instellation    = {self.instellation * solar_constant:.1f}           # Stellar flux at planet\'s orbital distance [W m-2].',
            26 : f'    input_star      = "res/stellar_spectra/{spectral_types[self.host_star_spectral_type]}"              # Path to stellar spectrum.',
            25 : f'    input_sf        = "res/spectral_files/Dayspring/{n_spectral_bands}/Dayspring.sf"   # Path to SOCRATES spectral file.',
            32 : f'    vmr_dict        = {vmr_dict}               # Volatile volume mixing ratios (=mole fractions).',
        }

        if self.tidally_locked:
            config_file_modifications[11] = '    s0_fact         = 1.0               # Stellar flux scale factor which accounts for planetary rotation (c.f. Cronin+13).'

        if T_iso is not None:

            config_file_modifications[58] = f'    initial_state   = ["iso", "{T_iso:.0f}"]     # Ordered list of requests describing the initial state of the atmosphere (see wiki).'
            T_surface = T_iso

        else:

            if n_levels != len(self.P):
                P_interpolator = CubicSpline(self.P, self.T)
                self.P = np.logspace(1, 5, num=n_levels)
                self.T = P_interpolator(self.P)

            n_levels = len(self.P)
            PT = pd.DataFrame({'P' : self.P, 'T' : self.P})
            PT_profile_file_path = f'{AGNI_path}/res/config/PT_initial.csv'

            PT.to_csv(PT_profile_file_path, index=False, header=False)

            config_file_modifications[58] = f'    initial_state   = ["csv", "{PT_profile_file_path}"]     # Ordered list of requests describing the initial state of the atmosphere (see wiki).'
            T_surface = self.T[-1]


        config_file_modifications[8] = f'    tmp_surf        = {T_surface:.1f}            # Surface temperature [kelvin].'
        config_file_modifications[43] = f'    num_levels      = {n_levels:.0f}                       # Number of model levels.'

        modify_file_by_lines(config_file_path, config_file_path_new, config_file_modifications)

        os.chdir(AGNI_path)
        subprocess.run(['./agni.jl', 'res/config/pl2.toml'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        os.chdir(wd)

        AGNI_output_file_path = f'{AGNI_path}/out/atm.nc'

        atm_nc = netCDF4.Dataset(AGNI_output_file_path)

        # print(atm_nc.variables.keys())

        x_gas_nc = np.array(atm_nc['x_gas'])
        x_gas = {}

        for i, l in enumerate(atm_nc['gases']):
            species = b''.join(l).decode('utf-8').rstrip()
            x_gas[species] = x_gas_nc[:, i]

        P = np.array(atm_nc['p'])
        T = np.array(atm_nc['tmp'])
        z = np.array(atm_nc['z'])
        mmw = np.array(atm_nc['mmw'])

        atm_nc.close()

        self.P = P
        self.T = T
        self.r = z + self.r_bottom

        self.rho = (P * mmw) / (IDEAL_GAS_CONSTANT * T)

        dr = - np.diff(self.r, prepend=self.r[0])
        dm = (4 * np.pi * (self.r ** 2) * self.rho) * dr
        self.m = (np.sum(dm) - np.cumsum(dm)) + self.m_bottom

        self.x_gas = x_gas
        self.mmw = mmw

        plt.plot(z, self.rho)
        plt.savefig('rhoz.png')
        plt.close()

        plt.plot(z, T)
        plt.savefig('Tz.png')
        plt.close()

        plt.plot(np.log10(z), T)
        plt.savefig('Tlogz.png')
        plt.close()

        plt.plot(z, np.log10(T))
        plt.savefig('logTz.png')
        plt.close()

        plt.plot(np.log10(z), np.log10(T))
        plt.savefig('logTlogz.png')
        plt.close()

        print('RADIATIVE CONVECTIVE EQULIBRIUM FOUND')
    
    def change_gas_species(self, modified_species, x_new):

        x_old = self.x_gas[modified_species]

        if type(x_new) is not np.ndarray:
            x_new = np.full_like(x_old, x_new)

        self.x_gas[modified_species] = x_new
        
        for species in self.x_gas.keys():
            if species != modified_species:
                self.x_gas[species] = self.x_gas[species] / (1 + (x_new - x_old))
    


