from layer import layer, M_earth, R_earth, G
from EOS.H2O import eos_water
from utils import modify_file_by_lines

import netCDF4
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import os

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

class atmo(layer):

    def __init__(self, m_bottom, r_bottom, P_surface, T_surface, atmosphere_type, instellation=1, spectral_type='G2', P_min=1):
        
        self.eos = None

        if atmosphere_type == 'HHe':
            pass
        elif atmosphere_type == 'H2O':
            self.eos = eos_water()
        elif atmosphere_type == 'CO2':
            pass
        elif atmosphere_type == 'N2':
            pass

        self.eos.make_interpolators()
        self.max_mass_step = None

        super().__init__(m_bottom * 1.1, m_bottom, r_bottom, P_surface, T_surface, self.eos, integrate_down=False, temp_profile='isothermal')

        P_mask = self.P > P_min

        self.m = self.m[P_mask]
        self.r = self.r[P_mask]
        self.P = self.P[P_mask]
        self.T = self.T[P_mask]
        self.rho = self.rho[P_mask]

        self.instellation = instellation
        self.host_star_spectral_type = spectral_type

        self.P_surface = self.P[0]
        self.r_surface = self.r[0]
        self.g_surface = (G * self.m[0]) / (self.r_surface ** 2)

    def radiative_transfer(self):
        
        P_initial = np.logspace(1, 5, num=25)

        P, T = self.run_AGNI(350, return_PT=True)

        new_P = np.logspace(1, 5, num=60)
        new_T = CubicSpline(P, T)(new_P)

        self.run_AGNI((new_P, new_T), high_spectral_res=True)

    def run_AGNI(self, PT_initial, high_spectral_res=False, return_PT=False):

        wd = os.getcwd()
        config_file_path = f'{AGNI_path}/res/config/default.toml'
        config_file_path_new = f'{AGNI_path}/res/config/pl2.toml'

        n_spectral_bands = 256 if high_spectral_res else 48

        config_file_modifications = {
            5 : 'title = "pl2"',
            55 : '    solution_type   = 3                         # Solution type (see wiki).',
            30 : f'    p_surf          = {self.P_surface / 1e5:.2f}                     # Total surface pressure [bar].',
            15 : f'    radius          = {self.r_surface:.3e}            # Planet radius at the surface [m].',
            16 : f'    gravity         = {self.g_surface:.2e}              # Gravitational acceleration at the surface [m s-2]',
            9 : f'    instellation    = {self.instellation * solar_constant:.1f}           # Stellar flux at planet\'s orbital distance [W m-2].',
            26 : f'    input_star      = "res/stellar_spectra/{spectral_types[self.host_star_spectral_type]}"              # Path to stellar spectrum.',
            56 : f'    solvers         = ["levenberg"]                        # Ordered list of solvers to apply (see wiki).',
            25 : f'    input_sf        = "res/spectral_files/Dayspring/{n_spectral_bands}/Dayspring.sf"   # Path to SOCRATES spectral file.',
            32 : '    vmr_dict        = { H2=1.0 }               # Volatile volume mixing ratios (=mole fractions).',
            60 : '    easy_start      = true                     # Initially down-scale convective/condensation fluxes, if initial guess is poor.',
            11 : '    s0_fact         = 1.0            # Stellar flux scale factor which accounts for planetary rotation (c.f. Cronin+13).',
            14 : '    albedo_s        = 0.5               # Grey surface albedo when material=greybody.'
        }

        if type(PT_initial) is float or type(PT_initial) is int:

            T_initial = PT_initial

            config_file_modifications[58] = f'    initial_state   = ["iso", "{T_initial:.0f}"]     # Ordered list of requests describing the initial state of the atmosphere (see wiki).'
            n_levels = 25
            T_surface = T_initial

        else:

            P_initial, T_initial = PT_initial

            n_levels = len(P_initial)
            PT = pd.DataFrame({'P' : P_initial, 'T' : T_initial})
            PT_profile_file_path = f'{AGNI_path}/res/config/PT_initial.csv'

            PT.to_csv(PT_profile_file_path, index=False, header=False)

            config_file_modifications[58] = f'    initial_state   = ["csv", "{PT_profile_file_path}"]     # Ordered list of requests describing the initial state of the atmosphere (see wiki).'
            T_surface = T_initial[-1]


        config_file_modifications[8] = f'    tmp_surf        = {T_surface:.1f}            # Surface temperature [kelvin].'
        config_file_modifications[43] = f'    num_levels      = {n_levels:.0f}                       # Number of model levels.'

        modify_file_by_lines(config_file_path, config_file_path_new, config_file_modifications)

        os.chdir(AGNI_path)
        os.system('./agni.jl res/config/pl2.toml')
        os.chdir(wd)

        AGNI_output_file_path = f'{AGNI_path}/out/atm.nc'

        atm_nc = netCDF4.Dataset(AGNI_output_file_path)

        # print(atm_nc.variables.keys())

        P = np.array(atm_nc['p'])
        T = np.array(atm_nc['tmp'])

        atm_nc.close()

        if return_PT:
            return (P, T)

