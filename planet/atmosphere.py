from layer import layer, M_earth, R_earth, G
from EOS.H2O import eos_water
from scipy.io import netcdf

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


    def run_AGNI(self):

        wd = os.getcwd()
        config_file_path = f'{AGNI_path}/res/config/default.toml'
        config_file_path_new = f'{AGNI_path}/res/config/pl2.toml'

        P_surface = self.P[0] / 1e5
        T_surface = self.T[0]
        r_surface = self.r[0]
        m_surface = self.m[0]

        g_surface = (G * m_surface) / (r_surface ** 2)

        config_file_modifications = {
            5 : 'title = "pl2"',
            55 : '    solution_type   = 3                         # Solution type (see wiki).',
            30 : f'    p_surf          = {P_surface:.2f}                     # Total surface pressure [bar].',
            8 : f'    tmp_surf        = {T_surface:.1f}            # Surface temperature [kelvin].',
            15 : f'    radius          = {r_surface:.3e}            # Planet radius at the surface [m].',
            16 : f'    gravity         = {g_surface:.2e}              # Gravitational acceleration at the surface [m s-2]',
            9 : f'    instellation    = {self.instellation * solar_constant:.1f}           # Stellar flux at planet\'s orbital distance [W m-2].',
            26 : f'    input_star      = "res/stellar_spectra/{spectral_types[self.host_star_spectral_type]}"              # Path to stellar spectrum.',
            #58 : f'    initial_state   = ["iso", "300"]     # Ordered list of requests describing the initial state of the atmosphere (see wiki).',
            56 : f'    solvers         = ["newton"]                        # Ordered list of solvers to apply (see wiki).'
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
            with open(config_file_path_new, 'w') as file:
                file.writelines(lines)

            print("File modified successfully.")
        
        except FileNotFoundError:
            print(f"Error: The file at '{config_file_path}' was not found.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        os.chdir(AGNI_path)
        os.system('./agni.jl res/config/pl2.toml')
        os.chdir(wd)
