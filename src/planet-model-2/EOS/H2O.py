from EOS.EOS import eos

import numpy as np
import pandas as pd

class eos_water(eos):

    def __init__(self):

        super().__init__()

        # Define column names based on the table header
        columns = [
            "press", "temp", "rho", "ad_grad", "s", "u", "c", "mmw", "x_ion", "x_d", "phase"
        ]

        # Specify column data types
        dtype_mapping = {
            "press": float, "temp": float, "rho": float, "ad_grad": float,
            "s": float, "u": float, "c": float, "mmw": float, "x_ion": float, "x_d": float, "phase": int
        }

        # Read the data, skipping metadata and header lines
        file_path_pt = "EOS-data/aqua_eos_pt_v1_0.dat"
        df_pt = pd.read_csv(file_path_pt, sep='\s+', names=columns, skiprows=19, engine='python', dtype=dtype_mapping)

        A2_P = np.reshape(df_pt['press'], (1093, 301))
        A2_T = np.reshape(df_pt['temp'], (1093, 301))
        A2_phase = np.reshape(df_pt['phase'], (1093, 301))

        A2_phase_simple = np.where(A2_phase == 5, -1, 3) # supercritical = -1
        A2_phase_simple = np.where(A2_phase == 3, 0, A2_phase_simple) # vapour = 0
        A2_phase_simple = np.where(A2_phase == 4, 1, A2_phase_simple) # liquid = 1
        A2_phase_simple = np.where(A2_phase == -1, 2, A2_phase_simple) # ice (low) = 2
        A2_phase_simple = np.where((A2_phase == -7) | (A2_phase == -10), 4, A2_phase_simple) # ice (high) = 4, ice (mid) = 3

        self.A2_phase = A2_phase_simple
        self.phase_names = {
            -1 : 'supercritical/superionic',
            0 : 'vapour',
            1 : 'liquid',
            2 : 'ice (I)',
            3 : 'ice (II, III, V, VI)',
            4 : 'ice (VII, X)'}
        
        self.A1_P = A2_P[:, 0]
        self.A1_T = A2_T[0, :]

        self.A2_rho_PT = np.reshape(df_pt['rho'], (1093, 301))
        self.A2_s_PT = np.reshape(df_pt['s'], (1093, 301))
        self.A2_u_PT = np.reshape(df_pt['u'], (1093, 301))

        super().make_interpolators()