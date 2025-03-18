from EOS import eos

import numpy as np
import pandas as pd

class eos_hydrogen_helium(eos):

    col_names = ['logT', 'logP', 'logRho', 'logU', 'logS', 'dlrho/dlT_P', 'dlrho/dlP_T', 'dlS/dlT_P', 'dlS/dlP_T', 'grad_ad']

    df_pt = pd.read_table('planet/EOS/Tables/TABLEEOS_2021_TP_Y0275_v1', skiprows=2, sep='\s+', names=col_names)

    print(df_pt.head())

    super().__init__()


if __name__ == '__main__':

    eos = eos_hydrogen_helium()

    