from layers.layer import layer
from EOS.H2O import eos_water
from utils import modify_file_by_lines
from constants import M_EARTH, R_EARTH, EARTH_ATM, ABSOLUTE_ZERO

import numpy as np
import os

PHREEQC_path = '/data/pt426/phreeqc/bin'

class hydrosphere(layer):

    def __init__(self, m, r, P, T, rho, eos=None, T_profile=None):
        
        super().__init__(m, r, P, T, rho, eos, T_profile)

        self.pH = 7
        self.molarity = {
            'C' : 0,
            'Ca' : 0,
            'Na' : 0,
            'Cl' : 0
        }