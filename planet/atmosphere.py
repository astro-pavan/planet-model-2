from layer import layer, M_earth, R_earth
from EOS.H2O import eos_water

class atmo(layer):

    def __init__(self, m_bottom, r_bottom, P_surface, T_surface, atmosphere_type, P_min=1):
        
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

        super().__init__(m_bottom * 1.5, m_bottom, r_bottom, P_surface, T_surface, self.eos, integrate_down=False, temp_profile='isothermal')

        P_mask = self.P > P_min

        self.m = self.m[P_mask]
        self.r = self.r[P_mask]
        self.P = self.P[P_mask]
        self.T = self.T[P_mask]
        self.rho = self.rho[P_mask]