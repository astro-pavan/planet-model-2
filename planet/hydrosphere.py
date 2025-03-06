from layer import layer, M_earth, R_earth
from EOS.H2O import eos_water

class hydrosphere(layer):

    def __init__(self, m_top, m_bottom, r_start, P_start, T_start, integrate_down=True, temp_profile='adiabatic'):

        eos_h2o = eos_water()
        eos_h2o.make_interpolators()
        
        super().__init__(m_top, m_bottom, r_start, P_start, T_start, eos_h2o, integrate_down, temp_profile)


if __name__ == '__main__':

    test = hydrosphere(M_earth * 5, M_earth * 4, R_earth * 2, 1e5, 300)