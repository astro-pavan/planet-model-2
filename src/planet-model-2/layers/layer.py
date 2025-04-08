from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import numpy as np

G = 6.674e-11 # m^3 kg^-1 s^-2
M_earth = 5.972e24 # kg
R_earth = 6371000 # m 
ideal_gas_constant = 8.314 # J K^-1 mol^-1

class layer:

    def __init__(self, m, r, P, T, rho, eos=None, T_profile=None):

        self.m = m
        self.r = r
        self.P = P
        self.T = T
        self.rho = rho

        self.eos = eos
        self.T_profile = T_profile

        self.r_bottom, self.r_top = np.min(self.r), np.max(self.r)
        self.m_bottom, self.m_top = np.max(self.m), np.min(self.m)
        self.layer_mass = self.m_top - self.m_bottom

    @classmethod
    def from_hydrostatic_equilibrium(cls, m_top, m_bottom, r_start, P_start, T_start, eos, integrate_down=True, temp_profile='adiabatic'):
        
        n = 64
        m = np.linspace(m_bottom, m_top, n)
        r, P, T, rho = np.zeros_like(m), np.zeros_like(m), np.zeros_like(m), np.zeros_like(m)

        max_mass_step = 0.01 * M_earth

        t_span = (m_top, m_bottom) if integrate_down else (m_bottom, m_top)

        y0 = [r_start, P_start]

        if temp_profile == 'adiabatic':
            T_profile = eos.generate_adiabat(P_start, T_start)
        elif temp_profile == 'isothermal':
            T_profile = lambda P: T_start
        
        def dr_dm(m, r, P):
            rho = eos.rho_PT(P, T_profile(P))[0]
            return 1 / (4 * np.pi * (r ** 2) * rho)
        
        def dP_dm(m, r, P):
            return - (G * m) / (4 * np.pi * (r ** 4))
        
        def f(t, y):
            dr = dr_dm(t, y[0], y[1])
            dP = dP_dm(t, y[0], y[1])
            return (dr, dP)
        
        solution = solve_ivp(f, t_span, y0, max_step=max_mass_step)

        m_solution = solution.t
        r_solution = solution.y[0, :]
        P_solution = solution.y[1, :]
        positive_mask = P_solution > 0

        m_solution = m_solution[positive_mask]
        r_solution = r_solution[positive_mask]
        log_P_solution = np.log10(P_solution[positive_mask])

        if np.all(positive_mask):

            try:
                r_interpolator = CubicSpline(m_solution, r_solution)
                log_P_interpolator = CubicSpline(m_solution, log_P_solution)
            except ValueError:
                r_interpolator = CubicSpline(m_solution[::-1], r_solution[::-1])
                log_P_interpolator = CubicSpline(m_solution[::-1], log_P_solution[::-1])

            r = r_interpolator(m)
            P = 10 ** log_P_interpolator(m)

        else:
            
            m = m_solution
            r = r_solution
            P = P_solution[positive_mask]

        T = T_profile(P) if temp_profile == 'adiabatic' else np.full_like(m, T_start)
        rho = eos.rho_PT(P, T)

        return cls(m, r, P, T, rho, eos, T_profile)


if __name__ == '__main__':

    from EOS.H2O import eos_water

    eos_h2o = eos_water()
    eos_h2o.make_interpolators()

    water_layer = layer(M_earth, 0.8 * M_earth, R_earth, 1e5, 300, eos_h2o)
