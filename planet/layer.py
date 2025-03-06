from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import numpy as np

G = 6.674e-11
M_earth = 5.972e24
R_earth = 6371000

class layer:

    def __init__(self, m_top, m_bottom, r_start, P_start, T_start, eos, integrate_down=True, temp_profile='adiabatic'):
        
        n = 64
        self.m = np.linspace(m_bottom, m_top, n)
        self.r, self.P, self.T, self.rho = np.zeros_like(self.m), np.zeros_like(self.m), np.zeros_like(self.m), np.zeros_like(self.m)
        self.eos = eos

        self.max_mass_step = 0.01 * M_earth

        t_span = (m_top, m_bottom) if integrate_down else (m_bottom, m_top)

        y0 = [r_start, P_start]

        if temp_profile == 'adiabatic':
            self.T_profile = self.eos.generate_adiabat(P_start, T_start)
        elif temp_profile == 'isothermal':
            self.T_profile = lambda P: T_start
        
        def dr_dm(m, r, P):
            rho = self.eos.rho_PT(P, self.T_profile(P))[0]
            return 1 / (4 * np.pi * (r ** 2) * rho)
        
        def dP_dm(m, r, P):
            return - (G * m) / (4 * np.pi * (r ** 4))
        
        def f(t, y):
            dr = dr_dm(t, y[0], y[1])
            dP = dP_dm(t, y[0], y[1])
            return (dr, dP)
        
        solution = solve_ivp(f, t_span, y0, max_step=self.max_mass_step)

        m_solution = solution.t
        r_solution = solution.y[0, :]
        P_solution = solution.y[1, :]
        positive_mask = P_solution > 0

        m_solution = m_solution[positive_mask]
        r_solution = r_solution[positive_mask]
        log_P_solution = np.log10(P_solution[positive_mask])

        try:
            r_interpolator = CubicSpline(m_solution, r_solution)
            log_P_interpolator = CubicSpline(m_solution, log_P_solution)
        except ValueError:
            r_interpolator = CubicSpline(m_solution[::-1], r_solution[::-1])
            log_P_interpolator = CubicSpline(m_solution[::-1], log_P_solution[::-1])

        self.r = r_interpolator(self.m)
        self.P = 10 ** log_P_interpolator(self.m)
        self.T = self.T_profile(self.P) if temp_profile == 'adiabatic' else np.full_like(self.m, T_start)

        print(self.P)
        print(self.T)

        self.rho = self.eos.rho_PT(self.P, self.T)


if __name__ == '__main__':

    from EOS.H2O import eos_water

    eos_h2o = eos_water()
    eos_h2o.make_interpolators()

    water_layer = layer(M_earth, 0.8 * M_earth, R_earth, 1e5, 300, eos_h2o)
