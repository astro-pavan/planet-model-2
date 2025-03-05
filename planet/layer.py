from scipy.integrate import solve_ivp
import numpy as np

G = 6.674e-11
M_earth = 5.972e24
R_earth = 6371000

class layer:

    def __init__(self, m_top, m_bottom, r_start, P_start, T_start, eos, integrate_down=True, temp_profile='adiabat'):

        n = 64

        t_span = (m_top, m_bottom) if integrate_down else (m_bottom, m_top)
        t_span = (m_top, m_bottom)
        t_eval = np.linspace(t_span[0], t_span[1], n, dtype=float)

        y0 = [r_start, P_start]

        if temp_profile == 'adiabat':
            self.T_profile = eos.generate_adiabat(P_start, T_start)
        elif temp_profile == 'isotherm':
            self.T_profile = lambda P: T_start

        def rho_eos(P):
            return eos.rho_PT(P, self.T_profile(P))
        
        def dr_dm(m, r, P):
            rho = rho_eos(P)[0]
            return 1 / (4 * np.pi * (r ** 2) * rho)
        
        def dP_dm(m, r, P):
            return - (G * m) / (4 * np.pi * (r ** 4))
        
        def f(t, y):
            dr = dr_dm(t, y[0], y[1])
            dP = dP_dm(t, y[0], y[1])
            return (dr, dP)

        print('Generating layer...')
        
        solution = solve_ivp(f, t_span, y0, max_step=0.01*M_earth)

        print('Layer generated')


if __name__ == '__main__':

    from EOS.H2O import eos_water

    eos_h2o = eos_water()
    eos_h2o.make_interpolators()

    water_layer = layer(M_earth, 0.8 * M_earth, R_earth, 1e5, 300, eos_h2o)
        
