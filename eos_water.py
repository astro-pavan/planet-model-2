import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, RegularGridInterpolator, CubicSpline
from contourpy import contour_generator
import matplotlib.pyplot as plt

from utils import make_into_pair_array

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


A1_P = A2_P[:, 0]
A1_T = A2_T[0, :]

A2_s = np.reshape(df_pt['s'], (1093, 301))

# Define boundaries
T_boundary_3_7 = np.linspace(300, 2250)
P_boundary_3_7 = 700e9 * np.ones_like(T_boundary_3_7)

T_boundary_5_7 = np.linspace(2250, 4000)
P_boundary_5_7 = 10 ** (np.log10(42e9) - np.log10(6) * (((T_boundary_5_7 / 1000) - 2) / 18))

T_boundary_6_7 = np.linspace(4000, 30000)
P_boundary_6_7 = 0.05e9 + (3e9 - 0.05e9) * (((T_boundary_6_7 / 1000) - 1) / 39)

P_boundary_5_7_isothermal = np.linspace(P_boundary_5_7[-1], P_boundary_6_7[0])
T_boundary_5_7_isothermal = 4000 * np.ones_like(P_boundary_5_7_isothermal)

P_boundary_3_7_isothermal = np.linspace(P_boundary_3_7[-1], P_boundary_5_7[0])
T_boundary_3_7_isothermal = 2250 * np.ones_like(P_boundary_5_7_isothermal)

# Concatenate all boundary points
T_mazevet_boundaries = np.concatenate([T_boundary_3_7, T_boundary_3_7_isothermal, 
                                T_boundary_5_7, T_boundary_5_7_isothermal, T_boundary_6_7])
P_mazevet_boundaries = np.concatenate([P_boundary_3_7, P_boundary_3_7_isothermal, 
                                P_boundary_5_7, P_boundary_5_7_isothermal, P_boundary_6_7])

# Sort boundary points in ascending order of T
sorted_indices = np.argsort(T_mazevet_boundaries)
T_mazevet_boundaries = T_mazevet_boundaries[sorted_indices]
P_mazevet_boundaries = P_mazevet_boundaries[sorted_indices]

# Interpolate boundary curve
PT_mazevet_boundary_interp = interp1d(T_mazevet_boundaries, P_mazevet_boundaries, bounds_error=False, fill_value=(P_mazevet_boundaries[0], P_mazevet_boundaries[-1]))
# PT_mazevet_boundary_interp = CubicSpline(T_mazevet_boundaries, P_mazevet_boundaries, bc_type='clamped', extrapolate=True)

s_PT_interpolator = RegularGridInterpolator(
    (np.log10(A1_P), np.log10(A1_T)),
    np.log10(A2_s),
    method='linear',
    fill_value=None,
    bounds_error=True
)

plt.contourf(A1_T, A1_P, A2_phase_simple, [-2, -1, 0, 1, 2, 3, 4], cmap='Set2')
plt.colorbar()
plt.xscale('log')
plt.yscale('log')
plt.savefig('phase.png', dpi=500)
plt.close()


def s(P, T):
    return 10 ** s_PT_interpolator(make_into_pair_array(np.log10(P), np.log10(T)))


def is_in_mazevet(P, T):
    
    # Check if (T, P) is inside the region
    return P >= PT_mazevet_boundary_interp(T)


def generate_adiabat(P, T):

    s_val = s(P, T)

    s_contour_generator = contour_generator(A2_T, A2_P, A2_s)

    s_contours = s_contour_generator.lines(s_val)

    adiabat = None

    for s_contour in s_contours:

        s_contour_P = s_contour[:, 1]
        s_contour_T = s_contour[:, 0]

        P_min, P_max = s_contour_P.min(), s_contour_P.max()

        if P_max > P > P_min:
            plt.plot(s_contour_T, s_contour_P, 'k-')
            adiabat = interp1d(s_contour_P, s_contour_T)
        else:
            plt.plot(s_contour_T, s_contour_P, 'k--')

        plt.plot(T_mazevet_boundaries, P_mazevet_boundaries, 'r-')

    
    return adiabat


plt.contourf(A1_T, A1_P, A2_phase_simple, [-2, -1, 0, 1, 2, 3, 4], cmap='Set2')
plt.colorbar()

generate_adiabat(1e5, 200)
generate_adiabat(1e5, 250)
generate_adiabat(1e5, 300)
generate_adiabat(1e5, 350)
generate_adiabat(1e5, 400)

plt.yscale('log')
plt.ylim([1e13, 1e4])
plt.xlim([100, 3000])
plt.savefig('contour.png', dpi=500)
plt.close()
