import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator
from contourpy import contour_generator


def make_into_pair_array(arr1, arr2):
    # turns two multidimensional numpy arrays into a form that can be used with the scipy interpolator

    arr1, arr2 = np.nan_to_num(arr1), np.nan_to_num(arr2)

    if type(arr1) is np.ndarray and type(arr2) is np.ndarray:

        if arr1.ndim == 0:
            return np.array([arr1[()], arr2[0]])
        if arr2.ndim == 0:
            return np.array([arr1[0], arr2[()]])

        try:
            assert np.all(arr1.shape == arr2.shape)
        except AssertionError:
            print(f'arr1 = {arr1} \n arr1 shape = {arr1.shape}')
            print(f'arr2 = {arr2} \n arr2 shape = {arr2.shape}')
            assert np.all(arr1.shape == arr2.shape)

        assert arr1.ndim == 1 or arr1.ndim == 2

        arr = np.array([arr1, arr2])

        if arr1.ndim == 1:
            return np.transpose(arr, axes=(1, 0))
        elif arr1.ndim == 2:
            return np.transpose(arr, axes=(1, 2, 0))

    else:

        if type(arr1) is np.ndarray:
            if arr1.ndim == 1:
                arr1 = arr1[0]

        if type(arr2) is np.ndarray:
            if arr2.ndim == 1:
                arr2 = arr2[0]

        return np.array([arr1, arr2])


class eos:

    def __init__(self):
        
        self.A1_P, self.A1_T = None, None

        self.phase_names = {}
        
        self.A2_rho_PT = None
        self.A2_s_PT = None
        self.A2_u_PT = None
        self.A2_phase = None

        self.A1_rho = None
        self.A2_P_rhoT = None
        self.A2_s_rhoT = None

        self.rho_PT_interpolator = None
        self.s_PT_interpolator = None
        self.u_PT_interpolator = None

    def make_interpolators(self):
        
        self.rho_PT_interpolator = RegularGridInterpolator(
            (np.log10(self.A1_P), np.log10(self.A1_T)),
            np.log10(self.A2_rho_PT),
            method='linear',
            fill_value=None,
            bounds_error=True
        )

        self.s_PT_interpolator = RegularGridInterpolator(
            (np.log10(self.A1_P), np.log10(self.A1_T)),
            np.log10(self.A2_s_PT),
            method='linear',
            fill_value=None,
            bounds_error=True
        )

        self.u_PT_interpolator = RegularGridInterpolator(
            (np.log10(self.A1_P), np.log10(self.A1_T)),
            np.log10(self.A2_u_PT),
            method='linear',
            fill_value=None,
            bounds_error=True
        )

    def rho_PT(self, P, T):
        return 10 ** self.rho_PT_interpolator(make_into_pair_array(np.log10(P), np.log10(T)))
    
    def s_PT(self, P, T):
        return 10 ** self.s_PT_interpolator(make_into_pair_array(np.log10(P), np.log10(T)))
    
    def u_PT(self, P, T):
        return 10 ** self.u_PT_interpolator(make_into_pair_array(np.log10(P), np.log10(T)))

    def generate_adiabat(self, P, T):

        s_val = self.s_PT(P, T)

        A2_T, A2_P = np.meshgrid(self.A1_T, self.A1_P)

        s_contour_generator = contour_generator(A2_P, A2_T, self.A2_s_PT)

        s_contours = s_contour_generator.lines(s_val)

        adiabat = None

        for s_contour in s_contours:

            s_contour_P = s_contour[:, 1]
            s_contour_T = s_contour[:, 0]

            P_min, P_max = s_contour_P.min(), s_contour_P.max()

            if P_max > P > P_min:
                adiabat = CubicSpline(s_contour_P, s_contour_T)

        return adiabat