from csaps import CubicSmoothingSpline
import numpy as np


def smooth_parameter(cop, x):
    dx = np.mean(np.diff(x))  # mean distance between two points
    lamda = (1 / (2 * cop * np.pi)) ** 4 / dx  # eg. 8 in Enting1987
    p = 1 / (1 + lamda)  # close to eq. 13 in Enting1987
    return p


def smooth_data(cop, ydata2smooth, xdata, new_grid):
    p = smooth_parameter(cop, xdata)
    sp = CubicSmoothingSpline(xdata, ydata2smooth, smooth=p)
    y_smooth = sp(new_grid)
    return y_smooth, new_grid

