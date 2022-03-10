import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import csaps
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'

cop1 = 1 / 200.

# Read the data --------------------------------------------------------------------------------------------------------
# d18O and corresponding ice_age
df = pd.read_excel(data_path, sheet_name='Sheet7')
depth_data = np.flipud(np.array(df[df.columns[1]]))
d18O_data = np.flipud(np.array(df[df.columns[2]])) / 1000.
ice_age = np.flipud(np.array(df[df.columns[0]])) * (-1)

# temperature and accumulation with corresponding depth
df2 = pd.read_excel(data_path, sheet_name='Sheet4')
temp = np.flipud(np.array(df2[df2.columns[4]]))
acc = np.flipud(np.array(df2[df2.columns[3]]))
depth_temp = np.flipud(np.array(df2[df2.columns[0]]))

# ice_age and corresponding depth
df3 = pd.read_excel(data_path, sheet_name='Sheet5')
ice_age_scales = np.flipud(np.array(df3[df3.columns[3]])) * (-1)
depth_scales = np.flipud(np.array(df3[df3.columns[0]]))


def smooth_parameter(cop, age_):
    dx = np.mean(np.diff(age_))  # mean distance between two points
    lamda = (1 / (2 * cop * np.pi)) ** 4 / dx  # eg. 8 in Enting1987
    p = 1 / (1 + lamda)  # close to eq. 13 in Enting1987
    return p


def cubic_func(x, a, b, c):
    return a * x ** 2 + b * x + c


def exponential_func(x, a_, b_, c_, d_):
    return a_ * np.exp(b_ * x**2 + c_ * x + d_)


# Get corresponding temperature values to d18O -------------------------------------------------------------------------
def interpolate_d18O_ice_age(depth_scales, ice_age_scales, ice_age, d18O_data, depth_temp):
    # interpolate depth - ice_age
    DIA = interpolate.interp1d(depth_scales, ice_age_scales, 'linear', fill_value='extrapolate')
    ice_age_T = DIA(depth_temp)

    # interpolate d18O and ice_age
    IAd18O = interpolate.interp1d(ice_age, d18O_data, 'linear', fill_value='extrapolate')
    d18O_T = IAd18O(ice_age_T)
    return d18O_T, ice_age_T


# Quadratic fit of d18O - temperature data -----------------------------------------------------------------------------
def fit_d18O_to_T(d18O_T, temp):
    list_d18O = list(d18O_T)
    list_temp = list(temp)
    list_d18O, list_temp = zip(*sorted(zip(list_d18O, list_temp)))
    popt_T, pcov_T = curve_fit(cubic_func, list_d18O, list_temp)
    a, b, c = popt_T[0], popt_T[1], popt_T[2]
    return a, b, c, np.array(list_d18O), np.array(list_temp)


# Exponential fit of d18O - accumulation data --------------------------------------------------------------------------
def fit_d18O_to_Acc(d18O_T, acc):
    list_d18O = list(d18O_T)
    list_acc = list(acc)
    list_d18O, list_acc = zip(*sorted(zip(list_d18O, list_acc)))
    popt_acc, pcov_acc = curve_fit(exponential_func, np.array(list_d18O), np.array(list_acc))
    a_, b_, c_, d_ = popt_acc[0], popt_acc[1], popt_acc[2], popt_acc[3]
    return a_, b_, c_, d_, np.array(list_d18O), np.array(list_acc)    # this is a sorted array of d18O values



d18O_T, ice_age_T = interpolate_d18O_ice_age(depth_scales, ice_age_scales, ice_age, d18O_data, depth_temp)
a, b, c, sorted_d18O, sorted_temp = fit_d18O_to_T(d18O_T[:], temp[:])
a_, b_, c_, d_, sorted_d18O_acc, sorted_acc = fit_d18O_to_Acc(d18O_T, acc)

plot_T = True
plot_acc = True

if plot_T:
    # Plot data / fit
    fig, axs = plt.subplots(2)  # , sharex=False, sharey=False
    fig.set_figheight(8)
    fig.set_figwidth(15)

    axs[0].plot(sorted_d18O, sorted_temp, 'bo', markersize=1, label='data')
    axs[0].plot(sorted_d18O, cubic_func(sorted_d18O, a, b, c), 'r',
                label='cubic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K' % (a, b, c), markersize=1)
    axs[0].set_xlabel('$\delta^{18}$O [‰]')
    axs[0].set_ylabel('Temperature [°C]')
    axs[0].grid(linestyle='--', color='gray', lw='0.5')
    axs[0].legend()

    axs[1].plot(ice_age_T, temp, 'bo', markersize=1, label='data')
    axs[1].plot(ice_age_T, cubic_func(d18O_T, a, b, c), 'r', label='fitted data')
    axs[1].set_xlabel('Time [yr]')
    axs[1].set_ylabel('Temperature [°C]')
    axs[1].grid(linestyle='--', color='gray', lw='0.5')
    axs[1].legend()
    plt.show()


if plot_acc:
    # Plot data / fit
    fig, axs = plt.subplots(2)  # , sharex=False, sharey=False
    fig.set_figheight(8)
    fig.set_figwidth(15)

    axs[0].plot(sorted_d18O_acc * 1000, sorted_acc, 'bo', markersize=1, label='data')
    axs[0].plot(sorted_d18O_acc * 1000, exponential_func(sorted_d18O_acc, a_, b_, c_, d_), 'r',
                label='exponential fit: a=%5.3f , b=%5.3f , c=%5.3f , d=%5.3f ' % (a_, b_, c_, d_), markersize=1)
    axs[0].set_xlabel('$\delta^{18}$O [‰]')
    axs[0].set_ylabel('Accumulation [m/yr] ice equ.')
    axs[0].grid(linestyle='--', color='gray', lw='0.5')
    axs[0].legend()

    axs[1].plot(ice_age_T, acc, 'bo', markersize=1, label='data')
    axs[1].plot(ice_age_T, exponential_func(d18O_T, a_, b_, c_, d_), 'r', label='fitted data')
    axs[1].set_xlabel('Time [yr]')
    axs[1].set_ylabel('Accumulation [m/yr] ice equ.')
    axs[1].grid(linestyle='--', color='gray', lw='0.5')
    axs[1].legend()
    plt.show()


'''
p = smooth_parameter(cop1, ice_age)
sp = csaps.CubicSmoothingSpline(ice_age, d18O_data, smooth=p)
d18O_smooth = sp(ice_age)

plt.plot(ice_age, d18O_data, 'bo', markersize=1)
plt.plot(ice_age, d18O_smooth, 'b')
plt.show()
'''
