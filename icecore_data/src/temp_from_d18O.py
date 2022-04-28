import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import csaps
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'


def smooth_parameter(cop, age_):
    dx = np.mean(np.diff(age_))  # mean distance between two points
    lamda = (1 / (2 * cop * np.pi)) ** 4 / dx  # eg. 8 in Enting1987
    p = 1 / (1 + lamda)  # close to eq. 13 in Enting1987
    return p


cop1 = 1 / 1000.

# Read the data --------------------------------------------------------------------------------------------------------
# d18O and corresponding ice_age
df = pd.read_excel(data_path, sheet_name='Sheet7')
depth_data = np.flipud(np.array(df[df.columns[1]])[500:])
d18O_data = np.flipud(np.array(df[df.columns[2]])[500:]) / 1000.
ice_age = np.flipud(np.array(df[df.columns[0]])[500:]) * (-1)
p = smooth_parameter(cop1, ice_age)
sp = csaps.CubicSmoothingSpline(ice_age, d18O_data, smooth=p)
d18O_smooth = sp(ice_age)

# temperature and accumulation with corresponding depth
df2 = pd.read_excel(data_path, sheet_name='Sheet4')
temp = np.flipud(np.array(df2[df2.columns[4]]))
acc = np.flipud(np.array(df2[df2.columns[3]]))
depth_temp = np.flipud(np.array(df2[df2.columns[0]]))

# ice_age and corresponding depth
df3 = pd.read_excel(data_path, sheet_name='Sheet5')
ice_age_scales = np.flipud(np.array(df3[df3.columns[3]])) * (-1)
depth_scales = np.flipud(np.array(df3[df3.columns[0]]))

def test_func():
    df = pd.read_excel(data_path, sheet_name='Sheet7')
    d18O_data = np.array(df[df.columns[2]])
    ice_age = np.array(df[df.columns[0]])
    p = smooth_parameter(1/1000., ice_age)
    print(p)
    sp = csaps.CubicSmoothingSpline(ice_age, d18O_data, smooth=p)
    d18O_smooth = sp(ice_age)
    plt.figure(0)
    plt.plot(ice_age, d18O_data)
    plt.plot(ice_age, d18O_smooth)
    plt.show()
    a = 39000/1000**2
    b = 5
    c = 100
    temp = a * d18O_smooth**2 + b * d18O_smooth + c
    df2 = pd.read_excel(data_path, sheet_name='Sheet4')
    temp_data = np.array(df2[df2.columns[4]])
    depth_temp = np.array(df2[df2.columns[0]])

    df3 = pd.read_excel(data_path, sheet_name='Sheet5')
    ice_age_scales = np.array(df3[df3.columns[3]])
    depth_scales = np.array(df3[df3.columns[0]])

    DIA = interpolate.interp1d(ice_age_scales, depth_scales, 'linear', fill_value='extrapolate')
    depth_interpolate = DIA(ice_age)

    # interpolate temperature
    TD = interpolate.interp1d(depth_temp, temp_data, 'linear', fill_value='extrapolate')
    T_interpolate = TD(depth_interpolate)



    plt.figure(1)
    plt.plot(ice_age, temp)
    plt.plot(ice_age, temp_data)
    plt.show()
    return 0


def cubic_func(x, a, b, c):
    return a * x ** 2 + b * x + c


def exponential_func(x, a_, b_, c_, d_):
    return a_ * np.exp(b_ * x ** 2 + c_ * x + d_)


def exponential_func2(x, a_, b_, c_):
    return a_ * np.exp(b_ * x + c_)

def exponential_func3(x, a_, b_, c_, d_, e_):
    return a_ * np.exp(b_ * x**3 + c_ * x**2 + d_ * x * e_)


# Get corresponding temperature values to d18O -------------------------------------------------------------------------
def interpolate_d18O_ice_age(depth_scales, ice_age_scales, ice_age, temp, depth_temp, acc):
    # interpolate depth - ice_age
    DIA = interpolate.interp1d(ice_age_scales, depth_scales, 'linear', fill_value='extrapolate')
    depth_interpolate = DIA(ice_age)

    # interpolate temperature
    TD = interpolate.interp1d(depth_temp, temp, 'linear', fill_value='extrapolate')
    T_interpolate = TD(depth_interpolate)

    # interpolate accumulation
    AD = interpolate.interp1d(depth_temp, acc, 'linear', fill_value='extrapolate')
    acc_interpolate = AD(depth_interpolate)

    return T_interpolate, acc_interpolate


# Quadratic fit of d18O - temperature data -----------------------------------------------------------------------------
def fit_d18O_to_T(d18O_data, T_interpolate):
    list_d18O = list(d18O_data)
    list_temp = list(T_interpolate)
    popt_T, pcov_T = curve_fit(cubic_func, list_d18O, list_temp)
    T_err = np.sqrt(np.diag(pcov_T))
    a, b, c = popt_T[0], popt_T[1], popt_T[2]
    return a, b, c, np.array(list_d18O), np.array(list_temp), T_err


# Exponential fit of d18O - accumulation data --------------------------------------------------------------------------
def fit_d18O_to_Acc(d18O_data, acc_interpolate):
    list_d18O = list(d18O_data)
    list_acc = list(acc_interpolate)
    list_d18O, list_acc = zip(*sorted(zip(list_d18O, list_acc)))
    popt_acc, pcov_acc = curve_fit(exponential_func, np.array(list_d18O), np.array(list_acc), maxfev=20000)
    acc_err = np.sqrt(np.diag(pcov_acc))
    a_, b_, c_, d_ = popt_acc[0], popt_acc[1], popt_acc[2], popt_acc[3]
    return a_, b_, c_, d_, np.array(list_d18O), np.array(list_acc), acc_err  # this is a sorted array of d18O values


def fit_d18O_to_Acc2(d18O_data, acc_interpolate):
    list_d18O = list(d18O_data)
    list_acc = list(acc_interpolate)
    list_d18O, list_acc = zip(*sorted(zip(list_d18O, list_acc)))
    popt_acc, pcov_acc = curve_fit(exponential_func2, np.array(list_d18O), np.array(list_acc), maxfev=20000)
    acc_err = np.sqrt(np.diag(pcov_acc))
    a_, b_, c_ = popt_acc[0], popt_acc[1], popt_acc[2]
    return a_, b_, c_, np.array(list_d18O), np.array(list_acc), acc_err  # this is a sorted array of d18O values


def fit_d18O_to_Acc3(d18O_data, acc_interpolate):
    list_d18O = list(d18O_data)
    list_acc = list(acc_interpolate)
    list_d18O, list_acc = zip(*sorted(zip(list_d18O, list_acc)))
    popt_acc, pcov_acc = curve_fit(exponential_func3, np.array(list_d18O), np.array(list_acc), maxfev=20000)
    acc_err = np.sqrt(np.diag(pcov_acc))
    a_, b_, c_, d_, e_ = popt_acc[0], popt_acc[1], popt_acc[2], popt_acc[3], popt_acc[4]
    return a_, b_, c_, d_, e_, np.array(list_d18O), np.array(list_acc), acc_err  # this is a sorted array of d18O values

'''
T_interpolate, acc_interpolate = interpolate_d18O_ice_age(depth_scales, ice_age_scales, ice_age, temp, depth_temp, acc)

a, b, c, sorted_d18O, sorted_temp, T_err = fit_d18O_to_T(d18O_smooth[:], T_interpolate[:])
a_, b_, c_, d_, sorted_d18O_acc, sorted_acc, acc_err = fit_d18O_to_Acc(d18O_smooth, acc_interpolate)
a_2, b_2, c_2, sorted_d18O_acc2, sorted_acc2, acc_err2 = fit_d18O_to_Acc2(d18O_smooth, acc_interpolate)
a_3, b_3, c_3, d_3, e_3, sorted_d18O_acc3, sorted_acc3, acc_err3 = fit_d18O_to_Acc3(d18O_smooth, acc_interpolate)


plot_T = False
plot_acc = False

if plot_T:
    # Plot data / fit
    fig, axs = plt.subplots(2)  # , sharex=False, sharey=False
    fig.set_figheight(8)
    fig.set_figwidth(15)

    axs[0].plot(sorted_d18O * 1000, sorted_temp, 'bo', markersize=1, label='data')
    axs[0].plot(sorted_d18O * 1000, cubic_func(sorted_d18O, a, b, c), 'r',
                label='cubic fit: a=(%5.3f $\pm$ %5.3f) K, b=(%5.3f $\pm$ %5.3f) K, c=(%5.3f $\pm$ %5.3f) K' % (
                    a, T_err[0], b, T_err[1], c, T_err[2]), markersize=1)
    axs[0].set_xlabel('$\delta^{18}$O [‰]')
    axs[0].set_ylabel('Temperature [°C]')
    axs[0].grid(linestyle='--', color='gray', lw='0.5')
    axs[0].legend()

    axs[1].plot(ice_age, T_interpolate, 'bo', markersize=1, label='data')
    axs[1].plot(ice_age, cubic_func(d18O_data, a, b, c), 'r', label='fitted data')
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
                label='exponential fit: a=(%5.3f $\pm$ %5.3f) m/yr , b=(%5.3f $\pm$ %5.3f) , c=(%5.3f $\pm$ %5.3f) , '
                      'd=(%5.3f $\pm$ %5.3f)' % (a_, acc_err[0], b_, acc_err[1], c_, acc_err[2], d_, acc_err[3]),
                markersize=1)
    axs[0].plot(sorted_d18O_acc2 * 1000, exponential_func2(sorted_d18O_acc2, a_2, b_2, c_2), 'orange', label='fitted data 2')
    axs[0].plot(sorted_d18O_acc3 * 1000, exponential_func3(sorted_d18O_acc3, a_3, b_3, c_3, d_3, e_3), 'green', label='fitted data 3')
    axs[0].set_xlabel('$\delta^{18}$O [‰]')
    axs[0].set_ylabel('Accumulation [m/yr] ice equ.')
    axs[0].grid(linestyle='--', color='gray', lw='0.5')
    axs[0].legend()

    axs[1].plot(ice_age, acc_interpolate, 'bo', markersize=1, label='data')
    axs[1].plot(ice_age, exponential_func(d18O_data, a_, b_, c_, d_), 'r', label='fitted data')
    axs[1].plot(ice_age, exponential_func2(d18O_data, a_2, b_2, c_2), 'orange', label='fitted data 2')
    axs[1].plot(ice_age, exponential_func3(d18O_data, a_3, b_3, c_3, d_3, e_3), 'orange', label='fitted data 3')
    axs[1].set_xlabel('Time [yr]')
    axs[1].set_ylabel('Accumulation [m/yr] ice equ.')
    axs[1].grid(linestyle='--', color='gray', lw='0.5')
    axs[1].legend()
    plt.show()


p = smooth_parameter(cop1, ice_age)
sp = csaps.CubicSmoothingSpline(ice_age, d18O_data, smooth=p)
d18O_smooth = sp(ice_age)

plt.plot(ice_age, d18O_data, 'bo', markersize=1)
plt.plot(ice_age, d18O_smooth, 'b')
plt.show()
'''
if __name__ == "__main__":
    test_func()