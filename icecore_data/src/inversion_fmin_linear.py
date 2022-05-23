from scipy.optimize import fmin
from read_d18O import *
from read_d15N import *
from read_temp_acc import *
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ----------------------------------------------------------------------------------------------------------------------
# Data

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
data_path2 = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'
model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_50_35kyr_300m_2yr_instant_acc.hdf5'
results_path = 'resultsFolder/2022-05-23_01_resultsInversion_fmin.h5'

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters

start_year_ = -44000  # start input year
end_year_ = -38500  # end input year
stpsPerYear = 0.5
S_PER_YEAR = 31557600.0

cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
cod_mode = 'cod'

optimizer = 'minimize'  # 'least_squares', 'minimize'
method = 'Nelder-Mead'  # 'BFGS', 'Nelder-Mead'
theta_0 = [0.42, 75]  # initial guess
N = 1000  # number of max iterations

d15n_age = 'ice_age'  # 'gas_age', 'ice_age'
frac_minimizer_interval = 0.5

# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------

depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                               ice_age_full,
                                                                               start_year_, end_year_,
                                                                               cop_)
d18O_interval_perm = d18O_interval * 1000
d18o_smooth = smooth_data(1 / 200., d18O_interval_perm, ice_age_interval, ice_age_interval)[0]

t = 1. / theta_0[0] * d18o_smooth + theta_0[1]

temp, temp_err = read_temp(data_path)
temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]
plt.plot(t)
plt.plot(temp_interval)
plt.show()

acc = read_acc(data_path)
acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
input_acc = np.array([ice_age_interval, acc_interval])
np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

years = (np.max(ice_age_interval) - np.min(ice_age_interval)) * 1.0
dt = S_PER_YEAR / stpsPerYear  # seconds per time step
stp = int(years * S_PER_YEAR / dt)  # -1       # total number of time steps, as integer
modeltime = np.linspace(start_year_, end_year_, stp + 1)[:-1]
minimizer_interval = int(np.shape(modeltime)[0] * frac_minimizer_interval)

opt_dict = {'count': np.zeros([N, 1], dtype=int),
            'a': np.zeros([N, 1]),
            'b': np.zeros([N, 1]),
            'c': np.zeros([N, 1]),
            'd': np.zeros([N, 1]),
            'd15N@cod': np.zeros([N, np.shape(modeltime)[0]]),
            'ice_age': np.zeros([N, np.shape(modeltime)[0]]),
            'gas_age': np.zeros([N, np.shape(modeltime)[0]]),
            'cost_function': np.zeros([N, 1])}


# ----------------------------------------------------------------------------------------------------------------------
# Define the function to optimize
# ---------------------------------

def fun(theta):
    count = int(np.max(opt_dict['count']))
    print('iteration', count)

    a = theta[0]
    b = theta[1]
    temperature = 1. / a * d18o_smooth + b
    input_temperature = np.array([ice_age_interval, temperature])

    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature, delimiter=",")
    os.chdir('../../CFM_main/')
    if count == 0:
        os.system('python3 main.py FirnAir_NGRIP.json -n')
    else:
        os.system('python3 main.py FirnAir_NGRIP.json')
    os.chdir('../icecore_data/src/')

    d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)
    ice_age_data_interv, gas_age_data_interv, d15N2_data_interv, d15N2_data_err_interv = \
        get_d15N_data_interval(data_path2, iceAge_model_)
    d15N2_model_interp, gasAge_model_interp = interpolate_d15Nmodel_2_d15Ndata(d15N2_model_, iceAge_model_,
                                                                               gasAge_model_, ice_age_data_interv)

    cost_func = 1 / (np.shape(d15N2_model_interp[-minimizer_interval:])[0] - 1) * np.sum((d15N2_model_interp
                                                                                          [-minimizer_interval:] -
                                                                                          d15N2_data_interv[
                                                                                          -minimizer_interval:]) ** 2)

    opt_dict['a'][count] = a
    opt_dict['b'][count] = b
    opt_dict['d15N@cod'][count, :] = d15N2_model_interp[:]
    opt_dict['ice_age'][count, :] = ice_age_data_interv[:]
    opt_dict['gas_age'][count, :] = gasAge_model_interp[:]
    opt_dict['cost_function'][count] = cost_func
    count += 1
    opt_dict['count'][count] = count

    return cost_func


def plot_fun(theta0, theta1):
    a0 = theta0[0]
    b0 = theta0[1]

    a1 = theta1[0]
    b1 = theta1[1]

    temperature0 = 1. / a0 * d18o_smooth + b0
    input_temperature0 = np.array([ice_age_interval, temperature0])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature0, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')
    d15N2_model0, iceAge_model0, gasAge_model0, deltaAge0 = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)

    temperature1 = 1. / a1 * d18o_smooth + b1
    input_temperature1 = np.array([ice_age_interval, temperature1])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature1, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')
    d15N2_model1, iceAge_model1, gasAge_model1, deltaAge1 = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)

    if d15n_age == 'ice_age':
        d15N2_data_, d15N2_err_ = get_d15N_data(data_path, iceAge_model1)
        plt.plot(iceAge_model0, d15N2_model0, label='Barnola first guess: a=%5.3f K, b=%5.3f K' % (
            theta0[0], theta0[1]), linewidth=0.9)
        plt.plot(iceAge_model1, d15N2_model1, label='Barnola best fit: a=%5.3f K, b=%5.3f K' % (
            theta1[0], theta1[1]), linewidth=0.9)
        plt.plot(iceAge_model1, d15N2_data_, 'ko', label='data', markersize=1)
        plt.plot(iceAge_model1, d15N2_data_ + d15N2_err_, 'k-.', linewidth=0.9, alpha=0.5)
        plt.plot(iceAge_model1, d15N2_data_ - d15N2_err_, 'k-.', linewidth=0.9, alpha=0.5)
        plt.fill_between(iceAge_model1, d15N2_data_ - d15N2_err_, d15N2_data_ + d15N2_err_, alpha=0.2, facecolor='k')
        plt.title('Least squares fit')
        plt.xlabel('GICC05modelext ice age [yr]')
    else:
        d15N2_data_, d15N2_err_ = get_d15N_data_gasage(data_path, gasAge_model1)
        plt.plot(gasAge_model0, d15N2_model0, label='Barnola first guess: a=%5.3f K, b=%5.3f K' % (
            theta0[0], theta0[1]), linewidth=0.9)
        plt.plot(gasAge_model1, d15N2_model1, label='Barnola best fit: a=%5.3f K, b=%5.3f K' % (
            theta1[0], theta1[1]), linewidth=0.9)
        plt.plot(gasAge_model1, d15N2_data_, 'ko', label='data', markersize=1)
        plt.plot(gasAge_model1, d15N2_data_ + d15N2_err_, 'k-.', linewidth=0.9, alpha=0.5)
        plt.plot(gasAge_model1, d15N2_data_ - d15N2_err_, 'k-.', linewidth=0.9, alpha=0.5)
        plt.fill_between(gasAge_model1, d15N2_data_ - d15N2_err_, d15N2_data_ + d15N2_err_, alpha=0.2, facecolor='k')
        plt.title('Least squares fit')
        plt.xlabel('GICC05modelext gas age [yr]')

    plt.ylabel('$\delta^{15}$N$_2$ [‰]')
    plt.legend()
    plt.grid(':')
    plt.show()
    return 0


# ----------------------------------------------------------------------------------------------------------------------
# MINIMIZE
# ----------------------

res_c = fmin(fun, theta_0)
entry_0 = np.where(opt_dict['count'] == 0)[0]
opt_dict['count'] = np.delete(opt_dict['count'], entry_0[1:])
opt_dict['count'] = opt_dict['count'][:-1]
max_int = np.shape(opt_dict['count'])[0]
with h5py.File(results_path, 'w') as f:
    for key in opt_dict:
        f[key] = opt_dict[key][:max_int]
f.close()

theta_c_1 = res_c[0]

print('----------------------------------------------')
print('|            INFO MINIMIZE                   |')
print('----------------------------------------------')
print(res_c.message)
print(res_c.success)
print('Theta1: ', theta_c_1)
