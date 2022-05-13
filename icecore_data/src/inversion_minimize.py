from scipy.optimize import minimize
from read_d18O import *
from read_d15N import *
from read_temp_acc import *
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters

start_year_ = -50000  # start input year
end_year_ = -45000  # end input year
stpsPerYear = 0.5
S_PER_YEAR = 31557600.0

cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
cod_mode = 'cod'

optimizer = 'minimize'                      # 'least_squares', 'minimize'
method = 'Nelder-Mead'                      # 'BFGS', 'Nelder-Mead'
theta_0 = [53013.175, 6612.959, 138.993]    # initial guess
N = 1000                                    # number of max iterations


# ----------------------------------------------------------------------------------------------------------------------
# Data

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_50_35kyr_300m_2yr_instant_acc.hdf5'


# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------

depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
depth_interval, d18O_interval, ice_age_interval = get_interval_data(depth_full, d18O_full,
                                                                    ice_age_full,
                                                                    start_year_, end_year_,
                                                                    time_grid_stp_,
                                                                    cop_)[:3]
d18o_smooth = smooth_data(1 / 200., d18O_interval, ice_age_interval, ice_age_interval)[0]

acc = read_acc(data_path)
acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
input_acc = np.array([ice_age_interval, acc_interval])
np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

years = (np.max(ice_age_interval) - np.min(ice_age_interval)) * 1.0
dt = S_PER_YEAR / stpsPerYear  # seconds per time step
stp = int(years * S_PER_YEAR / dt)  # -1       # total number of time steps, as integer
modeltime = np.linspace(start_year_, end_year_, stp + 1)[:-1]

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
    c = theta[2]
    temperature = a * d18o_smooth ** 2 + b * d18o_smooth + c
    input_temperature = np.array([ice_age_interval, temperature])

    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature, delimiter=",")
    os.chdir('../../CFM_main/')
    if count == 0:
        os.system('python3 main.py FirnAir_NGRIP.json -n')
    else:
        os.system('python3 main.py FirnAir_NGRIP.json')
    os.chdir('../icecore_data/src/')

    d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)
    d15N2_data_ = get_d15N_data(data_path, iceAge_model_, cop=1 / 200.)[0]
    cost_func = 1 / (np.shape(d15N2_model_)[0] - 1) * np.sum((d15N2_model_ - d15N2_data_) ** 2)

    opt_dict['a'][count] = a
    opt_dict['b'][count] = b
    opt_dict['c'][count] = c
    opt_dict['d15N@cod'][count, :] = d15N2_model_
    opt_dict['ice_age'][count, :] = iceAge_model_
    opt_dict['gas_age'][count, :] = gasAge_model_
    opt_dict['cost_function'][count] = cost_func
    count += 1
    opt_dict['count'][count] = count

    return cost_func


def plot_fun(theta0, theta1):
    a0 = theta0[0]
    b0 = theta0[1]
    c0 = theta0[2]

    a1 = theta1[0]
    b1 = theta1[1]
    c1 = theta1[2]
    temperature0 = a0 * d18o_smooth ** 2 + b0 * d18o_smooth + c0
    input_temperature0 = np.array([ice_age_interval, temperature0])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature0, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')
    d15N2_model0, iceAge_model0, gasAge_model0, deltaAge0 = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)

    temperature1 = a1 * d18o_smooth ** 2 + b1 * d18o_smooth + c1
    input_temperature1 = np.array([ice_age_interval, temperature1])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature1, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')
    d15N2_model1, iceAge_model1, gasAge_model1, deltaAge1 = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)

    d15N2_data_, d15N2_err_ = get_d15N_data(data_path, iceAge_model1, cop=1 / 200.)

    plt.plot(iceAge_model0, d15N2_model0, label='Barnola first guess: a=%5.3f K, b=%5.3f K, c=%5.3f K' % (
        theta0[0], theta0[1], theta0[2]), linewidth=0.9)
    plt.plot(iceAge_model1, d15N2_model1, label='Barnola best fit: a=%5.3f K, b=%5.3f K, c=%5.3f K' % (
        theta1[0], theta1[1], theta1[2]), linewidth=0.9)
    plt.plot(iceAge_model1, d15N2_data_, 'ko', label='data', markersize=1)
    plt.plot(iceAge_model1, d15N2_data_ + d15N2_err_, 'k-.', linewidth=0.9, alpha=0.5)
    plt.plot(iceAge_model1, d15N2_data_ - d15N2_err_, 'k-.', linewidth=0.9, alpha=0.5)
    plt.fill_between(iceAge_model1, d15N2_data_ - d15N2_err_, d15N2_data_ + d15N2_err_, alpha=0.2, facecolor='k')
    plt.title('Least squares fit')
    plt.xlabel('GICC05modelext ice age [yr]')
    plt.ylabel('$\delta^{15}$N$_2$ [‰]')
    plt.legend()
    plt.grid(':')
    plt.show()
    return 0


# ----------------------------------------------------------------------------------------------------------------------
# MINIMIZE
# ----------------------

if optimizer == 'minimize':
    res_c = minimize(fun, theta_0, method='Nelder-Mead', options={'maxiter': 70})
    entry_0 = np.where(opt_dict['count'] == 0)[0]
    opt_dict['count'] = np.delete(opt_dict['count'], entry_0[1:])
    opt_dict['count'] = opt_dict['count'][:-1]
    max_int = np.shape(opt_dict['count'])[0]
    with h5py.File('resultsFolder/resultsInversion_minimizer.h5', 'w') as f:
        for key in opt_dict:
            f[key] = opt_dict[key][:max_int]
    f.close()

    theta_c_1 = res_c.x

    print('----------------------------------------------')
    print('|            INFO MINIMIZE                   |')
    print('----------------------------------------------')
    print(res_c.message)
    print(res_c.success)
    print('Theta1: ', theta_c_1)

    plot_fun(theta_0, theta_c_1)
