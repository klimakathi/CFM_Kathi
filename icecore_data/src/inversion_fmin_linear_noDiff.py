from scipy.optimize import fmin
from read_d18O import *
from read_d15N import *
from read_temp_acc import *
from calculate_d15N import *
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ----------------------------------------------------------------------------------------------------------------------
# Data

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
data_path2 = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'
model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_50_35kyr_300m_2yr_instant_acc.hdf5'

spin_path = 'resultsFolder/2022-05-24_01_resultsInversion_fmin_noDiffusion_SPIN.h5'
results_path = 'resultsFolder/2022-05-24_01_resultsInversion_fmin_noDiffusion.h5'

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters

start_year_ = -44000            # start input year
end_year_ = -38500              # end input year
year_Spin = 1000                # years of Second Spin --> this is to find a good start temperature for the optimization
end_year_Spin = start_year_ + year_Spin

firnair_module = False          # this is to specify whether we use the firnair module in the CFM

stpsPerYear = 0.5
S_PER_YEAR = 31557600.0

cop_ = 1 / 200.                 # cut-off frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20             # step length time grid --> also for cubic smoothing spline
cod_mode = 'cod'                # 'cod', 'lid', '0_diff'

optimizer = 'minimize'          # 'least_squares', 'minimize'
method = 'Nelder-Mead'          # 'BFGS', 'Nelder-Mead'
theta_0 = [0.42, 75]            # initial guess
N = 1000                        # number of max iterations

d15n_age = 'ice_age'            # 'gas_age', 'ice_age'  NOTE: Until now 'gas_age' only works if firnair_module=True !!!
frac_minimizer_interval = 0.5   # fraction of ice_age/d15N interval where optimization is performed
no_points_minimize_Spin = 5     # first points of ice_age/d15N interval where spin optimization is performed

# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------

depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
temp, temp_err = read_temp(data_path)
acc = read_acc(data_path)

# For Spin run ---------------------------------------------------------------------------------------------------------
depth_interval_Spin, d18O_interval_Spin, ice_age_interval_Spin = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                              ice_age_full,
                                                                                              start_year_,
                                                                                              end_year_Spin)
d18O_interval_perm_Spin = d18O_interval_Spin * 1000
d18o_smooth_Spin = smooth_data(1 / 200., d18O_interval_perm_Spin, ice_age_interval_Spin, ice_age_interval_Spin)[0]

temp_interval_Spin = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_Spin)[0]
acc_interval_Spin = get_interval_acc(acc, ice_age_full, start_year_, end_year_Spin)
input_acc_Spin = np.array([ice_age_interval_Spin, acc_interval_Spin])
np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc_Spin, delimiter=",")

years_Spin = (np.max(ice_age_interval_Spin) - np.min(ice_age_interval_Spin)) * 1.0
dt_Spin = S_PER_YEAR / stpsPerYear  # seconds per time step
stp = int(years_Spin * S_PER_YEAR / dt_Spin)  # -1       # total number of time steps, as integer
modeltime_Spin = np.linspace(start_year_, end_year_Spin, stp + 1)[:-1]

opt_dict_Spin = {'count_Spin': np.zeros([N, 1], dtype=int),
                 'a_Spin': np.zeros([N, 1]),
                 'b_Spin': np.zeros([N, 1]),
                 'c_Spin': np.zeros([N, 1]),
                 'd_Spin': np.zeros([N, 1]),
                 'd15N@cod_Spin': np.zeros([N, np.shape(modeltime_Spin)[0]]),
                 'd15N_data_Spin': np.zeros([N, np.shape(modeltime_Spin)[0]]),
                 'd15N_data_err_Spin': np.zeros([N, np.shape(modeltime_Spin)[0]]),
                 'ice_age_Spin': np.zeros([N, np.shape(modeltime_Spin)[0]]),
                 'gas_age_Spin': np.zeros([N, np.shape(modeltime_Spin)[0]]),
                 'cost_function_Spin': np.zeros([N, 1])}

# For Main run ---------------------------------------------------------------------------------------------------------
depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                               ice_age_full,
                                                                               start_year_, end_year_)
d18O_interval_perm = d18O_interval * 1000
d18o_smooth = smooth_data(1 / 200., d18O_interval_perm, ice_age_interval, ice_age_interval)[0]

temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]
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
            'd15N_data': np.zeros([N, np.shape(modeltime)[0]]),
            'd15N_data_err': np.zeros([N, np.shape(modeltime)[0]]),
            'ice_age': np.zeros([N, np.shape(modeltime)[0]]),
            'gas_age': np.zeros([N, np.shape(modeltime)[0]]),
            'cost_function': np.zeros([N, 1])}


# ----------------------------------------------------------------------------------------------------------------------
# Define the function to optimize
# ---------------------------------
def fun_Spin(theta):
    count = int(np.max(opt_dict_Spin['count_Spin']))
    print('iteration', count)

    a = theta[0]
    b = theta[1]
    temperature_Spin = 1. / a * d18o_smooth_Spin + b
    input_temperature_Spin = np.array([ice_age_interval_Spin, temperature_Spin])

    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature_Spin, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP_noDiff.json -n')
    os.chdir('../icecore_data/src/')

    d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_path, mode=cod_mode,
                                                                           firnair=firnair_module, cop=1 / 200.)

    ice_age_data_interv, gas_age_data_interv, d15N2_data_interv, d15N2_data_err_interv = \
        get_d15N_data_interval(data_path2, iceAge_model_)

    d15N2_model_interp, gasAge_model_interp = interpolate_d15Nmodel_2_d15Ndata(d15N2_model_, iceAge_model_,
                                                                               gasAge_model_, ice_age_data_interv)

    cost_func = 1 / (np.shape(d15N2_model_interp[:no_points_minimize_Spin])[0] - 1) * \
                np.sum((d15N2_model_interp[:no_points_minimize_Spin] -
                        d15N2_data_interv[:no_points_minimize_Spin]) ** 2)

    opt_dict_Spin['a_Spin'][count] = a
    opt_dict_Spin['b_Spin'][count] = b
    opt_dict_Spin['d15N@cod_Spin'][count, :np.shape(d15N2_model_interp)[0]] = d15N2_model_interp[:]
    opt_dict_Spin['d15N_data_Spin'][count, :np.shape(d15N2_data_interv)[0]] = d15N2_data_interv[:]
    opt_dict_Spin['d15N_data_err_Spin'][count, :np.shape(d15N2_data_err_interv)[0]] = d15N2_data_err_interv[:]
    opt_dict_Spin['ice_age_Spin'][count, :np.shape(d15N2_model_interp)[0]] = ice_age_data_interv[:]
    opt_dict_Spin['gas_age_Spin'][count, :np.shape(d15N2_model_interp)[0]] = gasAge_model_interp[:]
    opt_dict_Spin['cost_function_Spin'][count] = cost_func
    count += 1
    opt_dict_Spin['count_Spin'][count] = count

    return cost_func


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
        os.system('python3 main.py FirnAir_NGRIP_noDiff.json -n')
    else:
        os.system('python3 main.py FirnAir_NGRIP_noDiff.json')
    os.chdir('../icecore_data/src/')

    d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_path, mode=cod_mode,
                                                                           firnair=firnair_module, cop=1 / 200.)

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
    opt_dict['d15N@cod'][count, :np.shape(d15N2_model_interp)[0]] = d15N2_model_interp[:]
    opt_dict['d15N_data'][count, :np.shape(d15N2_data_interv)[0]] = d15N2_data_interv[:]
    opt_dict['d15N_data_err'][count, :np.shape(d15N2_data_err_interv)[0]] = d15N2_data_err_interv[:]
    opt_dict['ice_age'][count, :np.shape(d15N2_model_interp)[0]] = ice_age_data_interv[:]
    opt_dict['gas_age'][count, :np.shape(d15N2_model_interp)[0]] = gasAge_model_interp[:]
    opt_dict['cost_function'][count] = cost_func
    count += 1
    opt_dict['count'][count] = count

    return cost_func


# ----------------------------------------------------------------------------------------------------------------------
# MINIMIZE
# ----------------------

# Spin run optimization
res_Spin = fmin(fun_Spin, theta_0)
entry_0 = np.where(opt_dict_Spin['count_Spin'] == 0)[0]
opt_dict_Spin['count_Spin'] = np.delete(opt_dict_Spin['count_Spin'], entry_0[1:])
opt_dict_Spin['count_Spin'] = opt_dict_Spin['count_Spin'][:-1]
max_int = np.shape(opt_dict_Spin['count_Spin'])[0]
with h5py.File(spin_path, 'w') as f_Spin:
    for key in opt_dict_Spin:
        f_Spin[key] = opt_dict_Spin[key][:max_int]
f_Spin.close()
theta_Spin = res_Spin

print('----------------------------------------------')
print('|            INFO MINIMIZE SPIN              |')
print('----------------------------------------------')
print('Theta_Spin: a= %5.3f ‰/K, b = %5.3f K' % (theta_Spin[0], theta_Spin[1]))

# Main run optimization
res_Main = fmin(fun, theta_Spin)
entry_0 = np.where(opt_dict['count'] == 0)[0]
opt_dict['count'] = np.delete(opt_dict['count'], entry_0[1:])
opt_dict['count'] = opt_dict['count'][:-1]
max_int = np.shape(opt_dict['count'])[0]
with h5py.File(results_path, 'w') as f:
    for key in opt_dict:
        f[key] = opt_dict[key][:max_int]
f.close()
theta_Main = res_Main[0]

print('----------------------------------------------')
print('|            INFO MINIMIZE MAIN              |')
print('----------------------------------------------')
print('Theta_Main: a= %5.3f ‰/K, b = %5.3f K' % (theta_Main[0], theta_Main[1]))
