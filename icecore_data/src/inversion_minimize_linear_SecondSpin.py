from scipy.optimize import minimize
from read_d18O import *
from read_d15N import *
from read_temp_acc import *
from calculate_d15N import *
from secondSpin import read_data_at_secondSpin, write_data_2_new_spinFile
import os
import json
import glob

# ----------------------------------------------------------------------------------------------------------------------
# Data & Model paths

data_path = '~/projects/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
data_path2 = '~/projects/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'

resultsFileName_Spin = 'CFMresults_NGRIP_Barnola_36-34kyr_300m_2yr_inversion-NM_SPIN2_2022-06-02_01.hdf5'
resultsFileName_Main = 'CFMresults_NGRIP_Barnola_36-34kyr_300m_2yr_inversion-NM_MAIN_2022-06-02_01.hdf5'

spin2_path = '../../CFM_main/resultsFolder/' + resultsFileName_Spin
model_path = '../../CFM_main/resultsFolder/' + resultsFileName_Main

finalResults_path_modelruns = '~/projects/finalResults/inversion/2022-06-01_01/'

json_SPIN = 'FirnAir_NGRIP.json'
json_MAIN = 'FirnAir_NGRIP_Spin2.json'

# optimization parameter files
results_minimizer_spin_path = 'resultsFolder/2022-06-02_01_resultsInversion_minimizer_SPIN.h5'
results_minimizer_main_path = 'resultsFolder/2022-06-02_01_resultsInversion_minimizer.h5'

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters

start_year_ = -36000  # start input year for the actual run (main run)
end_year_ = -34000  # end input year for the actual run (main run)
year_Spin = 3000  # Years of first Spin (with constant temperature and accumulation)
year_Spin2 = 8000  # Years of second Spin
overlap = 2000
start_year_Spin2 = start_year_ - (year_Spin2 - overlap)
end_year_Spin2 = start_year_ + overlap

firnair_module = True  # this is to specify whether we use the firnair module in the CFM

stpsPerYear = 0.5
S_PER_YEAR = 31557600.0

cop_ = 1 / 200.  # cut-off frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
cod_mode = 'cod'  # 'cod', 'lid', '0_diff'

optimizer = 'minimize'  # 'least_squares', 'minimize'
method = 'Nelder-Mead'  # 'BFGS', 'Nelder-Mead'
theta_0 = [0.37, 73]  # initial guess
N = 1000  # number of max iterations

d15n_age = 'ice_age'  # 'gas_age', 'ice_age'  NOTE: Until now 'gas_age' only works if firnair_module=True !!!
frac_minimizer_interval = 0.5  # fraction of ice_age/d15N interval where optimization is performed
no_points_minimize_Spin = 5  # first points of ice_age/d15N interval where spin optimization is performed

# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------

depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
temp, temp_err = read_temp(data_path)
acc = read_acc(data_path)

# For Spin run ---------------------------------------------------------------------------------------------------------
depth_interval_Spin, d18O_interval_Spin, ice_age_interval_Spin = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                              ice_age_full,
                                                                                              start_year_Spin2,
                                                                                              end_year_Spin2)
d18O_interval_perm_Spin = d18O_interval_Spin * 1000
d18o_smooth_Spin = smooth_data(1 / 200., d18O_interval_perm_Spin, ice_age_interval_Spin, ice_age_interval_Spin)[0]

t = 1. / theta_0[0] * d18o_smooth_Spin + theta_0[1]

temp, temp_err = read_temp(data_path)
temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_Spin2, end_year_Spin2)[0]
plt.plot(ice_age_interval_Spin, t)
plt.plot(ice_age_interval_Spin, temp_interval)
plt.show()


temp_interval_Spin = get_interval_temp(temp, temp_err, ice_age_full, start_year_Spin2, end_year_Spin2)[0]
acc_interval_Spin = get_interval_acc(acc, ice_age_full, start_year_Spin2, end_year_Spin2)
input_acc_Spin = np.array([ice_age_interval_Spin, acc_interval_Spin])
np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc_Spin, delimiter=",")

years_Spin = (np.max(ice_age_interval_Spin) - np.min(ice_age_interval_Spin)) * 1.0
dt_Spin = S_PER_YEAR / stpsPerYear  # seconds per time step
stp = int(years_Spin * S_PER_YEAR / dt_Spin)  # -1       # total number of time steps, as integer
modeltime_Spin = np.linspace(start_year_Spin2, end_year_Spin2, stp + 1)[:-1]

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
np.savetxt('../../CFM_main/CFMinput/optimize_acc_main.csv', input_acc, delimiter=",")

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
    print(ice_age_interval_Spin[0], ice_age_interval_Spin[-1])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature_Spin, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')

    d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(spin2_path, mode=cod_mode,
                                                                           firnair=firnair_module, cop=1 / 200.)

    ice_age_data_interv, gas_age_data_interv, d15N2_data_interv, d15N2_data_err_interv = \
        get_d15N_data_interval(data_path2, iceAge_model_)

    d15N2_model_interp, gasAge_model_interp = interpolate_d15Nmodel_2_d15Ndata(d15N2_model_, iceAge_model_,
                                                                               gasAge_model_, ice_age_data_interv)

    shape_optimize = int(np.shape(d15N2_model_interp)[0]/2)
    cost_func = 1 / (np.shape(d15N2_model_interp[-shape_optimize:])[0] - 1) * \
                np.sum((d15N2_model_interp[-shape_optimize:] -
                        d15N2_data_interv[-shape_optimize:]) ** 2)

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
    os.system('python3 main.py FirnAir_NGRIP_Spin2.json')
    os.chdir('../icecore_data/src/')

    d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_path, mode=cod_mode,
                                                                           firnair=firnair_module, cop=1 / 200.)

    ice_age_data_interv, gas_age_data_interv, d15N2_data_interv, d15N2_data_err_interv = \
        get_d15N_data_interval(data_path2, iceAge_model_)

    d15N2_model_interp, gasAge_model_interp = interpolate_d15Nmodel_2_d15Ndata(d15N2_model_, iceAge_model_,
                                                                               gasAge_model_, ice_age_data_interv)

    cost_func = 1 / (np.shape(d15N2_model_interp[:])[0] - 1) \
                * np.sum((d15N2_model_interp[:] - d15N2_data_interv[:]) ** 2)

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

# Spin prepare json
os.chdir('../../CFM_main/')
with open('FirnAir_NGRIP.json', 'r+') as json_file:
    cfm_params = json.load(json_file)
    cfm_params['TWriteStart'] = start_year_Spin2
    cfm_params['yearSpin'] = year_Spin
    # cfm_params['SecondSpin'] = False
    cfm_params['resultsFileName'] = resultsFileName_Spin
    json_file.seek(0)
    json.dump(cfm_params, json_file, indent=4)
    json_file.truncate()
    json_file.close()

# Spin run optimization
os.chdir('../icecore_data/src/')
res_Spin = minimize(fun_Spin, theta_0, method=method)
entry_0 = np.where(opt_dict_Spin['count_Spin'] == 0)[0]
opt_dict_Spin['count_Spin'] = np.delete(opt_dict_Spin['count_Spin'], entry_0[1:])
opt_dict_Spin['count_Spin'] = opt_dict_Spin['count_Spin'][:-1]
max_int = np.shape(opt_dict_Spin['count_Spin'])[0]
with h5py.File(results_minimizer_spin_path, 'w') as f_Spin:
    for key in opt_dict_Spin:
        f_Spin[key] = opt_dict_Spin[key][:max_int]
f_Spin.close()
theta_Spin = res_Spin.x

print('----------------------------------------------')
print('|            INFO MINIMIZE SPIN              |')
print('----------------------------------------------')
print(res_Spin.message)
print(res_Spin.success)
print('Theta Spin: ', theta_Spin)

os.chdir('../../CFM_main/')
model_path_2 = glob.glob('resultsFolder/CFMresults*.hdf5')[0]
spin_path_2 = glob.glob('resultsFolder/CFMspin*.hdf5')[0]

print('Reading Data at start year for main run ...')
dict_spin = read_data_at_secondSpin(model_path_2, spin_path_2, start_year_)
write_data_2_new_spinFile(spin_path_2, dict_spin)

os.system('mv %s %s' % (model_path_2, finalResults_path_modelruns))

# ----------------------------------------------------------------------------------------------------------------------
# Main run prepare json
with open('FirnAir_NGRIP_Spin2.json', 'r+') as json_file:
    cfm_params = json.load(json_file)
    cfm_params['TWriteStart'] = start_year_
    # cfm_params['SecondSpin'] = True
    cfm_params['resultsFileName'] = resultsFileName_Main
    json_file.seek(0)
    json.dump(cfm_params, json_file, indent=4)

# Main run optimization
os.chdir('../icecore_data/src/')
res_Main = minimize(fun, theta_Spin, method=method)
entry_0 = np.where(opt_dict['count'] == 0)[0]
opt_dict['count'] = np.delete(opt_dict['count'], entry_0[1:])
opt_dict['count'] = opt_dict['count'][:-1]
max_int = np.shape(opt_dict['count'])[0]
with h5py.File(results_minimizer_main_path, 'w') as f:
    for key in opt_dict:
        f[key] = opt_dict[key][:max_int]
f.close()
theta_Main = res_Main.x

print('----------------------------------------------')
print('|            INFO MINIMIZE MAIN              |')
print('----------------------------------------------')
print(res_Main.message)
print(res_Main.success)
print('Theta Main: ', theta_Main)

os.chdir('../../CFM_main/')
model_path_2 = glob.glob('resultsFolder/CFMresults*.hdf5')[0]
spin_path_2 = glob.glob('resultsFolder/CFMspin*.hdf5')[0]

os.system('mv %s %s' % (model_path_2, finalResults_path_modelruns))
os.system('mv %s %s' % (spin_path_2, finalResults_path_modelruns))
