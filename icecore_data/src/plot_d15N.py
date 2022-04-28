import matplotlib.pyplot as plt
from read_d18O import *
from read_d15N import *
from read_temp_acc import *
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters

start_year_ = -50000  # start input year
end_year_ = -10000  # end input year

cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
cod_mode = '0_diff'

# ----------------------------------------------------------------------------------------------------------------------
# Data

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_65_30kyr_300m_5yr.hdf5'

# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------
depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
depth_interval, d18O_interval, ice_age_interval = get_interval_data(depth_full, d18O_full,
                                                                    ice_age_full,
                                                                    start_year_, end_year_,
                                                                    time_grid_stp_,
                                                                    cop_)[:3]
acc = read_acc(data_path)
acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
temp, temp_err = read_temp(data_path)
temp_interval, temp_err_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)
input_temp = np.array([ice_age_interval, temp_interval])
input_acc = np.array([ice_age_interval, acc_interval])
print(temp_interval)
plt.plot(np.flipud(ice_age_interval * (-1)), np.flipud(temp_interval))
plt.grid()
plt.show()

np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")
np.savetxt('../../CFM_main/CFMinput/optimize_temp.csv', input_temp, delimiter=",")

os.chdir('../../CFM_main/')
os.system('python3 main.py FirnAir_NGRIP.json -n')
os.chdir('../icecore_data/src/')

d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)
d15N2_data_, d15N2_data_err_ = get_d15N_data(data_path, iceAge_model_, cop=1 / 200.)

plt.plot(iceAge_model_, d15N2_model_)
plt.plot(iceAge_model_, d15N2_data_)
plt.show()




