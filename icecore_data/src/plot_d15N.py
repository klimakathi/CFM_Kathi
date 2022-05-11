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
cod_mode = 'cod'


# ----------------------------------------------------------------------------------------------------------------------
# Data

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
path1 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/Barnola/CFMresults_NGRIP_Barnola_50_10kyr_300m_1yr.hdf5'
path2 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/Goujon/CFMresults_NGRIP_Goujon_50_10kyr_300m_1yr.hdf5'
path3 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/HLdynamic/CFMresults_NGRIP_HLdynamic_50_10kyr_300m_1yr.hdf5'
path4 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/HLSigfus/CFMresults_NGRIP_HLSigfus_50_10kyr_300m_1yr.hdf5'
path12 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/Barnola/CFMresults_NGRIP_Barnola_50_10kyr_300m_1yr_instant_acc.hdf5'
path22 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/Goujon/CFMresults_NGRIP_Goujon_50_10kyr_300m_1yr_instant_acc.hdf5'
path32 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/HLdynamic/CFMresults_NGRIP_HLdynamic_50_10kyr_300m_1yr_instant_acc.hdf5'
path42 = '../../CFM_main/results/Compare_d15Nmodel_to_d15Ndata_Kindler/HLSigfus/CFMresults_NGRIP_HLSigfus_50_10kyr_300m_1yr_instant_acc.hdf5'

model_paths = [path12, path22, path32, path42]

labels = ['Barnola', 'Goujon', 'HLdynamic', 'HLSigfus']

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
acc_smooth = smooth_data(cop_, acc_interval, ice_age_interval, ice_age_interval)

temp, temp_err = read_temp(data_path)
temp_interval, temp_err_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)
temp_smooth = smooth_data(cop_, temp_interval, ice_age_interval, ice_age_interval)

input_temp = np.array([ice_age_interval, temp_interval])
input_acc = np.array([ice_age_interval, acc_interval])

'''
np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")
np.savetxt('../../CFM_main/CFMinput/optimize_temp.csv', input_temp, delimiter=",")

os.chdir('../../CFM_main/')
os.system('python3 main.py FirnAir_NGRIP.json -n')
os.chdir('../icecore_data/src/')

model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_HLSigfus_50_10kyr_300m_1yr_instant_acc.hdf5'

d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_path, mode=cod_mode, cop=1 / 200.)
iceage = iceAge_model_
gasage = gasAge_model_
d15n = d15N2_model_

plt.plot(gasage, d15n)
d15N2_data_, d15N2_data_err_ = get_d15N_data(data_path, iceAge_model_, cop=1 / 200.)

plt.plot(gasage, d15N2_data_, 'k-', markersize=0.5, label='data Kindler')
plt.ylim([0.25, 0.6])
plt.show()
'''

plot = True
if plot:
    fig, axs = plt.subplots(2, sharex=False, sharey=False)
    fig.set_figheight(7)
    fig.set_figwidth(15)
    # fig.suptitle('', fontsize=16)

    axs[0].plot(ice_age_interval, temp_interval)
    axs[0].grid(linestyle='--', color='gray', lw='0.5')
    axs[0].set_xlabel('Age [yr]')
    axs[0].set_ylabel('Temperature [°C]')

    for i in range(len(model_paths)):
        d15N2_model_, iceAge_model_, gasAge_model_, deltaAge_ = get_d15N_model(model_paths[i], mode=cod_mode, cop=1/200.)
        iceage = iceAge_model_
        gasage = gasAge_model_
        d15n = d15N2_model_
        axs[1].plot(gasage, d15n, label=labels[i])

    d15N2_data_, d15N2_data_err_ = get_d15N_data_gasage(data_path, gasAge_model_, cop=1 / 200.)
    axs[1].plot(gasAge_model_, d15N2_data_, 'k-', markersize=0.5, label='data Kindler')
    axs[1].legend()
    axs[1].set_ylim([0.2, 0.6])
    axs[1].grid(linestyle='--', color='gray', lw='0.5')
    axs[1].set_xlabel('GICC05modelext gas age [yr]')
    axs[1].set_ylabel('$\delta^{15}$N$_2$ [‰]')
    plt.show()




