from scipy.optimize import least_squares
from read_d18O import *
from read_d15N import *
import os

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters

start_year_ = -80000  # start input year
end_year_ = -35000  # end input year

ice_age_cod_max = 1500
ice_age_cod_min = 500

cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline

# ----------------------------------------------------------------------------------------------------------------------
# Data

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'
model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_65_30kyr_300m_5yr.hdf5'

# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------
depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
depth_interval, d18O_interval, ice_age_interval, d18O_smooth_, time_ = get_interval_data(depth_full, d18O_full,
                                                                                         ice_age_full,
                                                                                         start_year_, end_year_,
                                                                                         time_grid_stp_,
                                                                                         cop_)


# ----------------------------------------------------------------------------------------------------------------------
# Define the function to optimize
# ---------------------------------

def fun(theta):
    a = theta[0]
    b = theta[1]
    c = theta[2]
    temperature = a * d18O_smooth_ ** 2 + b * d18O_smooth_ + c
    input_temperature = np.array([time_, temperature])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')
    d15N2_model, iceAge_model, gasAge_model, deltaAge = get_d15N_model(model_path, mode='cod', cop=1 / 200.)
    d15N2_data = get_d15N_data(data_path, iceAge_model, cop=1 / 200.)
    '''
    plt.plot(iceAge_model, d15N2_model, label='Barnola')
    plt.plot(iceAge_model, d15N2_data, 'ko', label='data', markersize=1)
    plt.xlabel('GICC05modelext ice age [yr]')
    plt.ylabel('$\delta^{15}$N [‰]')
    plt.legend()
    plt.grid(':')
    plt.show()
    plt.savefig('inversion.pdf')
    plt.close()
    '''

    return abs(d15N2_model - d15N2_data)


def plot_fun(theta0, theta1):
    a0 = theta0[0]
    b0 = theta0[1]
    c0 = theta0[2]

    a1 = theta1[0]
    b1 = theta1[1]
    c1 = theta1[2]
    temperature0 = a0 * d18O_smooth_ ** 2 + b0 * d18O_smooth_ + c0
    input_temperature0 = np.array([time_, temperature0])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature0, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')
    d15N2_model0, iceAge_model0, gasAge_model0, deltaAge0 = get_d15N_model(model_path, mode='cod', cop=1 / 200.)

    temperature1 = a1 * d18O_smooth_ ** 2 + b1 * d18O_smooth_ + c1
    input_temperature1 = np.array([time_, temperature1])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature1, delimiter=",")
    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP.json -n')
    os.chdir('../icecore_data/src/')
    d15N2_model1, iceAge_model1, gasAge_model1, deltaAge1 = get_d15N_model(model_path, mode='cod', cop=1 / 200.)

    d15N2_data = get_d15N_data(data_path, iceAge_model0, cop=1 / 200.)

    plt.plot(iceAge_model0, d15N2_model0, label='Barnola firs guess')
    plt.plot(iceAge_model1, d15N2_model1, label='Barnola optimization')
    plt.plot(iceAge_model0, d15N2_data, 'ko', label='data', markersize=1)
    plt.xlabel('GICC05modelext ice age [yr]')
    plt.ylabel('$\delta^{15}$N [‰]')
    plt.legend()
    plt.grid(':')
    plt.show()
    plt.savefig('inversion.pdf')
    plt.close()


# ----------------------------------------------------------------------------------------------------------------------
# Do the optimization
# ----------------------

theta0 = [25000, 6000, 150]
theta1 = [26874.41537224, 5855.47550825, 150.02516584]
theta_lb = [20000, 3000, 50]
theta_ub = [40000, 9000, 250]

plot_fun(theta0, theta1)


# res = least_squares(fun, theta0, bounds=(theta_lb, theta_ub),  method='trf', x_scale=(1000., 1000., 100.), diff_step=(0.01, 0.01, 0.001))
# print(res.x)

