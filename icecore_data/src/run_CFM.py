import os
import glob
import json
from read_d18O import read_data_d18O, get_interval_data_noTimeGrid
from smoothing_splines import smooth_data
from read_temp_acc import *
from read_d15N import *

if __name__ == '__main__':
    start_year_ = -39000  # start input year for the actual run (second main run)
    end_year_ = -33000  # end input year for the actual run (second main run)

    data_path = '~/projects/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
    results_path = '~/projects/finalResults/secondSpin/'
    resultsFileName_Spin = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_SPIN2_2022-05-31_01.hdf5'
    resultsFileName = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_2022-05-31_01.hdf5'

    json_compare = 'FirnAir_NGRIP_compare.json'

    depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)

    depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                   ice_age_full,
                                                                                   start_year_,
                                                                                   end_year_)

    temp, temp_err = read_temp(data_path)
    temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]
    input_temperature = np.array([ice_age_interval, temp_interval])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature, delimiter=",")

    acc = read_acc(data_path)
    acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
    input_acc = np.array([ice_age_interval, acc_interval])
    np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP_compare.json -n')

    model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_SPIN2_compare_2022-05-31_01.hdf5'
    os.chdir('../icecore_data/src/')
    data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'

    d15N2_model, iceAge_model, gasAge_model, deltaAge = get_d15N_model(model_path, mode='cod', firnair=True,
                                                                       cop=1 / 200.)
    # d15N2_data, d15N2_err = get_d15N_data(data_path, iceAge_model, cop=1 / 200.)
    plt.plot(gasAge_model, d15N2_model, label='model')

    plt.xlabel('GICC05modelext ice age [yr]')
    plt.ylabel('$\delta^{15}$N$_2$ [‰]')
    plt.legend()
    plt.grid(':')
    plt.plot()
    plt.show()