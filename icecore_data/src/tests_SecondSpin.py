from secondSpin import *
import os
import glob
import json
from read_d18O import read_data_d18O, get_interval_data_noTimeGrid
from smoothing_splines import smooth_data
from read_temp_acc import *


if __name__ == '__main__':

    start_year_ = -46700  # start input year for the actual run (second main run)
    end_year_ = -30000  # end input year for the actual run (second main run)
    year_Spin = 3000  # Years of first Spin (with constant temperature and accumulation)
    year_Spin2 = 6000  # Years of second Spin
    start_year_Spin2 = start_year_ - year_Spin2 / 2
    end_year_Spin2 = start_year_ + year_Spin2 / 2

    compare = True

    data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
    results_path = '~/projects/Thesis/finalResults/secondSpin/'
    resultsFileName_Spin = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_SPIN2_2022-05-30_01.hdf5'
    resultsFileName = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_2022-05-30_01.hdf5'

    json_SPIN = 'FirnAir_NGRIP.json'
    json_MAIN = 'FirnAir_NGRIP_Spin2.json'
    json_compare = 'FirnAir_NGRIP_compare.json'

    depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)

    # ------------------------------------------------------------------------------------------------------------------
    # SPIN -------------------------------------------------------------------------------------------------------------

    depth_interval_Spin, d18O_interval_Spin, ice_age_interval_Spin = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                                  ice_age_full,
                                                                                                  start_year_Spin2,
                                                                                                  end_year_)
    print(ice_age_interval_Spin[-1])
    d18O_interval_perm_Spin = d18O_interval_Spin * 1000
    d18o_smooth_Spin = smooth_data(1 / 200., d18O_interval_perm_Spin, ice_age_interval_Spin, ice_age_interval_Spin)[0]

    temp, temp_err = read_temp(data_path)
    temp_interval_Spin = get_interval_temp(temp, temp_err, ice_age_full, start_year_Spin2, end_year_)[0]
    input_temperature_Spin = np.array([ice_age_interval_Spin, temp_interval_Spin])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature_Spin, delimiter=",")

    acc = read_acc(data_path)
    acc_interval_Spin = get_interval_acc(acc, ice_age_full, start_year_Spin2, end_year_)
    input_acc_Spin = np.array([ice_age_interval_Spin, acc_interval_Spin])
    np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc_Spin, delimiter=",")

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

    os.system('python3 main.py FirnAir_NGRIP.json -n')

    model_path = glob.glob('resultsFolder/*.hdf5')[1]
    spin_path = glob.glob('resultsFolder/*.hdf5')[0]

    dict_spin = read_data_at_secondSpin(model_path, spin_path, start_year_)

    write_data_2_new_spinFile(spin_path, dict_spin)

    os.system('mv %s %s' % (model_path, results_path))

    # ------------------------------------------------------------------------------------------------------------------
    # MAIN -------------------------------------------------------------------------------------------------------------

    os.chdir('../icecore_data/src/')

    depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                   ice_age_full,
                                                                                   start_year_,
                                                                                   end_year_)
    print(ice_age_interval[0])
    d18O_interval_perm = d18O_interval * 1000
    d18o_smooth = smooth_data(1 / 200., d18O_interval_perm, ice_age_interval, ice_age_interval)[0]

    temp, temp_err = read_temp(data_path)
    temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]
    input_temperature = np.array([ice_age_interval, temp_interval])
    np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature, delimiter=",")

    acc = read_acc(data_path)
    acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
    input_acc = np.array([ice_age_interval, acc_interval])
    np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

    os.chdir('../../CFM_main/')

    with open('FirnAir_NGRIP_Spin2.json', 'r+') as json_file:
        cfm_params = json.load(json_file)
        cfm_params['TWriteStart'] = start_year_
        # cfm_params['SecondSpin'] = True
        cfm_params['resultsFileName'] = resultsFileName
        json_file.seek(0)
        json.dump(cfm_params, json_file, indent=4)

    os.system('python3 main.py FirnAir_NGRIP_Spin2.json')

    model_path = glob.glob('resultsFolder/*.hdf5')[1]
    spin_path = glob.glob('resultsFolder/*.hdf5')[0]
    os.system('mv %s %s' % (model_path, results_path))
    os.system('mv %s %s' % (spin_path, results_path))

    # ------------------------------------------------------------------------------------------------------------------
    # COMPARE ----------------------------------------------------------------------------------------------------------

    if compare:
        os.chdir('../icecore_data/src/')
        start_year_compare = -50000

        depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                       ice_age_full,
                                                                                       start_year_compare,
                                                                                       end_year_)
        d18O_interval_perm = d18O_interval * 1000
        d18o_smooth = smooth_data(1 / 200., d18O_interval_perm, ice_age_interval, ice_age_interval)[0]

        temp, temp_err = read_temp(data_path)
        temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_compare, end_year_)[0]
        input_temperature = np.array([ice_age_interval, temp_interval])
        np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature, delimiter=",")

        acc = read_acc(data_path)
        acc_interval = get_interval_acc(acc, ice_age_full, start_year_compare, end_year_)
        input_acc = np.array([ice_age_interval, acc_interval])
        np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

        os.chdir('../../CFM_main/')
        os.system('python3 main.py FirnAir_NGRIP_compare.json -n')
        model_path = glob.glob('resultsFolder/*.hdf5')[1]
        spin_path = glob.glob('resultsFolder/*.hdf5')[0]
        os.system('mv %s %s' % (model_path, results_path))
        os.system('mv %s %s' % (spin_path, results_path))

        os.chdir('../icecore_data/src/')

        depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                       ice_age_full,
                                                                                       start_year_,
                                                                                       end_year_)
        d18O_interval_perm = d18O_interval * 1000
        d18o_smooth = smooth_data(1 / 200., d18O_interval_perm, ice_age_interval, ice_age_interval)[0]

        temp, temp_err = read_temp(data_path)
        temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]
        input_temperature = np.array([ice_age_interval, temp_interval])
        np.savetxt('../../CFM_main/CFMinput/optimize_T.csv', input_temperature, delimiter=",")

        acc = read_acc(data_path)
        acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
        input_acc = np.array([ice_age_interval, acc_interval])
        np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

        os.chdir('../../CFM_main/')
        os.system('python3 main.py FirnAir_NGRIP_compare2.json -n')
        model_path = glob.glob('resultsFolder/*.hdf5')[1]
        spin_path = glob.glob('resultsFolder/*.hdf5')[0]
        os.system('mv %s %s' % (model_path, results_path))
        os.system('mv %s %s' % (spin_path, results_path))


