import numpy as np
import matplotlib.pyplot as plt
import h5py
from secondSpin import *
import os
from read_d18O import read_data_d18O, get_interval_data_noTimeGrid
from smoothing_splines import smooth_data
from read_temp_acc import get_interval_acc, read_acc


if __name__ == '__main__':
    start_year_ = -44000            # start input year for the actual run (second main run)
    end_year_ = -38500              # end input year for the actual run (second main run)
    year_Spin = 3000                # Years of first Spin (with constant temperature and accumulation)
    year_Spin2 = 2500               # Years of second Spin

    data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'

    depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
    depth_interval_Spin, d18O_interval_Spin, ice_age_interval_Spin = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                                  ice_age_full,
                                                                                                  start_year_,
                                                                                                  start_year_ + year_Spin2)
    d18O_interval_perm_Spin = d18O_interval_Spin * 1000
    d18o_smooth_Spin = smooth_data(1 / 200., d18O_interval_perm_Spin, ice_age_interval_Spin, ice_age_interval_Spin)[0]

    acc = read_acc(data_path)
    acc_interval_Spin = get_interval_acc(acc, ice_age_full, start_year_, start_year_ + year_Spin2)
    input_acc_Spin = np.array([ice_age_interval_Spin, acc_interval_Spin])
    np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc_Spin, delimiter=",")


    os.chdir('../../CFM_main/')
    os.system('python3 main.py FirnAir_NGRIP_Spin2.json -n')


