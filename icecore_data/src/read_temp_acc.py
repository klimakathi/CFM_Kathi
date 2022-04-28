import pandas as pd
import numpy as np
from read_d18O import find_start_end_index


def read_temp(path_data):
    df = pd.read_excel(path_data)
    temp = np.flipud(np.array(df[df.columns[3]]))
    temp_err = np.flipud(np.array(df[df.columns[4]]))
    return temp, temp_err


def read_acc(path_data):
    df = pd.read_excel(path_data)
    acc = np.flipud(np.array(df[df.columns[5]]))
    return acc


def get_interval_temp(temp, temp_err, ice_age_data, start_year, end_year):
    t_start_ind, t_end_ind = find_start_end_index(ice_age_data, start_year, end_year)
    temp_interval = temp[t_start_ind: t_end_ind]
    temp_err_interval = temp_err[t_start_ind: t_end_ind]
    return temp_interval, temp_err_interval


def get_interval_acc(acc, ice_age_data, start_year, end_year):
    t_start_ind, t_end_ind = find_start_end_index(ice_age_data, start_year, end_year)
    acc_interval = acc[t_start_ind: t_end_ind]
    return acc_interval


