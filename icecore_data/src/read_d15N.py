import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from smoothing_splines import *
import h5py
import matplotlib.pyplot as plt


def get_d15N_data(path_data, ice_age_d15N_model):
    # Data from ice core
    # ---------------------------------------------------------------------------------------------------
    df = pd.read_excel(path_data)
    depth_data = np.flipud(np.array(df[df.columns[0]]))
    d15n_data = np.flipud(np.array(df[df.columns[9]]))
    d15n_err = np.flipud(np.array(df[df.columns[10]]))
    ice_age_data = np.flipud(np.array(df[df.columns[1]])) * (-1)

    # Interpolate ice_age from d15N2 data to the ice_age from d15N_model
    # ---------------------------------------------------------------------------------------------------
    Data2Model = interpolate.interp1d(ice_age_data, depth_data, 'linear', fill_value='extrapolate')
    depth_regrid = Data2Model(ice_age_d15N_model)

    # get corresponding d15N values
    # ---------------------------------------------------------------------------------------------------
    ND = interpolate.interp1d(depth_data, d15n_data, 'linear', fill_value='extrapolate')
    d15n_data_regrid = ND(depth_regrid)
    ND_err = interpolate.interp1d(depth_data, d15n_err, 'linear', fill_value='extrapolate')
    d15n_err_regrid = ND_err(depth_regrid)

    # Apply cubic smoothing spline to d15N data
    # ---------------------------------------------------------------------------------------------------
    # ice_age_d15N_model = remove_negative_values(ice_age_d15N_model)  # remove values with negative np.diff(ice_age)
    # d15n_smooth = smooth_data(cop, d15n_data_regrid, ice_age_d15N_model, ice_age_d15N_model)[0]
    # d15n_err_smooth = smooth_data(cop, d15n_err_regrid, ice_age_d15N_model, ice_age_d15N_model)[0]

    return d15n_data_regrid, d15n_err_regrid


def get_d15N_data_gasage(path_data, gas_age_d15N_model):
    # Data from ice core
    # ---------------------------------------------------------------------------------------------------
    df = pd.read_excel(path_data)
    depth_data = np.flipud(np.array(df[df.columns[0]]))
    d15n_data = np.flipud(np.array(df[df.columns[9]]))
    d15n_err = np.flipud(np.array(df[df.columns[10]]))
    gas_age_data = np.flipud(np.array(df[df.columns[2]])) * (-1)

    # Interpolate ice_age from d15N2 data to the gas_age from d15N_model
    # ---------------------------------------------------------------------------------------------------
    Data2Model = interpolate.interp1d(gas_age_data, depth_data, 'linear', fill_value='extrapolate')
    depth_regrid = Data2Model(gas_age_d15N_model)

    # get corresponding d15N values
    # ---------------------------------------------------------------------------------------------------
    ND = interpolate.interp1d(depth_data, d15n_data, 'linear', fill_value='extrapolate')
    d15n_data_regrid = ND(depth_regrid)
    ND_err = interpolate.interp1d(depth_data, d15n_err, 'linear', fill_value='extrapolate')
    d15n_err_regrid = ND_err(depth_regrid)

    # Apply cubic smoothing spline to d15N data
    # ---------------------------------------------------------------------------------------------------
    # ice_age_d15N_model = remove_negative_values(gas_age_d15N_model)  # remove values with negative np.diff(ice_age)
    # d15n_smooth = smooth_data(cop, d15n_data_regrid, gas_age_d15N_model, gas_age_d15N_model)[0]
    # d15n_err_smooth = smooth_data(cop, d15n_err_regrid, gas_age_d15N_model, gas_age_d15N_model)[0]

    return d15n_data_regrid, d15n_err_regrid


def get_d15N_model(path_model, mode, cop):
    # Data from model
    # ---------------------------------------------------------------------------------------------------
    file = h5py.File(path_model, 'r')

    if mode == 'cod':
        # Get data at Martinerie close-off depth
        close_off_depth = file["BCO"][:, 2]
        depth_model = file['depth']
        d15n_model = file["d15N2"][:] - 1.
        gas_age_model = file["gas_age"][:]
        ice_age_model = file["age"][:]
        d15n_cod = np.ones_like(close_off_depth)
        gas_age_cod = np.ones_like(close_off_depth)
        ice_age_cod = np.ones_like(close_off_depth)

    if mode == 'lid':
        # Get data at LID
        close_off_depth = file["BCO"][:, 6]
        depth_model = file['depth']
        d15n_model = file["d15N2"][:] - 1.
        gas_age_model = file["gas_age"][:]
        ice_age_model = file["age"][:]
        d15n_cod = np.ones_like(close_off_depth)
        gas_age_cod = np.ones_like(close_off_depth)
        ice_age_cod = np.ones_like(close_off_depth)

    if mode == '0_diff':
        # Get data at depth where D_eff = 0
        diffusivity = file['diffusivity'][:]
        depth_model = file['depth']
        d15n_model = file["d15N2"][:] - 1.
        gas_age_model = file["gas_age"][:]
        ice_age_model = file["age"][:]
        index = np.zeros(np.shape(diffusivity)[0])
        ice_age_cod = np.zeros(np.shape(diffusivity)[0])
        gas_age_cod = np.zeros(np.shape(diffusivity)[0])
        d15n_cod = np.zeros(np.shape(diffusivity)[0])
        close_off_depth = np.zeros(np.shape(diffusivity)[0])
        for i in range(np.shape(diffusivity)[0]):
            index[i] = np.max(np.where(diffusivity[i, 1:] > 10 ** (-20))) + 1
            ice_age_cod[i] = ice_age_model[i, int(index[i])]
            d15n_cod[i] = d15n_model[i, int(index[i])]
            close_off_depth[i] = depth_model[i, int(index[i])]

    for i in range(depth_model.shape[0]):
        index = int(np.where(depth_model[i, 1:] == close_off_depth[i])[0])
        gas_age_cod[i] = gas_age_model[i, index]
        ice_age_cod[i] = ice_age_model[i, index]
        d15n_cod[i] = d15n_model[i, index]

    gas_age_cod = gas_age_cod[1:]
    ice_age_cod = ice_age_cod[1:]
    modeltime = depth_model[1:, 0]

    gas_age_cod_smooth = smooth_data(cop, gas_age_cod, modeltime, modeltime)[0]
    ice_age_cod_smooth = smooth_data(cop, ice_age_cod, modeltime, modeltime)[0]

    ice_age = modeltime - ice_age_cod_smooth  # signs seem to be the wrong way round because of negative modeltime
    delta_age = ice_age_cod_smooth - gas_age_cod_smooth
    gas_age = ice_age + delta_age
    d15n_cod = d15n_cod[1:] * 1000

    return d15n_cod, ice_age, gas_age, delta_age


# ----------------------------------------------------------------------------------------------------------------------
# Test the function
# -------------------
test = False

if test:
    model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_65_30kyr_300m_5yr.hdf5'
    data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'

    d15N2_model, iceAge_model, gasAge_model, deltaAge = get_d15N_model(model_path, mode='cod', cop=1/200.)
    d15N2_data, d15N2_err = get_d15N_data(data_path, iceAge_model, cop=1/200.)

    plt.plot(iceAge_model, d15N2_model, label='model')
    plt.plot(iceAge_model, d15N2_data, 'ko', markersize=1, label='data')
    plt.plot(iceAge_model, d15N2_data + d15N2_err, 'k-.', linewidth=0.9, alpha=0.5)
    plt.plot(iceAge_model, d15N2_data - d15N2_err, 'k-.', linewidth=0.9, alpha=0.5)
    plt.fill_between(iceAge_model, d15N2_data - d15N2_err, d15N2_data + d15N2_err, alpha=0.2, facecolor='k')
    plt.xlabel('GICC05modelext ice age [yr]')
    plt.ylabel('$\delta^{15}$N$_2$ [‰]')
    plt.legend()
    plt.grid(':')
    plt.plot()
    plt.show()



