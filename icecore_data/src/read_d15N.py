import pandas as pd
import scipy.interpolate as interpolate
from smoothing_splines import *
import h5py
import matplotlib.pyplot as plt


def get_d15N_data(path_data, ice_age_d15N_model, cop):
    # Data from ice core
    # ---------------------------------------------------------------------------------------------------
    df3 = pd.read_excel(path_data, sheet_name='Sheet6')
    depth_data1 = np.flipud(np.array(df3[df3.columns[2]]))
    d15n_data = np.flipud(np.array(df3[df3.columns[3]]))

    df2 = pd.read_excel(path_data, sheet_name='Sheet5')
    ice_age_data = np.flipud(np.array(df2[df2.columns[3]])) * (-1)
    depth_data2 = np.flipud(np.array(df2[df2.columns[0]]))

    # Interpolate ice_age from d15N2 data to the ice_age from d15N_model
    # ---------------------------------------------------------------------------------------------------
    Model2Data = interpolate.interp1d(ice_age_data, depth_data2, 'linear', fill_value='extrapolate')
    depth_regrid = Model2Data(ice_age_d15N_model)

    # get corresponding d15N values
    # ---------------------------------------------------------------------------------------------------
    ND = interpolate.interp1d(depth_data1, d15n_data, 'linear', fill_value='extrapolate')
    d15n_data_regrid = ND(depth_regrid)

    # Apply cubic smoothing spline to d15N data
    # ---------------------------------------------------------------------------------------------------
    d15n_smooth = smooth_data(cop, d15n_data_regrid, ice_age_d15N_model, ice_age_d15N_model)[0]

    return d15n_smooth


def get_d15N_model(path_model, mode, cop):
    # Data from model
    # ---------------------------------------------------------------------------------------------------
    file = h5py.File(path_model, 'r')

    if mode == 'cod':
        # Get data at Martinerie close-off depth
        close_off_depth = file["BCO"][:, 2]
        depth_model = file['depth']
        d15n_model = file["d15N2"][:] - 1.
        d15n_cod = np.ones_like(close_off_depth)
        gas_age_model = file["gas_age"][:]
        ice_age_model = file["age"][:]
        gas_age_cod = np.ones_like(close_off_depth)
        ice_age_cod = np.ones_like(close_off_depth)

    if mode == 'lid':
        # Get data at LID
        close_off_depth = file["BCO"][:, 6]
        depth_model = file['depth']
        d15n_model = file["d15N2"][:] - 1.
        d15n_cod = np.ones_like(close_off_depth)
        gas_age_model = file["gas_age"][:]
        ice_age_model = file["age"][:]
        gas_age_cod = np.ones_like(close_off_depth)
        ice_age_cod = np.ones_like(close_off_depth)

    else:
        pass

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

    ice_age = modeltime - ice_age_cod_smooth  # signs seem to be wrong way round because of negative modeltime

    delta_age = ice_age_cod_smooth - gas_age_cod_smooth
    gas_age = ice_age + delta_age

    d15n_cod = d15n_cod[1:] * 1000

    return d15n_cod, ice_age, gas_age, delta_age

model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_65_30kyr_300m_5yr.hdf5'
data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'

# d15N2_model, iceAge_model, gasAge_model, deltaAge = get_d15N_model(model_path, mode='cod', cop=1/200.)
# d15N2_data = get_d15N_data(data_path, iceAge_model, cop=1/200.)

# plt.plot(iceAge_model, d15N2_model * 1000)
# plt.plot(iceAge_model)
# plt.show()
# plt.plot(iceAge_model, d15N2_data, 'ko', markersize=1)
# plt.show()
