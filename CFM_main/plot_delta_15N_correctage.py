import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import csaps

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'
model_path = 'resultsFolder/CFMresults_NGRIP_HLdynamic_65_30kyr.hdf5'
cop1 = 1 / 200.  # cut-off period for d15N


def smooth_parameter(cop, age_):
    dx = np.mean(np.diff(age_))  # mean distance between two points
    lamda = (1 / (2 * cop * np.pi)) ** 4 / dx  # eg. 8 in Enting1987
    p = 1 / (1 + lamda)  # close to eq. 13 in Enting1987
    return p


def get_model_data_d15N(path_model, cod):
    # Data from model
    # ------------------------------------------------------------------------------------------------------
    file = h5py.File(path_model, 'r')

    if cod:
        # Get data at close-off density
        close_off_depth = file["BCO"][:, 2]
        depth_model = file['depth']
        d15n = file["d15N2"][:] - 1.
        d15n_cod = np.ones_like(close_off_depth)
        gas_age_model = file["gas_age"][:]
        ice_age_model = file["age"][:]
        gas_age_cod = np.ones_like(close_off_depth)
        ice_age_cod = np.ones_like(close_off_depth)

    else:
        close_off_depth = file["BCO"][:, 6]
        depth_model = file['depth']
        d15n = file["d15N2"][:] - 1.
        d15n_cod = np.ones_like(close_off_depth)
        gas_age_model = file["gas_age"][:]
        ice_age_model = file["age"][:]
        gas_age_cod = np.ones_like(close_off_depth)
        ice_age_cod = np.ones_like(close_off_depth)


    for i in range(depth_model.shape[0]):
        index = int(np.where(depth_model[i, 1:] == close_off_depth[i])[0])
        gas_age_cod[i] = gas_age_model[i, index]
        ice_age_cod[i] = ice_age_model[i, index]
        d15n_cod[i] = d15n[i, index]

    gas_age_cod = gas_age_cod[1:]
    ice_age_cod = ice_age_cod[1:]
    modeltime = depth_model[1:, 0]

    cop1 = 1. / 200
    p = smooth_parameter(cop1, modeltime)
    sp = csaps.CubicSmoothingSpline(modeltime, gas_age_cod, smooth=p)
    gas_age_cod_smooth = sp(modeltime)

    sp2 = csaps.CubicSmoothingSpline(modeltime, ice_age_cod, smooth=p)
    ice_age_cod_smooth = sp2(modeltime)

    ice_age = modeltime - ice_age_cod_smooth  # signs seem to be wrong way round because of negative modeltime
    delta_age = ice_age_cod_smooth - gas_age_cod_smooth
    gas_age = ice_age + delta_age

    d15n_cod = d15n_cod[1:]

    return gas_age, d15n_cod, delta_age, ice_age


def get_icecore_d15N_data(path_data, cop, gas_age, delta_age, ice_age):
    # Data from ice core
    # ---------------------------------------------------------------------------------------------------
    df3 = pd.read_excel(path_data, sheet_name='Sheet6')
    depth_data = np.flipud(np.array(df3[df3.columns[2]]))
    d15n_data = np.flipud(np.array(df3[df3.columns[3]]))

    df2 = pd.read_excel(path_data, sheet_name='Sheet5')
    gas_age_data = np.flipud(np.array(df2[df2.columns[4]])) * (-1)
    ice_age_data = np.flipud(np.array(df2[df2.columns[3]])) * (-1)
    depth = np.flipud(np.array(df2[df2.columns[0]]))

    # get depth from data corresponding to the modelled gas_age --> depth_regrid
    # ---------------------------------------------------------------------------------------------------
    MD = interpolate.interp1d(ice_age_data, depth, 'linear', fill_value='extrapolate')
    depth_regrid = MD(ice_age)

    # get corresponding d15N values
    # ---------------------------------------------------------------------------------------------------
    ND = interpolate.interp1d(depth_data, d15n_data, 'linear', fill_value='extrapolate')
    d15n_data_regrid = ND(depth_regrid)

    # Apply cubic smoothing spline to d15N data
    # ---------------------------------------------------------------------------------------------------
    p = smooth_parameter(cop, ice_age)
    sp = csaps.CubicSmoothingSpline(ice_age, d15n_data_regrid, smooth=p)
    d15n_smooth = sp(ice_age)
    # plt.plot(gas_age, d15n_data_regrid)
    # plt.plot(gas_age, d15n_smooth)
    # plt.show()

    return d15n_smooth


   # return 0


a = get_model_data_d15N(model_path, cod=False)
a2 = get_model_data_d15N(model_path, cod=True)
b = get_icecore_d15N_data(data_path, cop1, a[0], a[2], a[3])
f = h5py.File(model_path, 'r')


'''
fig, ax = plt.subplots()
ax.plot(t / 1000, temps, 'o', markersize=1, color='blue')
ax.plot(t_grid / 1000, temp_smooth, color='blue')
ax.set_xlabel("Age GICC05modelext [kyr]")
ax.set_ylabel("Temperature [K]", color="blue")
ax2 = ax.twinx()
ax2.plot(t / 1000, accs, 'o', markersize=1, color="orange")
ax2.plot(t_grid / 1000, accs_smooth, color="orange")
ax2.set_ylabel("Accumulation ice equivalent [m/yr]", color="orange")
plt.show()
'''

plt.plot(a[3], a[1] * 1000, label='model LID')
plt.plot(a2[3], a2[1] * 1000, label='model COD')
plt.plot(a[3], b, label='data')
plt.xlabel('GICC05modelext ice age [yr]')
plt.ylabel('$\delta^{15}$N [‰]')

#plt.plot(T[:, 0], (T[:, 2] - np.mean(T[:, 2]))/np.std(T[:, 2]), label='T')
plt.legend()
plt.show()







