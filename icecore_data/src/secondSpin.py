import h5py
import numpy as np


def find_second_spin_index(ice_age_data, second_spin_year):  # second_spin_year is the end year of 2nd spin
    t_second_spin_ind = np.min(np.where(ice_age_data >= second_spin_year))
    return t_second_spin_ind


def read_data_at_secondSpinIndex(path_model, t_second_spin_ind):
    f = h5py.File(path_model, 'r')
    dict_SecondSpin = {
        'age_spin2': f['age'][t_second_spin_ind, :],
        'density_spin2': f['density'][t_second_spin_ind, :],
        'depth_spin2': f['depth'][t_second_spin_ind, :],
        'temp_spin2': f['temperature'][t_second_spin_ind, :],
        'diffusivity_spin2': f['diffusivity'][t_second_spin_ind, :],
        'gas_age_spin2': f['gas_age'][t_second_spin_ind, :],
        'w_air_spin2': f['w_air'][t_second_spin_ind, :],
        'w_firn_spin2': f['w_firn'][t_second_spin_ind, :],
        'd15n2_spin2': f['d15N2'][t_second_spin_ind, :]
    }
    f.close()
    return dict_SecondSpin


def write_data_2_new_spinFile(path_spin, dict_SecondSpin):
    with h5py.File(path_spin, 'w') as f:
        for key in dict_SecondSpin:
            f[key] = dict_SecondSpin[key][:]
    f.close()
    return


if __name__ == '__main__':

    model_path = '../../CFM_main/resultsFolder/CFMresults_NGRIP_Barnola_50_35kyr_300m_2yr_instant_acc.hdf5'
    spin_path = '../../CFM_main/resultsFolder/CFMspin_NGRIP_Barnola_50_35kyr_300m_2yr_instant_acc.hdf5'

    dict_spin = read_data_at_secondSpinIndex(model_path, 524)
    write_data_2_new_spinFile(spin_path, dict_spin)


