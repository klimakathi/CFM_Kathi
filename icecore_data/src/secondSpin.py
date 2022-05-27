import h5py
import numpy as np


def find_index_from_year(time, year):
    year_ind = np.min(np.where(time >= year))
    return year_ind


def read_data_at_secondSpin(path_model, t_second_spin):
    f = h5py.File(path_model, 'r')
    time = f['depth'][1:, 0]
    t_second_spin_ind = find_index_from_year(time, t_second_spin)
    dict_SecondSpin = {
        'ageSpin2': f['age'][t_second_spin_ind, :],
        'densitySpin2': f['density'][t_second_spin_ind, :],
        'depthSpin2': f['depth'][t_second_spin_ind, :],
        'tempSpin2': f['temperature'][t_second_spin_ind, :],
        'diffusivitySpin2': f['diffusivity'][t_second_spin_ind, :],
        'gas_ageSpin2': f['gas_age'][t_second_spin_ind, :],
        'w_airSpin2': f['w_air'][t_second_spin_ind, :],
        'w_firnSpin2': f['w_firn'][t_second_spin_ind, :],
        'd15N2Spin2': f['d15N2'][t_second_spin_ind, :]
    }
    print(dict_SecondSpin['d15N2Spin2'])
    print(f['d15N2'][:])
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

    dict_spin = read_data_at_secondSpin(model_path, -39000)
    print(dict_spin)
    write_data_2_new_spinFile(spin_path, dict_spin)


