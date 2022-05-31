import matplotlib.pyplot as plt
import numpy as np
import h5py
from read_d15N import get_d15N_model


# Model Paths ----------------------------------------------------------------------------------------------------------
results_path = '../../../finalResults/secondSpin/first_try_secondSpin/2022-05-30_01/'

results_path_SPIN = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_SPIN2_2022-05-30_01.hdf5'
results_path_MAIN = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_2022-05-30_01.hdf5'
results_path_COMPARE = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_SPIN2_compare_2022-05-30_01.hdf5'
results_path_COMPARE2 = 'CFMresults_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_SPIN2_compare2_2022-05-30_01.hdf5'

spin_path_SPIN = 'CFMspin_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_2022-05-30_01.hdf5'
spin_path_COMPARE = 'CFMspin_NGRIP_Barnola_49_38kyr_300m_2yr_instant_acc_compare_2022-05-30_01.hdf5'


# Parameters -----------------------------------------------------------------------------------------------------------
cop = 1/200.
firnair = True
mode = 'cod'

start_year_ = -50000  # start input year for the actual run (second main run)
end_year_ = -38500  # end input year for the actual run (second main run)
year_Spin = 3000  # Years of first Spin (with constant temperature and accumulation)
year_Spin2 = 5000  # Years of second Spin
start_year_Spin2 = start_year_ - year_Spin2 / 2
end_year_Spin2 = start_year_ + year_Spin2 / 2

if __name__ == '__main__':
    # Read model outputs -----------------------------------------------------------------------------------------------
    f_spin2 = h5py.File(results_path + results_path_SPIN, 'r')
    f_main = h5py.File(results_path + results_path_MAIN, 'r')
    f_spin1 = h5py.File(results_path + spin_path_SPIN, 'r')
    f_compare = h5py.File(results_path + results_path_COMPARE, 'r')
    f_compare2 = h5py.File(results_path + results_path_COMPARE2, 'r')

    T_spin2 = f_spin2['Modelclimate'][:, 2]
    Acc_spin2 = f_spin2['Modelclimate'][:, 1]
    modeltime_spin2 = f_spin2['Modelclimate'][:, 0]

    T_spin1 = np.ones(year_Spin) * T_spin2[0]
    Acc_spin1 = np.ones(year_Spin) * Acc_spin2[0]
    modeltime_spin1 = np.arange(modeltime_spin2[0] - year_Spin, modeltime_spin2[0], 1)

    T_main = f_main['Modelclimate'][:, 2]
    Acc_main = f_main['Modelclimate'][:, 1]
    modeltime_main = f_main['Modelclimate'][:, 0]

    T_compare_main = f_compare['Modelclimate'][:, 2]
    Acc_compare_main = f_compare['Modelclimate'][:, 1]
    modeltime_compare_main = f_compare['Modelclimate'][:, 0]

    T_compare_spin = np.ones(year_Spin) * T_compare_main[0]
    Acc_compare_spin = np.ones(year_Spin) * Acc_compare_main[0]
    modeltime_compare_spin = np.arange(modeltime_compare_main[0] - year_Spin, modeltime_main[0], 1)

    d15n_spin2, ice_age_spin2, gas_age_spin2, delta_age_spin2 = \
        get_d15N_model(results_path + results_path_SPIN, mode, firnair, cop)

    d15n_main, ice_age_main, gas_age_main, delta_age_main = \
        get_d15N_model(results_path + results_path_MAIN, mode, firnair, cop)

    d15n_compare, ice_age_compare, gas_age_compare, delta_age_compare = \
        get_d15N_model(results_path + results_path_COMPARE, mode, firnair, cop)

    d15n_compare2, ice_age_compare2, gas_age_compare2, delta_age_compare2 = \
        get_d15N_model(results_path + results_path_COMPARE2, mode, firnair, cop)

    age_spin = f_spin1['ageSpin'][:]
    age_spin2 = f_spin1['ageSpin2'][:]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig, axs = plt.subplots(2, sharex=True, sharey=False)
    fig.set_figheight(7)
    fig.set_figwidth(15)
    # fig.suptitle('', fontsize=16)

    axs[0].plot(modeltime_spin1, T_spin1, label='Spin 1', color='blue')
    axs[0].plot(modeltime_spin2, T_spin2, label='Spin 2', color='orange')
    axs[0].plot(modeltime_main, T_main, label='Main', color='green')
    # axs[0].plot(modeltime_compare_spin, T_compare_spin, label='No second spin', color='red')
    # axs[0].plot(modeltime_compare_main, T_compare_main, color='red')
    axs[0].grid(linestyle='--', color='gray', lw='0.5')
    axs[0].set_xlabel('Age [yr]')
    axs[0].set_ylabel('Temperature [°C]')
    axs[0].legend()

    axs[1].plot(gas_age_spin2, d15n_spin2, label='Spin 2', color='orange')
    axs[1].plot(gas_age_main, d15n_main, label='Main', color='green')
    axs[1].plot(gas_age_compare2, d15n_compare2, label='No second spin', color='red')
    axs[1].plot(gas_age_compare, d15n_compare, 'k-', label='Long run')

    axs[1].legend()
    axs[1].grid(linestyle='--', color='gray', lw='0.5')
    axs[1].set_xlabel('GICC05modelext gas age [yr]')
    axs[1].set_ylabel('$\delta^{15}$N$_2$ [‰]')
    plt.show()

    # plt.plot(ice_age_spin2, d15n_spin2)
    # plt.plot(ice_age_main, d15n_main)
    # plt.plot(ice_age_compare, d15n_compare)
    # plt.show()

    # plt.plot(age_spin)
    # plt.plot(age_spin2)
    # plt.show()

    # plt.plot(modeltime_spin1, Acc_spin1)
    # plt.plot(modeltime_spin2, Acc_spin2)
    # plt.plot(modeltime_main, Acc_main)
    # plt.show()


