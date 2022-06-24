import numpy as np
from secondSpin import find_index_from_year
from read_d18O import *
from read_d15N import *
from read_temp_acc import *
import h5py
from matplotlib.ticker import MaxNLocator

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data_path = '~/projects/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'

path1 = '../../../finalResults/inversion/2022-06-06_01/2022-06-06_01_resultsInversion_minimizer.h5'
path2 = '../../../finalResults/inversion/2022-06-06_02/2022-06-06_02_resultsInversion_minimizer.h5'
path3 = '../../../finalResults/inversion/2022-06-06_03/2022-06-06_03_resultsInversion_minimizer.h5'
path4 = '../../../finalResults/inversion/2022-06-06_04/2022-06-06_04_resultsInversion_minimizer.h5'
path5 = '../../../finalResults/inversion/2022-06-06_05/2022-06-06_05_resultsInversion_minimizer.h5'
path_list = [path1, path2, path3, path4, path5]

start_year = -36000
end_year1 = -34000
end_year2 = -32000
end_year3 = -30000
end_year4 = -28000
end_year5 = -26000

end_year_list = [end_year1, end_year2, end_year3, end_year4, end_year5]


# plots      -----------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(6, sharex=True, sharey=False)
fig.set_figheight(20)
fig.set_figwidth(12)
N = len(path_list)
color_list = ['blue', 'orange', 'green', 'red', 'cyan']

for i in range(N):
    print(i)
    model_path_MAIN = path_list[i]

    # ------------------------------------------------------------------------------------------------------------------
    # Set parameters
    start_year_ = -36000  # start input year for the actual run (main run)
    end_year_ = end_year_list[i]  # end input year for the actual run (main run)

    stpsPerYear = 0.5
    S_PER_YEAR = 31557600.0

    cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
    time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
    cod_mode = 'cod'

    d15n_age = 'gas_age'  # 'gas_age', 'ice_age'

    # ----------------------------------------------------------------------------------------------------------------------
    # Read d18O data from NGRIP
    # ----------------------------
    depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
    depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                   ice_age_full,
                                                                                   start_year_, end_year_)
    d18O_interval_perm = d18O_interval * 1000
    d18o_smooth = smooth_data(1 / 200., d18O_interval_perm, ice_age_interval, ice_age_interval)[0]

    acc = read_acc(data_path)
    acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
    input_acc = np.array([ice_age_interval, acc_interval])
    np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

    years = (np.max(ice_age_interval) - np.min(ice_age_interval)) * 1.0
    dt = S_PER_YEAR / stpsPerYear  # seconds per time step
    stp = int(years * S_PER_YEAR / dt)  # -1       # total number of time steps, as integer
    modeltime = np.linspace(start_year_, end_year_, stp + 1)[:-1]

    f = h5py.File(model_path_MAIN, 'r')
    a = f['a'][:]
    b = f['b'][:]
    cost = f['cost_function'][:]
    d15n = f['d15N@cod'][:]
    d15n_data = f['d15N_data'][:]
    d15n_data_err = f['d15N_data_err'][:]
    ice_age = f['ice_age'][:]
    gas_age = f['gas_age'][:]

    # temperature ----------------------------------------------------------------------------------------------------------

    t0 = 1./a[0] * d18o_smooth + b[0]
    t1 = 1./a[-1] * d18o_smooth + b[-1]

    temp, temp_err = read_temp(data_path)
    temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]


    # Ice Age & Gas Age & d15N2 --------------------------------------------------------------------------------------------

    # MAIN

    entry0_best = np.where(ice_age[-1, :] == 0.)
    ice_age_best = np.delete(ice_age[-1, :], entry0_best)
    gas_age_best = np.delete(gas_age[-1, :], entry0_best)
    d15n_best = np.delete(d15n[-1, :], entry0_best)
    d15n_data_best = np.delete(d15n_data[-1, :], entry0_best)
    d15n_data_err_best = np.delete(d15n_data_err[-1, :], entry0_best)

    cost_list = []
    index = []
    print('length cost: ', np.shape(cost)[0])
    for j in range(np.shape(cost)[0] - 1):
        if len(cost_list) == 0:
            if cost[j + 1] < cost[j]:
                cost_list.append(cost[j + 1])
                index.append(j + 1)
        else:
            if cost[j + 1] < cost[j] and cost[j + 1] < cost_list[-1]:
                cost_list.append(cost[j + 1])
                index.append(j + 1)
    costs = np.array(cost_list)

    axs[0].plot(ice_age_interval, t1,  linestyle='-', linewidth=0.5, color=color_list[i],
                label='Main: Best fit: a = %5.3f ‰/K, b = %5.3f K' % (a[-1], b[-1]))
    axs[0].grid(linestyle='--', color='gray', lw='0.5')

    axs[0].legend(loc='lower right')
    duration = end_year_list[i] - start_year

    if d15n_age =='ice_age':
        axs[i + 1].plot(ice_age_best, d15n_data_best, 'k', linewidth=0.5, label='Kindler et al. 2014')
        axs[i + 1].plot(ice_age_best, d15n_data_best + d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
        axs[i + 1].plot(ice_age_best, d15n_data_best - d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
        axs[i + 1].fill_between(ice_age_best, d15n_data_best - d15n_data_err_best,
                                d15n_data_best + d15n_data_err_best, alpha=0.2, facecolor='k')
        axs[i + 1].plot(ice_age_best, d15n_best, color=color_list[i], linewidth=0.5, linestyle='-',
                        label='Best fit: %s yr' % duration)
        axs[i + 1].set_xlabel('GICC05modelext ice age [yr]')
    else:
        axs[i + 1].plot(gas_age_best, d15n_data_best, 'k', linewidth=0.5, label='Kindler et al. 2014')
        axs[i + 1].plot(gas_age_best, d15n_data_best + d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
        axs[i + 1].plot(gas_age_best, d15n_data_best - d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
        axs[i + 1].fill_between(gas_age_best, d15n_data_best - d15n_data_err_best,
                                d15n_data_best + d15n_data_err_best, alpha=0.2, facecolor='k')
        axs[i + 1].plot(gas_age_best, d15n_best, color=color_list[i], linewidth=0.5, linestyle='-',
                        label='Best fit: %s yr' % duration)

    axs[i + 1].set_xlim([-36500, -25500])
    axs[i + 1].set_ylim([0.25, 0.6])
    axs[i + 1].set_ylabel('$\delta^{15}$N$_2$ [‰]')
    axs[i + 1].grid(linestyle='--', color='gray', lw='0.5')
    axs[i + 1].legend(loc='lower right')
    # print(i)
axs[0].set_xlabel('Age [yr]')
axs[0].set_ylabel('Temperature [°C]')
axs[N].set_xlabel('GICC05modelext gas age [yr]')

plt.show()


# plot forcing and cost function
fig, axs = plt.subplots(2, sharex=False, sharey=False)
fig.set_figheight(7)
fig.set_figwidth(12)
N = len(path_list)
color_list = ['blue', 'orange', 'green', 'red', 'cyan']


for i in range(N):
    print(i)
    model_path_MAIN = path_list[i]

    # ------------------------------------------------------------------------------------------------------------------
    # Set parameters
    start_year_ = -36000  # start input year for the actual run (main run)
    end_year_ = end_year_list[i]  # end input year for the actual run (main run)

    stpsPerYear = 0.5
    S_PER_YEAR = 31557600.0

    cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
    time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
    cod_mode = 'cod'

    d15n_age = 'gas_age'  # 'gas_age', 'ice_age'

    # ----------------------------------------------------------------------------------------------------------------------
    # Read d18O data from NGRIP
    # ----------------------------
    depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)
    depth_interval, d18O_interval, ice_age_interval = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                   ice_age_full,
                                                                                   start_year_, end_year_)
    d18O_interval_perm = d18O_interval * 1000
    d18o_smooth = smooth_data(1 / 200., d18O_interval_perm, ice_age_interval, ice_age_interval)[0]

    acc = read_acc(data_path)
    acc_interval = get_interval_acc(acc, ice_age_full, start_year_, end_year_)
    input_acc = np.array([ice_age_interval, acc_interval])
    np.savetxt('../../CFM_main/CFMinput/optimize_acc.csv', input_acc, delimiter=",")

    years = (np.max(ice_age_interval) - np.min(ice_age_interval)) * 1.0
    dt = S_PER_YEAR / stpsPerYear  # seconds per time step
    stp = int(years * S_PER_YEAR / dt)  # -1       # total number of time steps, as integer
    modeltime = np.linspace(start_year_, end_year_, stp + 1)[:-1]

    f = h5py.File(model_path_MAIN, 'r')
    a = f['a'][:]
    b = f['b'][:]
    cost = f['cost_function'][:]
    d15n = f['d15N@cod'][:]
    d15n_data = f['d15N_data'][:]
    d15n_data_err = f['d15N_data_err'][:]
    ice_age = f['ice_age'][:]
    gas_age = f['gas_age'][:]

    # temperature ----------------------------------------------------------------------------------------------------------

    t0 = 1./a[0] * d18o_smooth + b[0]
    t1 = 1./a[-1] * d18o_smooth + b[-1]

    temp, temp_err = read_temp(data_path)
    temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]


    # Ice Age & Gas Age & d15N2 --------------------------------------------------------------------------------------------

    # MAIN

    entry0_best = np.where(ice_age[-1, :] == 0.)
    ice_age_best = np.delete(ice_age[-1, :], entry0_best)
    gas_age_best = np.delete(gas_age[-1, :], entry0_best)
    d15n_best = np.delete(d15n[-1, :], entry0_best)
    d15n_data_best = np.delete(d15n_data[-1, :], entry0_best)
    d15n_data_err_best = np.delete(d15n_data_err[-1, :], entry0_best)

    duration = end_year_list[i] - start_year

    cost_list = []
    index = []
    print('length cost: ', np.shape(cost)[0])
    for j in range(np.shape(cost)[0] - 1):
        if len(cost_list) == 0:
            if cost[j + 1] < cost[j]:
                cost_list.append(cost[j + 1])
                index.append(j + 1)
        else:
            if cost[j + 1] < cost[j] and cost[j + 1] < cost_list[-1]:
                cost_list.append(cost[j + 1])
                index.append(j + 1)
    costs = np.array(cost_list)

    axs[0].plot(ice_age_interval, t1,  linestyle='-', linewidth=0.5, color=color_list[i],
                label='Best fit: %s yr' % duration)
    axs[0].grid(linestyle='--', color='gray', lw='0.5')
    axs[0].set_xlabel('Age [yr]')
    axs[0].set_ylabel('Temperature [°C]')
    axs[0].legend(loc='lower right')

    duration = end_year_list[i] - start_year
    axs[1].plot(cost_list, linewidth=0.5, color=color_list[i], label='%s yr' % duration)
    axs[1].set_xlabel('Iteration')
    axs[1].set_ylabel('Cost function')
    # axs[2].set_xlim(0, 20)
    axs[1].grid(linestyle='--', color='gray', lw='0.5')
    axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    axs[1].legend()

axs[0].plot(ice_age_interval, temp_interval, 'k', linewidth=0.5, label='Kindler et al. 2014')
plt.show()


'''

'''