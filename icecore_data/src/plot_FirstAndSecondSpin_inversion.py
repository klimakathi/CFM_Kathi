import numpy as np
from secondSpin import find_index_from_year
from read_d18O import *
from read_d15N import *
from read_temp_acc import *
import h5py
from matplotlib.ticker import MaxNLocator

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
model_path_SPIN = '../../../finalResults/inversion/2022-05-31_01/2022-05-31_01_resultsInversion_minimizer_SPIN.h5'
model_path_MAIN = '../../../finalResults/inversion/2022-05-31_01/2022-05-31_01_resultsInversion_minimizer.h5'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters
start_year_ = -36000  # start input year for the actual run (main run)
end_year_ = -30000  # end input year for the actual run (main run)
year_Spin = 3000  # Years of first Spin (with constant temperature and accumulation)
year_Spin2 = 6000  # Years of second Spin
start_year_Spin2 = start_year_ - year_Spin2 / 2
end_year_Spin2 = start_year_ + year_Spin2 / 2

stpsPerYear = 0.5
S_PER_YEAR = 31557600.0

cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
cod_mode = 'cod'

d15n_age = 'gas_age'  # 'gas_age', 'ice_age'

# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------

# ----------------------------------------------------------------------------------------------------------------------
# SPIN
depth_full, d18O_full, ice_age_full = read_data_d18O(data_path)

depth_interval_Spin, d18O_interval_Spin, ice_age_interval_Spin = get_interval_data_noTimeGrid(depth_full, d18O_full,
                                                                                              ice_age_full,
                                                                                              start_year_Spin2,
                                                                                              end_year_Spin2)
d18O_interval_perm_Spin = d18O_interval_Spin * 1000
d18o_smooth_Spin = smooth_data(1 / 200., d18O_interval_perm_Spin, ice_age_interval_Spin, ice_age_interval_Spin)[0]

temp, temp_err = read_temp(data_path)
temp_interval_Spin = get_interval_temp(temp, temp_err, ice_age_full, start_year_Spin2, end_year_Spin2)[0]

years_Spin = (np.max(ice_age_interval_Spin) - np.min(ice_age_interval_Spin)) * 1.0
dt_Spin = S_PER_YEAR / stpsPerYear  # seconds per time step
stp = int(years_Spin * S_PER_YEAR / dt_Spin)  # -1       # total number of time steps, as integer
modeltime_Spin = np.linspace(start_year_Spin2, end_year_Spin2, stp + 1)[:-1]

f_Spin = h5py.File(model_path_SPIN, 'r')
a_Spin = f_Spin['a_Spin'][:]
b_Spin = f_Spin['b_Spin'][:]
cost_Spin = f_Spin['cost_function_Spin'][:]
d15n_Spin = f_Spin['d15N@cod_Spin'][:]
d15n_data_Spin = f_Spin['d15N_data_Spin'][:]
d15n_data_err_Spin = f_Spin['d15N_data_err_Spin'][:]
ice_age_Spin = f_Spin['ice_age_Spin'][:]
gas_age_Spin = f_Spin['gas_age_Spin'][:]


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
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
t0_Spin = 1./a_Spin[0] * d18o_smooth_Spin + b_Spin[0]
t1_Spin = 1./a_Spin[-1] * d18o_smooth_Spin + b_Spin[-1]

t0 = 1./a[0] * d18o_smooth + b[0]
t1 = 1./a[-1] * d18o_smooth + b[-1]

temp, temp_err = read_temp(data_path)
temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]

ice_age_Spin1 = np.linspace(ice_age_interval_Spin[0] - year_Spin, ice_age_interval_Spin[0] - 1, year_Spin)
t0_Spin1 = np.ones_like(ice_age_Spin1) * t0_Spin[0]
t1_Spin1 = np.ones_like(ice_age_Spin1) * t1_Spin[0]

# Ice Age & Gas Age & d15N2 --------------------------------------------------------------------------------------------
# SPIN
entry0_Spin_first = np.where(ice_age_Spin[0, :] == 0.)
ice_age_Spin_first = np.delete(ice_age_Spin[0, :], entry0_Spin_first)
gas_age_Spin_first = np.delete(gas_age_Spin[0, :], entry0_Spin_first)
d15n_Spin_first = np.delete(d15n_Spin[0, :], entry0_Spin_first)
d15n_data_Spin_first = np.delete(d15n_data_Spin[0, :], entry0_Spin_first)
d15n_data_err_Spin_first = np.delete(d15n_data_err_Spin[0, :], entry0_Spin_first)

entry0_Spin_best = np.where(ice_age_Spin[-1, :] == 0.)
ice_age_Spin_best = np.delete(ice_age_Spin[-1, :], entry0_Spin_best)
gas_age_Spin_best = np.delete(gas_age_Spin[-1, :], entry0_Spin_best)
d15n_Spin_best = np.delete(d15n_Spin[-1, :], entry0_Spin_best)
d15n_data_Spin_best = np.delete(d15n_data_Spin[-1, :], entry0_Spin_best)
d15n_data_err_Spin_best = np.delete(d15n_data_err_Spin[0, :], entry0_Spin_best)

# MAIN
entry0_first = np.where(ice_age[0, :] == 0.)
ice_age_first = np.delete(ice_age[0, :], entry0_first)
gas_age_first = np.delete(gas_age[0, :], entry0_first)
d15n_first = np.delete(d15n[0, :], entry0_first)
d15n_data_first = np.delete(d15n_data[0, :], entry0_first)

entry0_best = np.where(ice_age[-1, :] == 0.)
ice_age_best = np.delete(ice_age[-1, :], entry0_best)
gas_age_best = np.delete(gas_age[-1, :], entry0_best)
d15n_best = np.delete(d15n[-1, :], entry0_best)
d15n_data_best = np.delete(d15n_data[-1, :], entry0_best)
d15n_data_err_best = np.delete(d15n_data_err[-1, :], entry0_best)
shape_spin = int(np.shape(gas_age_Spin_best[:find_index_from_year(gas_age_Spin_best, gas_age_best[0])])[0])

# cost function --------------------------------------------------------------------------------------------------------
cost_list_Spin = []
index_Spin = []
print('length cost Spin: ', np.shape(cost_Spin)[0])
for i in range(np.shape(cost_Spin)[0] - 1):
    if len(cost_list_Spin) == 0:
        if cost_Spin[i + 1] < cost_Spin[i]:
            cost_list_Spin.append(cost_Spin[i + 1])
            index_Spin.append(i + 1)
    else:
        if cost_Spin[i + 1] < cost_Spin[i] and cost_Spin[i + 1] < cost_list_Spin[-1]:
            cost_list_Spin.append(cost_Spin[i + 1])
            index_Spin.append(i + 1)
costs_Spin = np.array(cost_list_Spin)

cost_list = []
index = []
print('length cost: ', np.shape(cost)[0])
for i in range(np.shape(cost)[0] - 1):
    if len(cost_list) == 0:
        if cost[i + 1] < cost[i]:
            cost_list.append(cost[i + 1])
            index.append(i + 1)
    else:
        if cost[i + 1] < cost[i] and cost[i + 1] < cost_list[-1]:
            cost_list.append(cost[i + 1])
            index.append(i + 1)
costs = np.array(cost_list)


# plots      -----------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3, sharex=False, sharey=False)
fig.set_figheight(10)
fig.set_figwidth(15)
# fig.suptitle('', fontsize=16)

axs[0].plot(ice_age_interval_Spin, temp_interval_Spin, 'k', linewidth=1)
axs[0].plot(ice_age_interval, temp_interval, 'k', linewidth=1, label='Kindler2014')
axs[0].plot(ice_age_interval_Spin, t0_Spin, linestyle='--', color='blue', linewidth=1,
            label='Second Spin: First guess: a = %5.3f ‰/K, b = %5.3f K' % (a_Spin[0], b_Spin[0]))
axs[0].plot(ice_age_Spin1, t0_Spin1, linestyle='--', color='blue', linewidth=1)
axs[0].plot(ice_age_interval_Spin, t1_Spin, linestyle='-', color='blue', linewidth=1,
            label='Second Spin: Best fit: a = %5.3f ‰/K, b = %5.3f K' % (a_Spin[-1], b_Spin[-1]))
axs[0].plot(ice_age_Spin1, t1_Spin1, linestyle='-', color='blue', linewidth=1)
axs[0].plot(ice_age_interval, t0, color='red', linestyle='--', linewidth=1,
            label='Main: First guess: a = %5.3f ‰/K, b = %5.3f K' % (a[0], b[0]))
axs[0].plot(ice_age_interval, t1, color='red', linestyle='-', linewidth=1,
            label='Main: Best fit: a = %5.3f ‰/K, b = %5.3f K' % (a[-1], b[-1]))
axs[0].grid(linestyle='--', color='gray', lw='0.5')
axs[0].set_xlim(start_year_Spin2 - year_Spin - 500, end_year_ + 500)
axs[0].set_xlabel('Age [yr]')
axs[0].set_ylabel('Temperature [°C]')
axs[0].legend(loc='lower left')

if d15n_age =='ice_age':
    axs[1].plot(ice_age_best, d15n_data_best, 'k', linewidth=1, label='Kindler2014')
    axs[1].plot(ice_age_best, d15n_data_best + d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
    axs[1].plot(ice_age_best, d15n_data_best - d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
    axs[1].fill_between(ice_age_best, d15n_data_best - d15n_data_err_best,
                        d15n_data_best + d15n_data_err_best, alpha=0.2, facecolor='k')
    axs[1].plot(ice_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin], 'k', linewidth=1)
    axs[1].plot(ice_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin] + d15n_data_err_Spin_best[:shape_spin], 'k', linewidth=0.1, alpha=0.5)
    axs[1].plot(ice_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin] - d15n_data_err_Spin_best[:shape_spin], 'k', linewidth=0.1, alpha=0.5)
    axs[1].fill_between(ice_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin] - d15n_data_err_Spin_best[:shape_spin],
                        d15n_data_Spin_best[:shape_spin] + d15n_data_err_Spin_best[:shape_spin], alpha=0.2, facecolor='k')
    axs[1].plot(ice_age_first, d15n_first, linewidth=1, linestyle='--', color='red', label='First guess')
    axs[1].plot(ice_age_best, d15n_best, linewidth=1, linestyle='-', color='red', label='Best fit')
    axs[1].plot(ice_age_Spin_first, d15n_Spin_first, linewidth=1, linestyle='--', color='blue', label='First guess')
    axs[1].plot(ice_age_Spin_best, d15n_Spin_best, linewidth=1, linestyle='-', color='blue', label='Best fit')
    axs[1].set_xlabel('GICC05modelext ice age [yr]')
else:
    axs[1].plot(gas_age_best, d15n_data_best, 'k', linewidth=1, label='Kindler2014')
    axs[1].plot(gas_age_best, d15n_data_best + d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
    axs[1].plot(gas_age_best, d15n_data_best - d15n_data_err_best, 'k', linewidth=0.1, alpha=0.5)
    axs[1].fill_between(gas_age_best, d15n_data_best - d15n_data_err_best,
                        d15n_data_best + d15n_data_err_best, alpha=0.2, facecolor='k')
    axs[1].plot(gas_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin], 'k', linewidth=1)
    axs[1].plot(gas_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin] + d15n_data_err_Spin_best[:shape_spin], 'k', linewidth=0.1, alpha=0.5)
    axs[1].plot(gas_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin] - d15n_data_err_Spin_best[:shape_spin], 'k', linewidth=0.1, alpha=0.5)
    axs[1].fill_between(gas_age_Spin_best[:shape_spin], d15n_data_Spin_best[:shape_spin] - d15n_data_err_Spin_best[:shape_spin],
                        d15n_data_Spin_best[:shape_spin] + d15n_data_err_Spin_best[:shape_spin], alpha=0.2, facecolor='k')
    axs[1].plot(gas_age_first, d15n_first, linewidth=1, linestyle='--', color='red', label='First guess')
    axs[1].plot(gas_age_best, d15n_best, linewidth=1, linestyle='-', color='red', label='Best fit')
    axs[1].plot(gas_age_Spin_first, d15n_Spin_first, linewidth=1, linestyle='--', color='blue', label='First guess')
    axs[1].plot(gas_age_Spin_best, d15n_Spin_best, linewidth=1, linestyle='-', color='blue', label='Best fit')
    axs[1].set_xlabel('GICC05modelext gas age [yr]')

axs[1].set_ylim([0.15, 0.55])
axs[1].set_xlim(start_year_Spin2 - year_Spin - 500, end_year_ + 500)
axs[1].set_ylabel('$\delta^{15}$N$_2$ [‰]')
axs[1].grid(linestyle='--', color='gray', lw='0.5')
axs[1].legend(loc='lower left')

axs[2].plot(costs_Spin, color='blue', linestyle='-', linewidth=1, label='Second spin')
axs[2].plot(costs, color='red', linestyle='-', linewidth=1, label='Main run')
axs[2].set_xlabel('Iteration')
axs[2].set_ylabel('Cost function')
# axs[2].set_xlim(0, 20)
axs[2].grid(linestyle='--', color='gray', lw='0.5')
axs[2].xaxis.set_major_locator(MaxNLocator(integer=True))
axs[2].legend()

plt.show()