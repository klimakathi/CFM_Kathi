from read_d18O import *
from read_d15N import *
from read_temp_acc import *
import h5py

data_path = '~/projects/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'
model_path = '../../../finalResults/inversion/2022-05-31_01/2022-05-31_01_resultsInversion_minimizer.h5'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ----------------------------------------------------------------------------------------------------------------------
# Set parameters

start_year_ = -44000  # start input year
end_year_ = -38500  # end input year
stpsPerYear = 0.5
S_PER_YEAR = 31557600.0

cop_ = 1 / 200.  # frequency for cubic smoothing spline (low pass filter)
time_grid_stp_ = 20  # step length time grid --> also for cubic smoothing spline
cod_mode = 'cod'

d15n_age = 'ice_age'  # 'gas_age'

# ----------------------------------------------------------------------------------------------------------------------
# Read d18O data from NGRIP
# -----------------------------

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


f = h5py.File(model_path, 'r')
a = f['a'][:]
b = f['b'][:]
cost = f['cost_function'][:]
d15n = f['d15N@cod'][:]
d15n_data = f['d15N_data'][:]
d15n_data_err = f['d15N_data_err'][:]
ice_age = f['ice_age'][:]
print(np.shape(ice_age[-1, :]))
gas_age = f['gas_age'][:]

# temperature ----------------------------------------------------------------------------------------------------------
t0 = 1./a[0] * d18o_smooth + b[0]
t1 = 1./a[-1] * d18o_smooth + b[-1]

temp, temp_err = read_temp(data_path)
temp_interval = get_interval_temp(temp, temp_err, ice_age_full, start_year_, end_year_)[0]

# cost function --------------------------------------------------------------------------------------------------------
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
print(costs[:10])


# plots      -----------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3, sharex=False, sharey=False)
fig.set_figheight(10)
fig.set_figwidth(15)
# fig.suptitle('', fontsize=16)

axs[0].plot(ice_age_interval, temp_interval, 'k', label='Kindler2014')
axs[0].plot(ice_age_interval, t0, label='First guess: a = %5.3f ‰/K, b = %5.3f K' % (a[0], b[0]))
axs[0].plot(ice_age_interval, t1, label='Best fit: a = %5.3f ‰/K, b = %5.3f K' % (a[-1], b[-1]))
axs[0].grid(linestyle='--', color='gray', lw='0.5')
axs[0].set_xlabel('Age [yr]')
axs[0].set_ylabel('Temperature [°C]')
axs[0].legend()

if d15n_age =='ice_age':
    axs[1].plot(ice_age[-1, :], d15n_data[-1, :], 'k', label='Kindler2014')
    axs[1].plot(ice_age[-1, :], d15n_data[-1, :], + d15n_data_err[-1, :], 'k', linewidth=0.1, alpha=0.5)
    axs[1].plot(ice_age[-1, :], d15n_data[-1, :], - d15n_data_err[-1, :], 'k', linewidth=0.1, alpha=0.5)
    axs[1].fill_between(ice_age[-1, :], d15n_data[-1, :] - d15n_data_err[-1, :],
                        d15n_data[-1, :] + d15n_data_err[-1, :], alpha=0.2, facecolor='k')
    axs[1].plot(ice_age[0, :], d15n[0, :], label='first guess')
    axs[1].plot(ice_age[-1, :], d15n[-1, :], label='Best fit')
    axs[1].set_xlabel('GICC05modelext ice age [yr]')
else:
    axs[1].plot(gas_age[-1, :], d15n_data[-1, :], 'k', label='Kindler2014')
    axs[1].plot(gas_age[-1, :], d15n_data[-1, :] + d15n_data_err[-1, :], 'k', linewidth=0.1, alpha=0.5)
    axs[1].plot(gas_age[-1, :], d15n_data[-1, :] - d15n_data_err[-1, :], 'k', linewidth=0.1, alpha=0.5)
    axs[1].fill_between(gas_age[-1, :], d15n_data[-1, :] - d15n_data_err[-1, :],
                        d15n_data[-1, :] + d15n_data_err[-1, :], alpha=0.2, facecolor='k')
    axs[1].plot(gas_age[0, :], d15n[0, :], label='first guess')
    axs[1].plot(gas_age[-1, :], d15n[-1, :], label='Best fit')
    axs[1].set_xlabel('GICC05modelext gas age [yr]')
axs[1].set_ylim([0.15, 0.55])
axs[1].set_ylabel('$\delta^{15}$N$_2$ [‰]')
axs[1].grid(linestyle='--', color='gray', lw='0.5')
axs[1].legend()

axs[2].plot(costs)
axs[2].set_xlabel('Iteration')
axs[2].set_ylabel('Cost function')
axs[2].grid(linestyle='--', color='gray', lw='0.5')

plt.show()