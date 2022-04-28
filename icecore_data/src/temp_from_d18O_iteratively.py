import numpy as np

from temp_from_d18O import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

start_year_ = -80000
end_year_ = -35000

interval_ = 2000  # in years
step_width_ = 100  # in years


def find_start_end_index(ice_ages, start_year, end_year):
    t_start_ind = np.min(np.where(ice_ages >= start_year))
    t_end_ind = np.min(np.where(ice_ages >= end_year))
    return t_start_ind, t_end_ind


def calc_interval_step(step_width, ice_ages):
    interval_stp = int(step_width / np.diff(ice_ages)[0])
    return interval_stp


def calc_interval_int(interval, ice_ages):
    interval_int = int(interval / np.diff(ice_ages)[0])
    return interval_int


def number_iterations(end_year, start_year, interval, step_width):
    N = abs(int((abs(end_year - start_year) - interval) / step_width)) # number of iterations
    return N


def one_curve_fit(t_start_ind, t_end_ind, d18O_data, t_interpolate, ice_age):
    d18O_interval = d18O_data[t_start_ind: t_end_ind]
    T_interval = t_interpolate[t_start_ind: t_end_ind]
    ice_age_interval = ice_age[t_start_ind: t_end_ind]
    a_t, b_t, c_t, sorted_d18O_t, sorted_t, t_err = fit_d18O_to_T(d18O_interval, T_interval)
    temp_fit = cubic_func(d18O_interval, a_t, b_t, c_t)
    plt.plot(ice_age_interval, T_interval, 'bo', markersize=1)
    plt.plot(ice_age_interval, temp_fit)
    plt.show()
    return 0


def iterate_curve_fits(t_start_ind, interval_stp, interval_int, N, d18O_data, T_interpolate):
    curve_fit_temp_a = np.zeros([N, np.shape(d18O_data[t_start_ind: t_start_ind + (N-1) * interval_stp + interval_int])[0]])
    curve_fit_temp_b = np.zeros_like(curve_fit_temp_a)
    curve_fit_temp_c = np.zeros_like(curve_fit_temp_a)
    curve_fit_acc_a = np.zeros_like(curve_fit_temp_a)
    curve_fit_acc_b = np.zeros_like(curve_fit_temp_a)
    curve_fit_acc_c = np.zeros_like(curve_fit_temp_a)
    curve_fit_acc_d = np.zeros_like(curve_fit_temp_a)
    for i in range(N):
        ice_age_interval = ice_age[t_start_ind + i * interval_stp: t_start_ind + i * interval_stp + interval_int]
        d18O_interval = d18O_data[t_start_ind + i * interval_stp: t_start_ind + i * interval_stp + interval_int]
        # print(t_start_ind + i*interval_stp)
        T_interval = T_interpolate[t_start_ind + i * interval_stp: t_start_ind + i * interval_stp + interval_int]
        acc_interval = acc_interpolate[t_start_ind + i * interval_stp: t_start_ind + i * interval_stp + interval_int]

        a_t, b_t, c_t, sorted_d18O_t, sorted_t, t_err = fit_d18O_to_T(d18O_interval, T_interval)
        curve_fit_temp_a[i, i * interval_stp: i * interval_stp + interval_int] += a_t
        curve_fit_temp_b[i, i * interval_stp: i * interval_stp + interval_int] += b_t
        curve_fit_temp_c[i, i * interval_stp: i * interval_stp + interval_int] += c_t

        a_a, b_a, c_a, d_a, sorted_d18O_a, sorted_a, a_err = fit_d18O_to_Acc(d18O_interval, acc_interval)
        curve_fit_acc_a[i, i * interval_stp: i * interval_stp + interval_int] += a_a
        curve_fit_acc_b[i, i * interval_stp: i * interval_stp + interval_int] += b_a
        curve_fit_acc_c[i, i * interval_stp: i * interval_stp + interval_int] += c_a
        curve_fit_acc_d[i, i * interval_stp: i * interval_stp + interval_int] += d_a

    return curve_fit_temp_a, curve_fit_temp_b, curve_fit_temp_c, curve_fit_acc_a, curve_fit_acc_b, curve_fit_acc_c, curve_fit_acc_d


def best_fit(eval_type, d18O_data, ice_age, curve_fit_temp_a, curve_fit_temp_b, curve_fit_temp_c, curve_fit_acc_a, curve_fit_acc_b, curve_fit_acc_c, curve_fit_acc_d):
    if eval_type == 'mean':
        a_mean = np.true_divide(curve_fit_temp_a.sum(0), (curve_fit_temp_a != 0.).sum(0))
        b_mean = np.true_divide(curve_fit_temp_b.sum(0), (curve_fit_temp_b != 0.).sum(0))
        c_mean = np.true_divide(curve_fit_temp_c.sum(0), (curve_fit_temp_c != 0.).sum(0))
        a_acc_mean = np.true_divide(curve_fit_acc_a.sum(0), (curve_fit_acc_a != 0.).sum(0))
        b_acc_mean = np.true_divide(curve_fit_acc_b.sum(0), (curve_fit_acc_b != 0.).sum(0))
        c_acc_mean = np.true_divide(curve_fit_acc_c.sum(0), (curve_fit_acc_c != 0.).sum(0))
        d_acc_mean = np.true_divide(curve_fit_acc_d.sum(0), (curve_fit_acc_d != 0.).sum(0))
        # print('x: ', np.shape(d18O_data))
        # print('a: ', np.shape(a_mean))
    else:
        pass
    return a_mean, b_mean, c_mean, a_acc_mean, b_acc_mean, c_acc_mean, d_acc_mean


def plot_best_fit(d18O_data, ice_age, t_interpolate, acc_interpolate, a_mean, b_mean, c_mean, a_acc_mean, b_acc_mean, c_acc_mean, d_acc_mean):
    a_t, b_t, c_t, sorted_d18O_t, sorted_t, t_err = fit_d18O_to_T(d18O_data, t_interpolate)
    temp_fit = cubic_func(d18O_data, a_t, b_t, c_t)
    a_a, b_a, c_a, d_a, sorted_d18O_t, sorted_acc, acc_err = fit_d18O_to_Acc(d18O_data, acc_interpolate)
    acc_fit = exponential_func(d18O_data, a_a, b_a, c_a, d_a)

    temp_best_fit = cubic_func(d18O_data[:-5], a_mean, b_mean, c_mean)
    acc_best_fit = exponential_func(d18O_data[:-5], a_acc_mean, b_acc_mean, c_acc_mean, d_acc_mean)

    fig, axs = plt.subplots(2)
    fig.set_figheight(8)
    fig.set_figwidth(15)
    axs[0].plot(ice_age, temp_fit, label='fit', color='red')
    axs[0].plot(ice_age[:-5], temp_best_fit, label='best fit: mean', color='orange')
    axs[0].plot(ice_age, t_interpolate, 'ko', markersize=1, label='data Kindler')
    axs[0].set_xlabel('ice age [yr]')
    axs[0].set_ylabel('temperature [°C]')
    axs[0].legend()

    axs[1].plot(ice_age, acc_fit, label='fit', color='red')
    # axs[1].plot(ice_age[:-5], acc_best_fit, label='best fit: mean', color='orange')
    axs[1].plot(ice_age, acc_interpolate, 'ko', markersize=1, label='data Kindler')
    axs[1].set_xlabel('ice age [yr]')
    axs[1].set_ylabel('accumulation [m/yr]')
    axs[1].legend()
    plt.show()
    return 0





t_start_ind_, t_end_ind_ = find_start_end_index(ice_age, start_year_, end_year_)
print(t_start_ind_, t_end_ind_)
interval_stp_ = calc_interval_step(step_width_, ice_age)
print(interval_stp_)
interval_int_ = calc_interval_int(interval_, ice_age)
print(interval_int_)
n = number_iterations(end_year_, start_year_, interval_, step_width_)
print(n)
A, B, C, A_, B_, C_, D_ = iterate_curve_fits(t_start_ind_, interval_stp_, interval_int_, n, d18O_data, T_interpolate)

A_mean, B_mean, C_mean, A_mean_, B_mean_, C_mean_, D_mean_ = best_fit('mean', d18O_data[t_start_ind_: t_end_ind_], ice_age[t_start_ind_: t_end_ind_], A, B, C, A_, B_, C_, D_)
plot_best_fit(d18O_data[t_start_ind_:t_end_ind_], ice_age[t_start_ind_:t_end_ind_], T_interpolate[t_start_ind_: t_end_ind_], acc_interpolate[t_start_ind_:t_end_ind_], A_mean, B_mean, C_mean, A_mean_, B_mean_, C_mean_, D_mean_)
# one_curve_fit(t_start_ind_, t_end_ind_, d18O_data, T_interpolate, ice_age)
# print(np.shape(d18O_data))
# print(np.shape(A))

