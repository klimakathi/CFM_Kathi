import matplotlib.pyplot as plt
from scipy.optimize import least_squares, basinhopping, minimize, leastsq
from read_d18O import read_data_d18O, get_interval_data
from read_temp_acc import *
from smoothing_splines import smooth_data

data_path = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/interpolated_data.xlsx'


def quadratic_func(d18O_data, theta):
    a = theta[0]
    b = theta[1]
    c = theta[2]
    T = a * d18O_data ** 2 + b * d18O_data + c
    return T


def cubic_func(d18O_data, theta):
    a = theta[0]
    b = theta[1]
    c = theta[2]
    d = theta[3]
    T = a * d18O_data ** 3 + b * d18O_data ** 2 + c * d18O_data + d
    return T


def fun_quadratic(theta):
    a = theta[0]
    b = theta[1]
    c = theta[2]
    temps_model = a * d18o_smooth ** 2 + b * d18o_smooth + c
    return 1 / (np.shape(temps_interval)[0] - 1) * np.sqrt(np.sum((temps_interval - temps_model) ** 2))


def fun_quadratic_v(theta):
    a = theta[0]
    b = theta[1]
    c = theta[2]
    temps_model = a * d18o_smooth ** 2 + b * d18o_smooth + c
    return temps_interval - temps_model


def fun_cubic(theta):
    a = theta[0]
    b = theta[1]
    c = theta[2]
    d = theta[3]
    temps_model = a * d18o_smooth ** 3 + b * d18o_smooth ** 2 + c * d18o_smooth + d
    return 1 / (np.shape(temps_interval)[0] - 1) * np.sqrt(np.sum((temps_interval - temps_model) ** 2))


def fun_cubic_v(theta):
    a = theta[0]
    b = theta[1]
    c = theta[2]
    d = theta[3]
    temps_model = a * d18o_smooth ** 3 + b * d18o_smooth ** 2 + c * d18o_smooth + d
    return temps_interval - temps_model


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    mode = 'minimize'
    quadratic = True
    cubic = True

    start_year_ = -50000
    end_year_ = -10000
    time_grid_stp_ = 10
    cop_ = 1 / 200.

    theta_c_0 = [1000, 1000, 1000, 1000]
    theta_q_0 = [1000, 1000, 1000]

    depth_full, d18o_full, ice_age_full = read_data_d18O(data_path)
    depth_interval, d18o_interval, ice_age_interval, d18o_smooth_, time_ = get_interval_data(depth_full, d18o_full,
                                                                                             ice_age_full,
                                                                                             start_year_, end_year_,
                                                                                             time_grid_stp_,
                                                                                             cop_)
    d18o_smooth = smooth_data(cop_, d18o_interval, ice_age_interval, ice_age_interval)[0]
    temps, temps_err = read_temp(data_path)
    temps_interval, temps_err_interval = get_interval_temp(temps, temps_err, ice_age_full, start_year_, end_year_)
    # temps_interval = smooth_data(cop_, temps_interval, ice_age_interval, ice_age_interval)[0]
    # temps_err_interval = smooth_data(cop_, temps_err_interval, ice_age_interval, ice_age_interval)[0]

    plot_compare_minimizers = True
    if plot_compare_minimizers:
        fig, axs = plt.subplots(4, sharex=True, sharey=False)
        fig.set_figheight(15)
        fig.set_figwidth(15)
        fig.suptitle('Temperature fit from $\delta^{18}O$', fontsize=16)

        mode = 'least_squares'
        if mode == 'least_squares':
            if cubic and quadratic:
                res_c = least_squares(fun_cubic_v, theta_c_0, method='trf')
                theta_c_1 = res_c.x
                f_c = cubic_func(d18o_smooth, theta_c_1)
                sigma_c = 1 / np.shape(f_c)[0] * np.sum((temps_interval - f_c)**2)

                res_q = least_squares(fun_quadratic_v, theta_q_0, method='trf')
                theta_q_1 = res_q.x
                f_q = quadratic_func(d18o_smooth, theta_q_1)
                sigma_q = 1 / np.shape(f_q)[0] * np.sum((temps_interval - f_q)**2)

                axs[0].plot(ice_age_interval, temps_interval, 'bo', markersize=1, label='data')
                axs[0].plot(ice_age_interval, temps_interval + temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[0].plot(ice_age_interval, temps_interval - temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[0].fill_between(ice_age_interval, temps_interval - temps_err_interval,
                                    temps_interval + temps_err_interval,
                                    alpha=0.2, facecolor='b')
                axs[0].plot(ice_age_interval, quadratic_func(d18o_smooth, theta_q_1), 'r',
                            label='quadratic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K' % (
                                theta_q_1[0], theta_q_1[1], theta_q_1[2]))
                axs[0].plot(ice_age_interval, cubic_func(d18o_smooth, theta_c_1), 'orange',
                            label='cubic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K, d=%5.3f K' % (
                                theta_c_1[0], theta_c_1[1], theta_c_1[2], theta_c_1[3]))
                axs[0].set_title('Least squares fit "TRF" - $\sigma_q^2$ = %5.3f, $\sigma_c^2 $ = %5.3f' % (sigma_q, sigma_c))
                # axs[0].set_xlabel('GICC05modelext ice age [yr]')
                axs[0].set_ylabel('Temperature [°C]')
                axs[0].legend()
                axs[0].grid(linestyle='--', color='gray', lw='0.5')

                res_c = least_squares(fun_cubic_v, theta_c_0, method='lm')
                theta_c_1 = res_c.x
                f_c = cubic_func(d18o_smooth, theta_c_1)
                sigma_c = 1 / np.shape(f_c)[0] * np.sum((temps_interval - f_c)**2)

                res_q = least_squares(fun_quadratic_v, theta_q_0, method='lm')
                theta_q_1 = res_q.x
                f_q = quadratic_func(d18o_smooth, theta_q_1)
                sigma_q = 1 / np.shape(f_q)[0] * np.sum((temps_interval - f_q)**2)

                axs[1].plot(ice_age_interval, temps_interval, 'bo', markersize=1, label='data')
                axs[1].plot(ice_age_interval, temps_interval + temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[1].plot(ice_age_interval, temps_interval - temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[1].fill_between(ice_age_interval, temps_interval - temps_err_interval,
                                    temps_interval + temps_err_interval,
                                    alpha=0.2, facecolor='b')
                axs[1].plot(ice_age_interval, quadratic_func(d18o_smooth, theta_q_1), 'r',
                            label='quadratic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K' % (
                                theta_q_1[0], theta_q_1[1], theta_q_1[2]))
                axs[1].plot(ice_age_interval, cubic_func(d18o_smooth, theta_c_1), 'orange',
                            label='cubic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K, d=%5.3f K' % (
                                theta_c_1[0], theta_c_1[1], theta_c_1[2], theta_c_1[3]))
                axs[1].set_title('Least squares fit "LM" - $\sigma_q^2$ = %5.3f, $\sigma_c^2$ = %5.3f' % (sigma_q, sigma_c))
                # axs[0].set_xlabel('GICC05modelext ice age [yr]')
                axs[1].set_ylabel('Temperature [°C]')
                axs[1].legend()
                axs[1].grid(linestyle='--', color='gray', lw='0.5')

        mode = 'minimize_BFGS'
        if mode == 'minimize_BFGS':
            if cubic and quadratic:
                # minimizer_kwargs = {"method": "BFGS", "jac": True}
                res_c = minimize(fun_cubic, theta_c_0, method='BFGS')
                theta_c_1 = res_c.x
                f_c = cubic_func(d18o_smooth, theta_c_1)
                sigma_c = 1 / np.shape(f_c)[0] * np.sum((temps_interval - f_c)**2)

                res_q = minimize(fun_quadratic, theta_q_0, method='BFGS')
                theta_q_1 = res_q.x
                f_q = quadratic_func(d18o_smooth, theta_q_1)
                sigma_q = 1 / np.shape(f_q)[0] * np.sum((temps_interval - f_q)**2)

                axs[2].plot(ice_age_interval, temps_interval, 'bo', markersize=1, label='data')
                axs[2].plot(ice_age_interval, temps_interval + temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[2].plot(ice_age_interval, temps_interval - temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[2].fill_between(ice_age_interval, temps_interval - temps_err_interval,
                                    temps_interval + temps_err_interval,
                                    alpha=0.2, facecolor='b')
                axs[2].plot(ice_age_interval, quadratic_func(d18o_smooth, theta_q_1), 'r',
                            label='quadratic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K' % (
                                theta_q_1[0], theta_q_1[1], theta_q_1[2]))
                axs[2].plot(ice_age_interval, cubic_func(d18o_smooth, theta_c_1), 'orange',
                            label='cubic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K, d=%5.3f K' % (
                                theta_c_1[0], theta_c_1[1], theta_c_1[2], theta_c_1[3]))
                axs[2].set_title('Minimize "BFGS" fit - $\sigma_q^2$ = %5.3f, $\sigma_c^2$ = %5.3f' % (sigma_q, sigma_c))
                # axs[1].set_xlabel('GICC05modelext ice age [yr]')
                axs[2].set_ylabel('Temperature [°C]')
                axs[2].legend()
                axs[2].grid(linestyle='--', color='gray', lw='0.5')

        mode = 'basinhopping'
        if mode == 'basinhopping':
            if cubic and quadratic:
                minimizer_kwargs = {"method": "BFGS"}

                res_c = basinhopping(fun_cubic, theta_c_0, minimizer_kwargs=minimizer_kwargs)
                theta_c_1 = res_c.x
                f_c = cubic_func(d18o_smooth, theta_c_1)
                sigma_c = 1 / np.shape(f_c)[0] * np.sum((temps_interval - f_c)**2)

                res_q = basinhopping(fun_quadratic, theta_q_0, minimizer_kwargs=minimizer_kwargs)
                theta_q_1 = res_q.x
                f_q = quadratic_func(d18o_smooth, theta_q_1)
                sigma_q = 1 / np.shape(f_q)[0] * np.sum((temps_interval - f_q)**2)

                axs[3].plot(ice_age_interval, temps_interval, 'bo', markersize=1, label='data')
                axs[3].plot(ice_age_interval, temps_interval + temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[3].plot(ice_age_interval, temps_interval - temps_err_interval, 'b-', linewidth=0.5, alpha=0.5)
                axs[3].fill_between(ice_age_interval, temps_interval - temps_err_interval,
                                    temps_interval + temps_err_interval,
                                    alpha=0.2, facecolor='b')
                axs[3].plot(ice_age_interval, quadratic_func(d18o_smooth, theta_q_1), 'r',
                            label='quadratic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K' % (
                                theta_q_1[0], theta_q_1[1], theta_q_1[2]))
                axs[3].plot(ice_age_interval, cubic_func(d18o_smooth, theta_c_1), 'orange',
                            label='cubic fit: a=%5.3f K, b=%5.3f K, c=%5.3f K, d=%5.3f K' % (
                                theta_c_1[0], theta_c_1[1], theta_c_1[2], theta_c_1[3]))
                axs[3].set_title('Basinhopping "BFGS" fit - $\sigma_q^2$ = %5.3f, $\sigma_c^2$ = %5.3f' % (sigma_q, sigma_c))
                axs[3].set_xlabel('GICC05modelext ice age [yr]')
                axs[3].set_ylabel('Temperature [°C]')
                axs[3].legend()
                axs[3].grid(linestyle='--', color='gray', lw='0.5')

        plt.tight_layout
        plt.show()




