from __future__ import division
import numpy as np
import scipy as sp
from scipy.optimize import leastsq
from matplotlib import animation
import matplotlib.pyplot as plt
import h5py
import time
import os.path
import string
import shutil
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


class CfmPlotter():

    def __init__(self, fpath=None):
        # fpath = "./DO_results/DO_tests_vary_tr_time/cfm_DO_trtime_1500/Goujon_DO_trtime_1500.hdf5"
        self.fpath = fpath
        f = h5py.File(fpath)
        self.fs = os.path.split(fpath)[1]
        # self.fs = string.replace(self.fs, "_", "-")
        print(self.fs)
        print(f.keys())
        self.z = f["depth"][:]
        self.rho = f["density"][:]
        self.temperature = f["temperature"][:]
        self.age = f["age"][:]
        self.climate = f["Modelclimate"][:]
        # self.iso_sigmaD = f["iso_sigmaD"][:]
        # self.iso_sigma18 = f["iso_sigma18"][:]
        # self.iso_sigma17 = f["iso_sigma17"][:]
        self.model_time = np.array(([a[0] for a in self.z[:]]))
        self.close_off_depth = f["BCO"][:, 2]
        # self.model_time = np.array(([a[0] for a in self.rho[:]]))
        self.d15N = (f["d15N2"][:] - 1.) * 1000

        f.close()
        return

    def init_plot(self):
        self.f0 = plt.figure(num=0, figsize=(11, 7), dpi=200)
        self.f0.tight_layout(pad=2.8)
        # self.f0.suptitle("CFM diffusion" , fontsize=12)
        self.ax01 = plt.subplot2grid((2, 3), (0, 1))
        self.ax02 = plt.subplot2grid((2, 3), (1, 1))
        self.ax03 = plt.subplot2grid((2, 3), (0, 0))
        self.ax04 = plt.subplot2grid((2, 3), (1, 0))
        self.ax05 = plt.subplot2grid((2, 3), (0, 2))
        self.ax06 = plt.subplot2grid((2, 3), (1, 2))

        self.ax01.set_ylim((250, 0))
        self.ax02.set_ylim((250, 0))
        self.ax03.set_ylim((225, 250))
        self.ax04.set_ylim((0.15, 0.23))
        self.ax02.set_xlim((225, 250))
        self.ax03.set_xlim((self.model_time[0], self.model_time[-1]))
        self.ax04.set_xlim((self.model_time[0], self.model_time[-1]))
        self.ax05.set_ylim((250, 0))
        self.ax06.set_ylim((130, 60))
        self.ax06.set_xlim((self.model_time[0], self.model_time[-1]))

        self.ax01.set_ylabel(r"Depth [m]")
        self.ax01.set_xlabel(r"Density [$\mathrm{kgm}^{-3}$]", labelpad=-1.5, fontsize=9)  # labelpad=-1
        self.ax02.set_ylabel(r"Depth [m]")
        self.ax02.set_xlabel(r"Temperature [K]", labelpad=-1.5, fontsize=9)
        self.ax03.set_ylabel(r"Temperature Forcing [K]")
        self.ax03.set_xlabel(r"Model Time [y]", labelpad=-1.5, fontsize=9)
        self.ax04.set_ylabel(r"Accumulation Forcing [$\mathrm{my}^{-1}$ ice eq.]")
        self.ax04.set_xlabel(r"Model Time [y]", labelpad=-1.5, fontsize=9)
        self.ax05.set_ylabel(r"Depth [m]")
        self.ax05.set_xlabel(r"$\delta^{15}N$ [‰]", labelpad=-1.5, fontsize=9)
        self.ax06.set_ylabel(r"Close-off depth [m]")
        self.ax06.set_xlabel(r"Model Time [y]", labelpad=-1.5, fontsize=9)

        #self.ax01.set_title('Density profile', fontsize=10)
        #self.ax02.set_title('CFM diffusion', fontsize=10)
        #self.ax03.set_title('Temperature Forcing', fontsize=10)
        #self.ax04.set_title("Accumulation Forcing", fontsize=10)
        #self.ax05.set_title('Isotope Signal', fontsize=10)
        #self.ax06.set_title('Close-off depth', fontsize=10)

        self.ax01.grid(linestyle=':', color='gray', lw='0.3')
        self.ax02.grid(linestyle=':', color='gray', lw='0.3')
        self.ax03.grid(linestyle=':', color='gray', lw='0.3')
        self.ax04.grid(linestyle=':', color='gray', lw='0.3')
        self.ax05.grid(linestyle=':', color='gray', lw='0.3')
        self.ax06.grid(linestyle=':', color='gray', lw='0.3')


        self.p0011, = self.ax01.plot(self.rho[500][1:], self.z[500][1:], 'b-', linewidth=0.7, label='$t_0$')
        self.p0112, = self.ax01.plot(self.rho[1500][1:], self.z[1500][1:], 'b-', alpha=0.7, linewidth=0.7, label='$t_1$')
        self.p0113, = self.ax01.plot(self.rho[2000][1:], self.z[2000][1:], 'b-', alpha=0.6, linewidth=0.7, label='$t_2$')
        self.p0114, = self.ax01.plot(self.rho[2500][1:], self.z[2500][1:], 'b-', alpha=0.5, linewidth=0.7, label='$t_3$')
        self.p0115, = self.ax01.plot(self.rho[3000][1:], self.z[3000][1:], 'b-', alpha=0.3, linewidth=0.7, label='$t_4$')
        self.p0116, = self.ax01.plot(self.rho[4000][1:], self.z[4000][1:], 'b-', alpha=0.1, linewidth=0.7, label='$t_5$')
        self.p0117, = self.ax01.plot(self.rho[4990][1:], self.z[4990][1:], 'b-', alpha=0.1, linewidth=0.7, label='$t_6$')
        self.ax01.legend(loc='lower left', fontsize=8)

        self.p0121, = self.ax02.plot(self.temperature[500][1:], self.z[500][1:], 'b-', linewidth=0.7, label='$t_0$')
        self.p0122, = self.ax02.plot(self.temperature[1500][1:], self.z[1500][1:], 'b-', alpha=0.7, linewidth=0.7, label='$t_1$')
        self.p0123, = self.ax02.plot(self.temperature[2000][1:], self.z[2000][1:], 'b-', alpha=0.6, linewidth=0.7, label='$t_2$')
        self.p0124, = self.ax02.plot(self.temperature[2500][1:], self.z[2500][1:], 'b-', alpha=0.5, linewidth=0.7, label='$t_3$')
        self.p0125, = self.ax02.plot(self.temperature[3000][1:], self.z[3000][1:], 'b-', alpha=0.3, linewidth=0.7, label='$t_4$')
        self.p0126, = self.ax02.plot(self.temperature[4000][1:], self.z[4000][1:], 'b-', alpha=0.2, linewidth=0.7, label='$t_5$')
        self.p0127, = self.ax02.plot(self.temperature[4990][1:], self.z[4990][1:], 'b-', alpha=0.1, linewidth=0.7, label='$t_6$')
        self.ax02.legend(loc='lower left', fontsize=8)


        # self.p021, = self.ax03.plot(self.climate[0, 0], self.climate[0, 1], 'k-')
        # self.p022, = self.ax04.plot(self.climate[0, 0], self.climate[0, 2], 'k-')
        # TODO: changed to self.climate[:, 0] to select the whole column...
        self.p021, = self.ax03.plot(self.climate[:, 0], self.climate[:, 2], 'k-')
        self.ax03.axvline(x=-4500, linestyle='--', linewidth=0.5)
        self.ax03.text(x=-4400, y=227, s='$t_0$')
        self.ax03.axvline(x=-3500, linestyle='--', linewidth=0.5)
        self.ax03.text(x=-3400, y=227, s='$t_1$')
        self.ax03.axvline(x=-3000, linestyle='--', linewidth=0.5)
        self.ax03.text(x=-2900, y=227, s='$t_2$')
        self.ax03.axvline(x=-2500, linestyle='--', linewidth=0.5)
        self.ax03.text(x=-2400, y=227, s='$t_3$')
        self.ax03.axvline(x=-2000, linestyle='--', linewidth=0.5)
        self.ax03.text(x=-1900, y=227, s='$t_4$')
        self.ax03.axvline(x=-1000, linestyle='--', linewidth=0.5)
        self.ax03.text(x=-900, y=227, s='$t_5$')
        self.ax03.axvline(x=-30, linestyle='--', linewidth=0.5)
        self.ax03.text(x=-300, y=227, s='$t_6$')

        self.ax06.axvline(x=-4500, linestyle='--', linewidth=0.5)
        self.ax06.text(x=-4400, y=127, s='$t_0$')
        self.ax06.axvline(x=-3500, linestyle='--', linewidth=0.5)
        self.ax06.text(x=-3400, y=127, s='$t_1$')
        self.ax06.axvline(x=-3000, linestyle='--', linewidth=0.5)
        self.ax06.text(x=-2900, y=127, s='$t_2$')
        self.ax06.axvline(x=-2500, linestyle='--', linewidth=0.5)
        self.ax06.text(x=-2400, y=127, s='$t_3$')
        self.ax06.axvline(x=-2000, linestyle='--', linewidth=0.5)
        self.ax06.text(x=-1900, y=127, s='$t_4$')
        self.ax06.axvline(x=-1000, linestyle='--', linewidth=0.5)
        self.ax06.text(x=-900, y=127, s='$t_5$')
        self.ax06.axvline(x=-30, linestyle='--', linewidth=0.5)
        self.ax06.text(x=-300, y=127, s='$t_6$')


        self.p022, = self.ax04.plot(self.climate[:, 0], self.climate[:, 1], 'k-')

        self.p0231, = self.ax05.plot(self.d15N[500][1:], self.z[500][1:], 'b-', linewidth=0.7, label='$t_0$')
        self.p0232, = self.ax05.plot(self.d15N[1500][1:], self.z[1500][1:], 'b-', alpha=0.7, linewidth=0.7, label='$t_1$')
        self.p0233, = self.ax05.plot(self.d15N[2000][1:], self.z[2000][1:], 'b-', alpha=0.6, linewidth=0.7, label='$t_2$')
        self.p0234, = self.ax05.plot(self.d15N[2500][1:], self.z[2500][1:], 'b-', alpha=0.5, linewidth=0.7, label='$t_3$')
        self.p0235, = self.ax05.plot(self.d15N[3000][1:], self.z[3000][1:], 'b-', alpha=0.3, linewidth=0.7, label='$t_4$')
        self.p0236, = self.ax05.plot(self.d15N[4000][1:], self.z[4000][1:], 'b-', alpha=0.2, linewidth=0.7, label='$t_5$')
        self.p0237, = self.ax05.plot(self.d15N[4990][1:], self.z[4990][1:], 'b-', alpha=0.1, linewidth=0.7, label='$t_6$')

        self.ax05.legend(loc='lower left', fontsize=8)
        # self.iso_sigma18_co = np.array((self.iso_sigma18[0][1:][self.rho[0][1:]>804.3][0],))
        self.p0251, = self.ax06.plot(self.z[:, 0], self.close_off_depth, 'b-', alpha=1, linewidth=0.7, label='$t_0$')

        plt.tight_layout()
        # TODO: added plt.tight_layout; this helps a bit, however labels and headers still overlap
        plt.savefig('results/densification/ramp_Goujon.pdf')
        plt.show()
        return

    def update_data(self, i):
        self.p011.set_data(self.rho[i][1:], self.z[i][1:])
        self.p012.set_data(self.temperature[i][1:], self.z[i][1:])
        self.p021.set_data(self.climate[:i, 0], self.climate[:i, 2])
        self.p022.set_data(self.climate[:i, 0], self.climate[:i, 1])
        self.p023.set_data(self.rho[i][1:], self.z[i][1:])
        #self.p024.set_data(self.rho[i][1:], self.z[i][1:])
        #self.iso_sigma18_co = np.append(self.iso_sigma18_co, self.iso_sigma18[i][1:][self.rho[i][1:] > 804.3][0])
        #self.p025.set_data(self.climate[:i, 0], self.iso_sigma18_co[:i])
        self.p025.set_data(self.z[:i, 0], self.close_off_depth[:i])
        #self.f0.suptitle(r"Model Phys: %s  Model time: %0.2f" % (self.fs, self.z[i][0]))
        #self.f0.suptitle('Model Phys ' + str(self.fs) + 'Model time: ' + str(self.z[i][0]))
        #TODO: include the modeltime and model physics... this really freaks me out
        self.f0.tight_layout(pad=2.8)
        return self.p011, self.p012, self.p021, self.p022, self.p023, self.p025  #,self.p024,

    def plot_final(self):
        self.init_plot()
        self.p011, = self.ax01.plot(self.rho[0][1:], self.z[0][1:], 'b-')
        self.p012, = self.ax02.plot(self.temperature[0][1:], self.z[0][1:], 'k-')
        self.p021, = self.ax03.plot(self.climate[:, 0], self.climate[:, 2], 'k-')
        self.p022, = self.ax04.plot(self.climate[:, 0], self.climate[:, 1], 'k-')  #
        self.p023, = self.ax05.plot(self.rho[0][1:], self.z[0][1:], 'r-')
        #self.p024, = self.ax05.plot(self.rho[0][1:], self.z[0][1:], 'b-')
        # self.iso_sigmaD_co = np.array((self.iso_sigmaD[0][self.rho[0]>804.3][0],))
        #iso_sigma18_co = np.array(([self.iso_sigma18[i][1:][self.rho[i][1:] >= 804.3][0] for i in range(len(self.z))]))
        # self.iso_sigmaD_co = np.array(([self.iso_sigmaD[j][self.rho[j]>804.3][0] for j in np.arange(size(self.))]))
        #self.p025, = self.ax06.plot(self.climate[:, 0], iso_sigma18_co, 'b-')  #
        self.p025, = self.ax05.plot(self.z[:, 0], self.close_off_depth, 'b-')
        plt.show()
        return

    # def export_ascii(self, fout = None):
    #     """
    #     exports ascii with temp, accum history, close off depth and age as well as diffusion
    #     lengths at CO
    #     """
    #     model_time = np.array(([a[0] for a in self.z[:]]))
    #     temp_forcing = self.climate[:,2]
    #     accum_forcing = self.climate[:,1]
    #
    #     depth_co = np.array(([self.z[i][1:][self.rho[i][1:]>=804.3][0] for i in range(len(self.z))]))
    #     age_co = np.array(([self.age[i][1:][self.rho[i][1:]>=804.3][0] for i in range(len(self.z))]))
    #     sigma17_co = np.array(([self.iso_sigma17[i][1:][self.rho[i][1:]>=804.3][0] for i in range(len(self.z))]))
    #     sigma18_co = np.array(([self.iso_sigma18[i][1:][self.rho[i][1:]>=804.3][0] for i in range(len(self.z))]))
    #     sigmaD_co = np.array(([self.iso_sigmaD[i][1:][self.rho[i][1:]>=804.3][0] for i in range(len(self.z))]))
    #
    #     dataout = np.transpose(np.vstack((model_time, temp_forcing, accum_forcing, depth_co, age_co, sigma17_co,\
    #         sigma18_co, sigmaD_co)))
    #     dataout = np.delete(dataout, 0, 0)
    #
    #     if fout:
    #         f = open(fout, "w")
    #     else:
    #         f = open("./ramp_animation.out", "w")
    #     f.write("model_time\ttemp_forcing\tacuum_forcing\tdepth_co\tage_co\tsigma17_co\tsigma18_co\tsigmaD_co\n")
    #     np.savetxt(f, dataout, fmt = "%0.1f\t%0.2f\t%0.4e\t%0.3f\t%0.1f\t%0.6e\t%0.6e\t%0.6e")
    #     f.close()
    #
    #
    #     return


t1 = time.time()
# #plt.ion()
# folder_name = "./DO_results"
# list_of_files = os.listdir(folder_name)
#
# for j in list_of_files:
#     condition_1 = ("RC_tau" in j) & ("hdf5" in j) & ("Goujon" in j)  &\
#     ("acc300" in j)
#     if condition_1:
#         try:
#             t2 = time.time()
#             #Ascii out block
#             print("\n\nReading file %s" %folder_name + "/" + j)
#             plot_class = CfmPlotter(fpath = folder_name + "/" + j)
#             plot_class.export_ascii()
#             os.rename("./ramp_animation.out", "./" + os.path.splitext(j)[0] + ".out")
#             print("Saving file %s" %("./" + os.path.splitext(j)[0] + ".out"))
#
#             #Plotting block
#             plot_class.init_plot()
#             print len(plot_class.z)
#             anim = animation.FuncAnimation(plot_class.f0, plot_class.update_data, frames = np.int(len(plot_class.z)), interval = 10, blit = True)
#             anim.save('./basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
#             os.rename("./basic_animation.mp4", "./" + os.path.splitext(j)[0] + ".mp4")
#             print("Saving file %s" %("./" + os.path.splitext(j)[0] + ".mp4"))
#             print("Loop time %0.1f minutes." %((time.time() - t2)/60.))
#
#
#         except IOError:
#             print("\nError reading %s " %folder_name + "/" + j)
#             continue
#     else:
#         print("Not processing file  %s" %folder_name + "/" + j)
#         continue


plt.close("all")

plot_class = CfmPlotter(fpath="results/densification/CFMresults_ramp_T_Goujon.hdf5")
# plot_class.export_ascii()

plot_class.init_plot()

# print len(plot_class.z)
#anim = animation.FuncAnimation(plot_class.f0, plot_class.update_data, frames=(len(plot_class.z)), interval=5, blit=True)

#anim.save('./basic_animation.mp4', fps=50, extra_args=['-vcodec', 'libx264'])

plt.close("all")
#anim.save('./basic_animation.mp4', fps=50)



#plt.show()

print("Script time %0.1f minutes." % ((time.time() - t1) / 60.))

