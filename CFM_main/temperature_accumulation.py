import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

stps = 1000
t = np.arange(-stps, 1, 1)
base = 240
f1 = 1./20000
f2 = 1./5000
f3 = 1./2000
f4 = 1./500
f5 = 1./100
f6 = 1./1000

temps = 7 * np.sin(f1 * t * 2*np.pi ) + 7 * np.cos(f2 * t* 2*np.pi) + 10 * np.sin(f3 * t* 2*np.pi - 170) + \
        7 * np.cos(f4 * t* 2*np.pi + 100) + 2 * np.sin(f5 * t* 2*np.pi) - 10 * np.cos(f6 * t * 2*np.pi + 150) + base
acc = np.exp(-21.492 + 0.0811 * temps)

input_temps = np.array([t, temps])
input_acc = np.array([t, acc])

np.savetxt('CFMinput/solver_test_T.csv', input_temps, delimiter=",")
np.savetxt('CFMinput/solver_test_Acc.csv', input_acc, delimiter=",")


plot = True
if plot:
    fig, ax = plt.subplots()
    ax.plot(t, temps, color='blue')
    ax.set_xlabel("time [yr]")
    ax.set_ylabel("Temperature [K]", color="blue")
    ax2 = ax.twinx()
    ax2.plot(t, acc, color="orange")
    ax2.set_ylabel("Accumulation ice equivalent [m/yr]", color="orange")
    plt.show()
