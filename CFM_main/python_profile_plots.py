import matplotlib.pyplot as plt
import numpy as np

times_ns = np.loadtxt('all_times.txt')
times_s = times_ns * 10 ** (-9)


plt.plot(times_s)
plt.xlabel('calls')
plt.ylabel('time per call [s]')
plt.title('time for w')
plt.text(x=1.4e-3, y=0.045, s='total time: ' + str(round(sum(times_s), 2)) + ' s')
plt.show()

