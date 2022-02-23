import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

start_year = 100000

# ----------------------------------------------------------------------------------------------------------------------
# Import temperature and accumulation from NGRIP

file_location = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'
df1 = pd.read_excel(file_location, sheet_name='Sheet4')
df2 = pd.read_excel(file_location, sheet_name='Sheet5')
df3 = pd.read_excel(file_location, sheet_name='Sheet6')

age0 = np.array(df1[df1.columns[2]])
t_start_ind = np.min(np.where(age0 >= start_year))
print(t_start_ind)

depth = np.array(df1[df1.columns[0]])[:t_start_ind]
age = np.array(df1[df1.columns[2]])[:t_start_ind]
t = np.flipud(age) * (-1)
accs = np.flipud(np.array(df1[df1.columns[3]]))[:t_start_ind]
temps = np.flipud(np.array(df1[df1.columns[4]]) + 273.15)[:t_start_ind]
temps_err = np.flipud(np.array(df1[df1.columns[5]]))[:t_start_ind]
temps_lo = temps - temps_err
temps_up = temps + temps_err


# Save in CFM format
input_temps = np.array([t, temps])
input_temps_lo = np.array([t, temps_lo])
input_temps_up = np.array([t, temps_up])
input_acc = np.array([t, accs])

np.savetxt('../../CFM_main/CFMinput/NGRIP_T_100kyr.csv', input_temps, delimiter=",")
np.savetxt('../../CFM_main/CFMinput/NGRIP_T_lo_100kyr.csv', input_temps_lo, delimiter=",")
np.savetxt('../../CFM_main/CFMinput/NGRIP_T_up_100kyr.csv', input_temps_up, delimiter=",")
np.savetxt('../../CFM_main/CFMinput/NGRIP_Acc_100kyr.csv', input_acc, delimiter=",")


# Plot temperature and accumulation
plot = True
if plot:
    fig, ax = plt.subplots()
    ax.plot(t/1000, temps, color='blue')
    ax.set_xlabel("Age GICC05modelext [kyr]")
    ax.set_ylabel("Temperature [K]", color="blue")
    ax2 = ax.twinx()
    ax2.plot(t/1000, accs, color="orange")
    ax2.set_ylabel("Accumulation ice equivalent [m/yr]", color="orange")
    plt.show()

