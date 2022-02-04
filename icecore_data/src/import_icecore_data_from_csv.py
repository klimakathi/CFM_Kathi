import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ----------------------------------------------------------------------------------------------------------------------
# Import temperature and accumulation from NGRIP

file_location = '~/projects/Thesis/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'
df1 = pd.read_excel(file_location, sheet_name='Sheet4')
df2 = pd.read_excel(file_location, sheet_name='Sheet5')
df3 = pd.read_excel(file_location, sheet_name='Sheet6')

depth = np.array(df1[df1.columns[0]])
age = np.array(df1[df1.columns[2]])
t = np.flipud(age) * (-1)
accs = np.flipud(np.array(df1[df1.columns[3]]))
temps = np.flipud(np.array(df1[df1.columns[4]]) + 273.15)


# Save in CFM format
input_temps = np.array([t, temps])
input_acc = np.array([t, accs])

np.savetxt('../../../../CFM_main/CFMinput/NGRIP_T.csv', input_temps, delimiter=",")
np.savetxt('../../../../CFM_main/CFMinput/NGRIP_Acc.csv', input_acc, delimiter=",")


# Plot temperature and accumulation
plot = False
if plot:
    fig, ax = plt.subplots()
    ax.plot(t/1000, temps, color='blue')
    ax.set_xlabel("Age GICC05modelext [kyr]")
    ax.set_ylabel("Temperature [K]", color="blue")
    ax2 = ax.twinx()
    ax2.plot(t/1000, accs, color="orange")
    ax2.set_ylabel("Accumulation ice equivalent [m/yr]", color="orange")
    plt.show()

