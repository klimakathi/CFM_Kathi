import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = h5py.File('resultsFolder/CFMresults_NGRIP_HLd.hdf5', 'r')

d15n_model = file["d15N2"][:] - 1.
depth_model = file['depth']

file_location = '~/projects/CFM_Kathi/icecore_data/data/NGRIP/supplement.xlsx'
df3 = pd.read_excel(file_location, sheet_name='Sheet6')

depth_data = np.array(df3[df3.columns[2]])[:] * (-1)
d15n_data = np.array(df3[df3.columns[3]])[:]

plt.plot(depth_model, d15n_model, 'bo')
#plt.plot(depth_data, d15n_data)
plt.show()