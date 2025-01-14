import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg'
dataDir = "Main_Dynamic"
free_file = os.path.join(dataDir, "DFree.out")
time_series1 = np.loadtxt(free_file)
Drift = (time_series1[:, 1]) # 3.0 is the height between node 1 and node 2
Time = (time_series1[:, 0])
print(Drift)

plt.plot(Time, Drift, color='blue', linewidth=2, label='TimeSeries')
plt.xlabel('Time (s)')
plt.ylabel("Displacement at the Top (m)", fontsize=14)
plt.title("Time-History", fontsize=16)
plt.show()
exit()