# %%

import numpy as np
from scipy import signal as sig
from matplotlib import pyplot as plt

# %%

timestamp = '2024-01-09T09:02:05.974216'

ddl_col2name = ["Setpoint", "Metrology", "Error", "Stage", "Piezo"]
vib_col2name = ["Setpoint"]
ddl_dtype = [("Setpoint", float), ("Metrology", float), ("Error", float), ("Stage", float), ("Piezo", float)]

# %%

ddl_data = np.loadtxt(f"data/tacRecord_ddl3_{timestamp}.dat", usecols=(1, 2, 3, 4, 5), dtype=ddl_dtype)
vib_data = np.loadtxt(f"data/tacRecord_vib_{timestamp}.dat", usecols=(1,))

# %%

freq, Cxy1 = sig.coherence(ddl_data['Setpoint'], ddl_data['Metrology'], fs=4000.0)
freq, Txy1 = sig.csd(ddl_data['Setpoint'], ddl_data['Metrology'], fs=4000.0)
freq, Txx1 = sig.csd(ddl_data['Setpoint'], ddl_data['Setpoint'], fs=4000.0)
tf1 = Txy1/Txx1

freq, Cxy2 = sig.coherence(ddl_data['Setpoint'], np.roll(ddl_data['Metrology'], 1), fs=4000.0)
freq, Txy2 = sig.csd(ddl_data['Setpoint'], np.roll(ddl_data['Metrology'], 1), fs=4000.0)
freq, Txx2 = sig.csd(ddl_data['Setpoint'], ddl_data['Setpoint'], fs=4000.0)
tf2 = Txy2/Txx2

# %%

fig, axarr = plt.subplots(3, 1, sharex=True)
axarr[1].axhline(-90, color='k', ls='--', lw=1)
axarr[0].plot(freq, np.abs(tf1))
axarr[1].plot(freq, np.rad2deg(np.angle(tf1)), label="Without Manhattan delay")
axarr[2].plot(freq, Cxy1)
axarr[0].plot(freq, np.abs(tf2))
axarr[1].plot(freq, np.rad2deg(np.angle(tf2)), label="With Manhattan delay")
axarr[2].plot(freq, Cxy2)
axarr[2].set_xlabel("Frequency [Hz]")
axarr[0].set_ylabel("TF amplitude")
axarr[1].set_ylabel("TF phase [deg]")
axarr[2].set_ylabel("TF coherence")
axarr[0].set_xlim(0, 800)
axarr[1].legend()
for i in range(3):
    axarr[i].axvline(185, ls='--', color='C0')
    axarr[i].axvline(160, ls='--', color='C1')
axarr[1].text(185, 100, ' 185 Hz', ha='left', color='C0')
axarr[1].text(160, 100, '160 Hz ', ha='right', color='C1')
fig.savefig("BCDDL_TF.pdf")
plt.show()

# %%

fig, axarr = plt.subplots(5, 1, sharex=True)
for i in range(ddl_data.shape[1]):
    axarr[i].plot(ddl_data[:, i], label=ddl_col2name[i])
    axarr[i].legend(loc='upper right')
plt.show()

# %%

