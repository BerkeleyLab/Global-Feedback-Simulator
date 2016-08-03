#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt

# bw = 10.^[3:.05:4.7];  # Hz, zero-dB feedback
bw = np.logspace(3, 4.7, num=34)  # Hz, zero-dB feedback
T = 1/(12*bw)  # s, allowed feedback delay
T0 = 1e-6  # s, hardware delay
fbw = 1/(T-T0)/2/np.pi  # Hz, low-pass shaping filter bandwidth in feedback
fnbw = fbw*np.pi/2  # Hz, noise bandwidth of above
pgain = bw/16  # proportional (amplitude) gain
drive = 0.04  # rms amplitude of drive
npd = np.square(drive)/np.square(pgain)/fnbw

plt.semilogx(bw*1e-3, np.log10(npd)*10, linewidth=3)
plt.xlim([1, 50])
# plt.title(title, fontsize=40, y=1.02)
plt.ylabel('required NPD [dBc/Hz]', fontsize=40)
plt.xlabel('zero-dB crossing [kHz]', fontsize=40)
plt.legend(loc='upper right')
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.rc('font', **{'size': 20})


bw = np.array([20, 30, 40])*1e3
T = 1/(12*bw)
fbw = 1/(T-T0)/2/np.pi
fnbw = fbw*np.pi/2
pgain = bw/16
drive = 0.04
npd = np.square(drive)/np.square(pgain)/fnbw
print npd

plt.semilogx(bw*1e-3, np.log10(npd)*10, 'o', markersize=15)
plt.text(bw[0]*1e-3*1.0, np.log10(npd[0])*10+20, 'nominal', verticalalignment='top',
         horizontalalignment='center', rotation=90, fontsize=50)
plt.text(bw[1]*1e-3*1.0, np.log10(npd[1])*10+21.5, 'HoBiCaT', verticalalignment='top',
         horizontalalignment='center', rotation=90, fontsize=50)
plt.text(bw[2]*1e-3*1.0, np.log10(npd[2])*10+23, 'high-gain', verticalalignment='top',
         horizontalalignment='center', rotation=90, fontsize=50)

plt.show()
