import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import scipy.stats as stats
import scipy
import pycbc
from pycbc import events
import event_distr_functions as func
import itertools
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib import cm
from colorspacious import cspace_converter
from collections import OrderedDict
from matplotlib.lines import Line2D
plt.rcParams.update({
    "text.usetex": True})

class trigger_params:
	def __init__(self, ind, snr):
		self.ind = ind
		self.snr = snr

def get_og_indices(a,w):
	avg = []
	for i in range(int(len(a)/w)):
		temp = np.average(a[i*w:(i+1)*w])
		avg.append(temp)
	index = np.argsort(avg)[::-1]
	return index 	

def get_triggers(a, threshold):
	trigger_indices = []
	trigger_snrs = []
	for i in range(len(a)):
		if (a[i] >= threshold):
			trigger_indices.append(i)
			trigger_snrs.append(a[i])
	return trigger_indices, trigger_snrs

w = 16
axes = gridspec.GridSpec(3,1)
plt.figure()
threshold = 4.0
srate = 2048.0
xaxis = np.linspace(0.0, 128.0, 262144)

a = np.loadtxt('zoomed/snr_1_0')
b = np.loadtxt('zoomed/avg_16_0')
c = np.loadtxt('zoomed/triggers_16_0')

ax1 = plt.subplot(axes[0])
ax1.grid()
ax1.plot(xaxis, a)
trigger_indices, trigger_snrs = get_triggers(a, threshold)
index_times = [x/srate for x in trigger_indices]
ax1.set_xlim([58.984, 59.082])
ax1.set_ylim([0, 5.0])
ax1.axvline(x=(trigger_indices[48]-w/2)/srate, color='lime', linestyle = '--')
ax1.axvline(x=(trigger_indices[48]+w/2)/srate, color = 'lime', linestyle = '--')
ax1.plot(index_times, trigger_snrs, 'x', color = 'red')
print(trigger_indices[48])

avg = []
for i in range(len(b)):
	for k in range(w):
		avg.append(b[i])

ax2 = plt.subplot(axes[1])
ax2.grid()
ax2.set_xlim([58.984, 59.082])
ax2.set_ylim([0, 5.0])
ax2.plot(xaxis, avg)
ind = np.zeros(w)
temp_avg = np.zeros(w)
for i in range(w):
	ind[i] = trigger_indices[48] - int(w/2) + i
	temp_avg[i] = avg[int(ind[i])]
ax2.plot(ind/srate, temp_avg, color='red')
ax2.set_ylabel('SNR')

zoomed = np.zeros(w) 
for i in range(w):
	zoomed[i] = c[int(ind[i])]

ax3 = plt.subplot(axes[2])
ax3.grid()
ax3.set_xlim([58.984, 59.082])
ax3.set_ylim([0, 5.0])
ax3.plot(ind/srate, zoomed)
ax3.axvline(x=(trigger_indices[48]-w/2)/srate, color='lime', linestyle = '--')
ax3.axvline(x=(trigger_indices[48]+w/2)/srate, color = 'lime', linestyle = '--')

plt.xlabel('time (s)')
plt.savefig('zoomed.png', dpi = 600)
plt.show()
