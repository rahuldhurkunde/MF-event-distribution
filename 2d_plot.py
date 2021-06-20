import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import scipy.stats as stats
import pycbc
from pycbc import events
import event_distr_functions as func
import itertools

def fit_curve(cutoff, target, window):
	z = np.polyfit(cutoff, triggers, 2)
	p = np.poly1d(z)
	xp = np.linspace(0.0, 3.0, 500)
	plt.plot(target, cutoff, label = window)	
	#plt.plot(cutoff, triggers, '.', xp, p(xp), '--')
	#plt.yscale('log')
	#plt.show()
	return p

def compute_cost(pca, N, templates, window, cutoff, target, triggers, marker):
	cost = []
	for i in range(len(triggers)):
		temp = (pca * N * templates/window + pca * window * triggers[i])/10**9 
		#temp = pca * N * templates / window
		#temp = pca *  triggers[i] * window
		cost.append(temp)
	#plt.plot(target, cutoff, linestyle = '', marker = next(marker), label=window)
	#plt.plot(cutoff, cost, linestyle = '', marker = next(marker), label=window)
	plt.plot(target, cost, linestyle = '', marker = next(marker), label=window)
	print("Min cost for window", window, np.min(cost))
	#plt.plot(target, cost, label=window)

window_array = [2, 4, 8, 16, 32]
marker = itertools.cycle(('+', 'o', '*')) 
for k in window_array:
	window = k
	pca = 556
	N = 2048*128
	templates = 1
	filename = "conv_SNR_%d.txt" % window
	
	cutoff = []
	target = []
	triggers = []
	a = np.loadtxt(filename)
	for i in range(len(a)):
		cutoff.append(a[i][0])
		target.append(a[i][1])
		triggers.append(a[i][2])

	#p = fit_curve(cutoff, target, window)
	compute_cost(pca, N, templates, window, cutoff, target, triggers, marker)

#plt.yscale('log')
plt.xlabel('Target SNR')
#plt.xlabel('First-step cutoff')
#plt.ylabel('First-step cutoff')
plt.ylabel('# operations / $10^9$')
plt.legend()
#plt.savefig('Costs_vs_target_snr.png', dpi = 600)
plt.savefig('temp.png', dpi = 600)
plt.show()
