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


def fit_curve(cutoff, target):
	z = np.polyfit(cutoff[140:], target[140:], 2)
	p = np.poly1d(z)
	xp = np.linspace(1, 5, 500)
	plt.plot(cutoff, target, '*', xp, p(xp), '--')	
	#plt.show()
	print("YOLO",p(4.0))

def cutoff_vs_triggers(cutoff, counts, window, marker, colors, N, index):
	#value, indices = np.unique(target, return_index=True)
	plt.plot(cutoff, counts/N*2048, linestyle=" ", color = colors[index], marker=marker[index], label = window)	
	#plt.plot(cutoff[indices], value, linestyle=" ", marker=marker[index], label = window)
	#plt.ylabel('No. of 1st stage triggers/sec')
	plt.ylabel('$f (\omega,$'r'$\rho_{cut}) sec^{-1}$')
	plt.xlabel(r'$\rho_{cut}$')
	plt.yscale('log')


def cutoff_vs_target(cutoff, target, window, marker, colors, index):
	value, indices = np.unique(target, return_index=True)
	#plt.plot(cutoff, target, linestyle=" ", marker=marker[index], label = window)	
	plt.plot(cutoff[indices], value, linestyle=" ", color = colors[index], marker=marker[index], label = window)
	plt.ylim(1.8, 6)
	plt.ylabel('Target SNR')
	plt.xlabel(r'$\rho_{cut}$')

def compute_cost(axes, pca, N, templates, window, cutoff, target, triggers, marker, colors, index):
	total = []
	second = []
	first = []
	fft_cost = 2*templates*N*(np.log(N)+1)/10**15

	for i in range(len(triggers)):
		ein = (pca*N*templates/window + 2*pca*N*np.log(N) )/10**15 
		zwei = (pca * window * triggers[i] * templates)/10**15 
		ganz = (pca * N * templates/window + pca * window * triggers[i] * templates + pca*N*np.log(N) + 5*pca*N)/10**15 
		first.append(ein)
		second.append(zwei)
		total.append(ganz)
	
	ax1 = plt.subplot(axes[0,:])
	ax1.grid()
	plt.plot(target, first, linestyle = '-', color = colors[index], marker = marker[index])
	#if (window == 16):
		#plt.axhline(y = fft_cost, color = 'deeppink', label = 'FFT')

	ax2 = plt.subplot(axes[0,:])
	#ax2.set_yscale('log')
	#ax2.grid()
	plt.plot(target, second, linestyle = '--', color = colors[index])
	if (window == 16):
		plt.axhline(y = fft_cost, color = 'darkorange', label = 'FFT')

	ax3 = plt.subplot(axes[1,:])
	ax3.grid()
	ax3.set_yscale('log')
	plt.plot(target, total, linestyle = '', color = colors[index], marker = marker[index], label=window)
	if (window == 16):
		plt.axhline(y = fft_cost, color = 'darkorange', label = 'FFT')
	print("Min cost for window", window, np.min(total))
	print("FFT cost", fft_cost)
	print("Reduction of", np.round((1-abs(np.min(total)/fft_cost))*100, 2))	
	
	ax1.set_ylabel('# operations / $10^{15}$')
	ax3.set_ylabel('# operations / $10^{15}$')
	ax3.legend()
	plt.xlabel('Target SNR')
	return ax1, ax2, ax3


def read_conv_SNR(filename):
	df = pd.read_csv(filename, sep = " ")
	cutoff = df["cutoff"]
	target = df['target']
	count = df['count']
	return cutoff, target, count	


def legend_for_stack_plot(ax2, ax3):
	line_style = ['-', '--']
	line_marker = ['o']
	lines=[]
	for i in range(2):
		if (i>0):
			lines.append(Line2D([0], [0], color='blue', linewidth=1, linestyle=line_style[i]))
		else:
			lines.append(Line2D([0], [0], color='blue', marker = line_marker[0], linewidth=1, linestyle=line_style[i]))
	labels = ['First stage', 'Second stage']
	ax2.legend(lines, labels, loc = 9)
	ax1.grid()
	ax3.grid()



axes = gridspec.GridSpec(2, 2)
plt.figure()
window_array = [ 2, 4, 8, 16 ]
#window_array = [16]
colors = [cm.Set1(x) for x in range(len(window_array))]
no_realizations = 5000
marker = ['+', 'x', 'o', '*', '^'] 
index = 0
for k in window_array:
	window = k
	pca = 254 
	N = 2048*128*no_realizations
	templates = 6250 
	filename = "conv_SNRs/best_bank/%s/conv_SNR_%s_%s.csv" % (window, window, no_realizations)
	
	cutoff, target, count = read_conv_SNR(filename)

	#cutoff_vs_triggers(cutoff, count, window, marker, colors, N, index)
	#cutoff_vs_target(cutoff, target, window, marker, colors, index)
	ax1, ax2, ax3 = compute_cost(axes, pca, N, templates, window, cutoff, target, count, marker, colors, index)
	index += 1
	#fit_curve(cutoff, target)

legend_for_stack_plot(ax2, ax3)

#plt.grid()
plt.legend()
#plt.savefig('Costs_vs_target_snr.png', dpi = 600)
#plt.savefig('HIGH_cutoff_vs_target.png', dpi = 600)
#plt.savefig('Triggers_vs_cutoff.png', dpi=600)
plt.savefig('Costs_stack.png', dpi=600)
plt.show()
