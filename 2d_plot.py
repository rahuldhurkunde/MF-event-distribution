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
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter


def flop_to_speed(y1, y2):
	return 1.0/y2, 1.0/y1

def convert_axflop_to_speed(ax3, ax4):
	"""
	Update second axis according with first axis.
	"""
	y1, y2 = ax3.get_ylim()
	y1n, y2n = flop_to_speed(y1, y2)
	print('YOOOOO', y1n, y2n)
	ax4.set_ylim(y1n, y2n)
	ax4.invert_yaxis()

def fit_curve(cutoff, target):
	z = np.polyfit(cutoff[140:], target[140:], 2)
	p = np.poly1d(z)
	xp = np.linspace(1, 5, 500)
	plt.plot(cutoff, target, '*', xp, p(xp), '--')	
	#plt.show()
	print("YOLO",p(4.0))

def cutoff_vs_triggers(cutoff, counts, window, marker, colors, N, index):
	#value, indices = np.unique(target, return_index=True)
	plt.plot(cutoff, counts/N*2048, linestyle=" ", color = colors[index], marker=marker[index], label = 'w = %s' %window)	
	plt.ylabel('$f (\omega,$'r'$\rho_{cut}) sec^{-1}$')
	plt.xlabel(r'$\rho_{I}$')
	plt.yscale('log')
	plt.grid()


def cutoff_vs_target(cutoff, target, window, marker, colors, index):
	value, indices = np.unique(target, return_index=True)
	plt.plot(cutoff[indices], value, linestyle=" ", color = colors[index], marker=marker[index], label = 'w = %s' %window)
	plt.ylim(1.8, 6)
	plt.ylabel('Target SNR')
	plt.xlabel(r'$\rho_{I}$')
	plt.grid()

def compute_cost(axes, pca, N, templates, window, cutoff, target, triggers, marker, colors, index):
	total = []
	second = []
	first = []
	speedup = []
	fft_cost = templates*N*(5*np.log(N)+6)/10**15

	for i in range(len(triggers)):
		ein = templates*N/window*(5*np.log(N/window)+6+2/window)/10**15/fft_cost
		#ein = (4*pca*N*templates/window/2 + 5*pca*N*np.log(N) + 26*pca*N) / 10**15 
		zwei = (4*pca*window*triggers[i]*templates) / 10**15 /fft_cost 
		ganz = (ein + zwei) 
		first.append(ein)
		second.append(zwei)
		total.append(ganz)
		speedup.append(np.power(ganz, -1))
	
	ax1 = plt.subplot(axes[0,:])
	ax1.grid()
	plt.plot(target, first, linestyle = 'dashdot', color = colors[index])
	#if (window == 2):
		#plt.axhline(y = fft_cost, color = 'deeppink', label = 'FFT')

	ax2 = plt.subplot(axes[0,:])
	#ax2.set_yscale('log')
	#ax2.grid()
	#plt.plot(target, second, linestyle = 'dashed', color = colors[index])
	#if (window == 16):
		#plt.axhline(y = 1.0, color = 'darkorange', linewidth = 2,label = 'Template filtering')

	ax3 = plt.subplot(axes[0,:])
	ax3.plot(target, total, linestyle = '--', color = colors[index], marker = marker[index], markersize = 4,label= 'w = %s' %window)
	if (window == 16):
		ax3.axhline(y = 1.0, color = 'darkorange', label = 'Template filtering')
	
	ax4 = ax3.twinx()
	convert_axflop_to_speed(ax3,ax4)
	ax4.set_ylabel('Speedup factor')
	ax3.set_yscale('log')
	ax4.set_yscale('log')
	#ax3.yaxis.set_minor_locator(AutoMinorLocator())
	#ax3.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))	
	
	snr = 6.0 
	target = np.array(target)
	idx = (np.abs(target - snr)).argmin()
	print("min cost for window", window, total[idx])
	print("FFT cost", fft_cost)
	print("Improvement of", np.round((1/abs(total[idx])), 2), 'at SNR', snr)	
	
	ax1.set_ylabel('FLOP / $10^{15}$')
	ax3.set_ylabel('FLOP / $10^{15}$')
	ax3.legend(loc=1, prop={"size":10})
	ax3.set_xlabel('Target SNR')
	return ax1, ax2, ax3

def read_conv_SNR(filename):
	df = pd.read_csv(filename, sep = " ")
	cutoff = df["cutoff"]
	target = df['target']
	count = df['count']
	return cutoff, target, count	


def legend_for_stack_plot(ax2, ax3):
	line_style = ['dashdot', 'dashed']
	line_marker = ['o']
	lines=[]
	for i in range(2):
		if (i>0):
			lines.append(Line2D([0], [0], color='blue', linewidth=1, linestyle=line_style[i]))
		else:
			#lines.append(Line2D([0], [0], color='blue', marker = line_marker[0], linewidth=1, linestyle=line_style[i]))
			lines.append(Line2D([0], [0], color='blue', linewidth=1, linestyle=line_style[i]))
	labels = ['First stage', 'Second stage']
	ax2.legend(lines, labels, loc = 3)
	ax1.grid()
	ax3.grid()



axes = gridspec.GridSpec(1, 1)
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

#legend_for_stack_plot(ax2, ax3)

plt.grid()
#plt.legend()
#plt.savefig('Costs_vs_target_snr.png', dpi = 600)
#plt.savefig('HIGH_cutoff_vs_target.png', dpi = 600)
#plt.savefig('Triggers_vs_cutoff.png', dpi=600)
plt.savefig('Final_costs.png', dpi=600)
plt.show()
