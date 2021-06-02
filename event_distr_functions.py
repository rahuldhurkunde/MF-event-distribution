import matplotlib.pyplot as plt
import scipy.stats as stats
import pycbc
from pycbc import events
import numpy as np

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def calculate_events_above_threshold(series, values):
	events = np.zeros(len(values))
	for i in range(len(values)):
		events[i] = np.sum(series >= values[i])	
	return events

def accumulate_triggers(SNR, filename):
	a = np.loadtxt(filename)
	#temp = np.zeros(len(a))
	for i in range(len(a)):
		SNR.append(a[i])
	#SNR.append(temp)
	#return SNR

def non_hierarchical_distribution(ax, SNR, no_realizations):
	SNR = np.array(SNR)
	#thresholds = np.linspace(np.min(SNR), np.max(SNR), 1000)
	snrs  = np.sort(SNR)
	counts = np.arange(len(snrs), 0, -1)
	ax.plot(snrs, counts, 'x', color='black', label = 'Without hierarchical')
	ax.set_yscale('log')
	plt.xlabel("SNR")
	plt.ylabel("No. of events per sec")
	plt.grid()
	#plt.savefig("FAP_only_noise.png", dpi=600)
	#np.savetxt('first_step_nlouder_data.txt', [thresholds, n_louder/128])

def hierarchical_distribution(ax, SNR, no_realizations, cutoff):
	SNR = np.array(SNR)
	snrs  = np.sort(SNR)
	counts = np.arange(len(snrs), 0, -1)
	ax.plot(snrs, counts, '+', label='Hierarchical %0.1f' %cutoff)
	ax.set_yscale('log')
	plt.xlabel("SNR")
	plt.ylabel("No. of events per sec")
	plt.grid()
	#np.savetxt('hierarchical_nlouder_data.txt', [thresholds, n_louder/128])

def theoretical_distribution(ax, no_realizations):
	#x1 = np.random.normal(0,1,1000000)
	#x2 = np.random.normal(0,1,1000000)
	#SNR = np.sqrt(x1**2 + x2**2)
	x = np.random.chisquare(2, 262144 * no_realizations)
	SNR = np.sqrt(x)
	#print('Full theoretical -', sum(i > 2 for i in SNR)/262144, sum(i > 2 for i in SNR))
	#plt.hist(SNR, bins=50)
	thresholds = np.linspace(np.min(SNR), np.max(SNR), 10000)
	dec = np.ones(len(SNR))
	n_louder = events.coinc.calculate_n_louder(SNR, thresholds, dec, skip_background=True)	
	ax.plot(thresholds, n_louder/128/no_realizations, color='blue', label = 'Theoretical') 
	#plt.savefig("FAP_theoretical.png", dpi=600)
	#np.savetxt('theoretical_nlouder.txt', [thresholds, n_louder/128])

def theoretical_distribution_reduced(ax, no_realizations, window):
	x = np.random.chisquare(2, 262144 * no_realizations)
	SNR = np.sqrt(x)
	bins = 262144/window
	SNR_red = np.zeros(int(bins))
	for i in range(int(bins)-1):
		temp = 0.0
		for k in range(window):
			temp += SNR[i*window + k]
		temp /= window
		SNR_red[i] = temp
	print("Red", sum(i > 2.5 for i in SNR_red)/262144*window, sum(i > 2.5 for i in SNR_red))
	#plt.hist(SNR, bins=50)
	thresholds = np.linspace(np.min(SNR_red), np.max(SNR_red), 10000)
	dec = np.ones(len(SNR_red))
	n_louder = events.coinc.calculate_n_louder(SNR_red, thresholds, dec, skip_background=True)	
	ax.plot(thresholds, n_louder/128, color='blue', label = 'Theoretical_red') 
	#plt.savefig("FAP_theoretical.png", dpi=600)
	#np.savetxt('theoretical_nlouder.txt', [thresholds, n_louder/128])
