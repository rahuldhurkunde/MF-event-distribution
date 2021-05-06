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
	temp = np.zeros(len(a))
	for i in range(len(a)):
		SNR.append(a[i][1])
	#SNR.append(temp)
	#return SNR

def first_step_distribution(ax, SNR, no_realizations):
	SNR = np.array(SNR)
	thresholds = np.linspace(np.min(SNR), np.max(SNR), 1000)
	dec = np.ones(len(SNR))
	n_louder = events.coinc.calculate_n_louder(SNR, thresholds, dec, skip_background=True)	
	ax.plot(thresholds, n_louder/128, color='black', label = 'noise background')
	ax.set_yscale('log')
	plt.xlabel("SNR")
	plt.ylabel("No. of events per sec")
	plt.grid()
	#plt.savefig("FAP_only_noise.png", dpi=600)
	np.savetxt('first_step_nlouder.txt', [thresholds, n_louder])

def second_step_distribution(ax, SNR, no_realizations):
	SNR = np.array(SNR) 
	thresholds = np.linspace(np.min(SNR), np.max(SNR), 100)
	dec = np.ones(len(SNR))
	n_louder = events.coinc.calculate_n_louder(SNR, thresholds, dec, skip_background=True)	
	ax.plot(thresholds, n_louder/128, color='black', label = 'reduced-noise background')
	ax.set_yscale('log')
	plt.xlabel("SNR")
	plt.ylabel("No. of events per sec")
	plt.grid()
	np.savetxt('second_step_nlouder.txt', [thresholds, n_louder])

def theoretical_distribution(ax, no_realizations):
	#x1 = np.random.normal(0,1,1000000)
	#x2 = np.random.normal(0,1,1000000)
	#SNR = np.sqrt(x1**2 + x2**2)
	x = np.random.chisquare(2, 262144 * no_realizations)
	SNR = np.sqrt(x)
	#plt.hist(SNR, bins=50)
	thresholds = np.linspace(np.min(SNR), np.max(SNR), 10000)
	dec = np.ones(len(SNR))
	n_louder = events.coinc.calculate_n_louder(SNR, thresholds, dec, skip_background=True)	
	ax.plot(thresholds, n_louder/128, color='blue', label = 'Theoretical') 
	#plt.savefig("FAP_theoretical.png", dpi=600)
	np.savetxt('theoretical_nlouder.txt', [thresholds, n_louder])
