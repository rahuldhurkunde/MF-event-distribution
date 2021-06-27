import matplotlib.pyplot as plt
import scipy
import scipy.stats as stats
import pycbc
from pycbc import events
import numpy as np
import time
import sys
import h5py
import time
import pandas as pd


def stats_indices(SNR, bins):
	indices = []
	for i in range(len(bins)):
		idx = (np.abs(SNR - bins[i])).argmin()
		indices.append(idx)
	return indices

def find_target_SNR(og_snr, red_snr, og_counts, red_counts, first_step_cutoff, tolerance, min_cutoff_ind):
	print("Counts", red_counts[len(red_counts)-1] , og_counts[len(og_counts)-1-min_cutoff_ind], "SNRs", red_snr[len(red_snr)-1] , og_snr[len(og_counts)-1-min_cutoff_ind])
	if( np.abs(1 - red_counts[len(red_counts)-1] / og_counts[len(og_counts)-1] ) > tolerance):
		print("SNRs did not converge for SNR %f" %first_step_cutoff)
		print("Returning zero")
		conv_SNR = 0 
		conv_triggers = 0

	else:
		for i in range(len(red_snr)):
			ind_og = len(og_snr) - 1 - i
			ind_red = len(red_snr) - 1 - i
			print("Ezzz", red_counts[ind_red] , og_counts[ind_og])
			if( np.abs(1 - red_counts[ind_red]/og_counts[ind_og]) > tolerance):
				conv_SNR = og_snr[ind_og + 1]
				conv_triggers = og_counts[ind_og + 1]
				break

	return conv_SNR, conv_triggers


def find_convergent_SNR(og_snr, red_snr, og_counts, red_counts, first_step_cutoff, tolerance):
	if (round(np.max(red_snr),4) > round(np.max(og_snr),4)):
		print("Error for cutoff", first_step_cutoff, "CHECK FLOATING PRECESSION OF MAX. SNRS")
		print(np.max(red_snr), np.max(og_snr))
		conv_SNR = 0 
		conv_triggers = 0
		exit()	
	else:
		x_new = np.linspace(np.min(red_snr), np.max(red_snr), 300)
		f_red = scipy.interpolate.interp1d(red_snr, red_counts)
		f_og = scipy.interpolate.interp1d(og_snr, og_counts, fill_value="extrapolate")	
		#print("YOLO", np.max(red_snr), np.max(og_snr))
		for i in range(len(x_new)):
			if (  ( 1 - f_red(x_new[i])/f_og(x_new[i]) ) < tolerance):
				#plt.axvline(x = x_new[i], color='red')
				conv_SNR = x_new[i]
				conv_triggers = f_red(x_new[i])
				break
			elif (i == (len(x_new)) -1):
				print("SNRs did not converge for SNR %f" %first_step_cutoff)
				print("Returning zero")
				conv_SNR = 0 
				conv_triggers = 0
			
	return conv_SNR, conv_triggers 


def calculate_events_above_threshold(series, values):
	events = np.zeros(len(values))
	for i in range(len(values)):
		events[i] = np.sum(series >= values[i])	
	return events

def accumulate_triggers_hierarchical(SNR, filename, index, window):
	a = np.loadtxt(filename)
	for i in range(0, int(index*window)):
		SNR.append(a[i])

def accumulate_triggers(SNR, filename):
	a = np.loadtxt(filename)
	for i in range(len(a)):
		SNR.append(a[i])
	#SNR.append(temp)
	#return SNR

def non_hierarchical_distribution(ax, SNR):
	SNR = np.array(SNR)
	#thresholds = np.linspace(np.min(SNR), np.max(SNR), 1000)
	snrs  = np.sort(SNR)
	counts = np.arange(len(snrs), 0, -1)
	#ax.plot(snrs, counts,  color='black', label = 'Without hierarchical')
	#plt.xlabel("SNR")
	#plt.ylabel("No. of events per sec")
	#plt.grid()
	#plt.savefig("FAP_only_noise.png", dpi=600)
	#np.savetxt('first_step_nlouder_data.txt', [thresholds, n_louder/128])
	return snrs, counts

def hierarchical_distribution(SNR, og_snr, og_counts, cutoff, conv_tol, min_cutoff_ind):
	SNR = np.array(SNR)
	snrs  = np.sort(SNR)
	counts = np.arange(len(snrs), 0, -1)

	conv_SNR, conv_triggers = find_convergent_SNR(og_snr, snrs, og_counts, counts, cutoff, conv_tol)
	#conv_SNR, conv_triggers = find_target_SNR(og_snr, snrs, og_counts, counts, cutoff, conv_tol, min_cutoff_ind)

	#ax.plot(snrs, counts, label='Hierarchical')
#	plt.xlabel("SNR")
#	plt.ylabel("No. of events per sec")
#	plt.grid()
	#np.savetxt('hierarchical_nlouder_data.txt', [thresholds, n_louder/128])
	return conv_SNR, conv_triggers

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

def save_hdf(snr, filename):
	with h5py.File(filename, 'w') as f:
	    f.create_dataset("snr", data=snr)	
	f.close()
	
def read_hdf(filename):
	hf = h5py.File(filename, 'r')
	temp = hf.get('snr')
	snr = np.array(temp)
	return snr

#def accumulate_hdfs(no_realizations, window, part):
def accumulate_hdfs(no_realizations, window):
	avg_SNR = []
	SNR_hierarchical = []
	#for k in range(part*no_realizations, (part+1)*no_realizations):
	for k in range(no_realizations):
		print(k)
		avg_hdf = "hierarchical_matches/best_bank/hdf_%s/avg_%s_%s.hdf5" %(window, window, k)
		avg = read_hdf(avg_hdf)
		avg_SNR = np.append(avg_SNR, avg)

		triggers_hdf = "hierarchical_matches/best_bank/hdf_%s/triggers_%s_%s.hdf5" %(window, window, k)
		triggers = read_hdf(triggers_hdf)
		SNR_hierarchical = np.append(SNR_hierarchical, triggers)
	return avg_SNR, SNR_hierarchical

