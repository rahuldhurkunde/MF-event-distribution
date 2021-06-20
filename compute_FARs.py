import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import scipy.stats as stats
import pycbc
from pycbc import events
import event_distr_functions as func
import scipy
import time
import pandas as pd
#from tqdm import tqdm

def write_non_hierarchical(no_realizations):
	SNR = []
	og_snr = []
	og_counts = []
	for i in range(no_realizations):
		print(i)
		#filename = "non_hierarchical_matches/only_noise/snr_1_%s" % (i)
		filename = "non_hierarchical_matches/best_bank/snr_1_%s" % (i)
		func.accumulate_triggers(SNR, filename)
	
	og_snr, og_counts = func.non_hierarchical_distribution(ax, SNR)
	start = time.time()
	save_non_hierarchical(og_snr, og_counts, no_realizations)
	end = time.time()
	print('Max SNR in the series', np.max(og_snr), "Writing time", end-start)
	

def save_non_hierarchical(og_snr, og_counts, no_realizations):
	snr = np.array(og_snr)[:, np.newaxis]
	counts = np.array(og_counts)[:, np.newaxis]
	data = np.concatenate((snr, counts), axis = 1)
	df = pd.DataFrame(data, columns=["snr", "counts"])
	#filename = 'non_hierarchical_matches/ntriggers_non_hierarchical_%s.csv' % no_realizations
	filename = 'non_hierarchical_matches/best_bank/ntriggers_non_hierarchical_%s.csv' % no_realizations
	df.to_csv(filename)

def read_non_hierarchical(filename):
	df = pd.read_csv(filename)
	og_snr = df['snr']
	og_counts = df['counts']
	return og_snr, og_counts

def save_conv_SNRs(first_cutoff, conv_SNR, first_triggers, window, no_realizations):
	cutoff = np.array(first_cutoff)[:, np.newaxis]
	target = np.array(conv_SNR)[:, np.newaxis]
	count = np.array(first_triggers)[:, np.newaxis]
	data = np.concatenate((cutoff, target, count), axis = 1)
	df = pd.DataFrame(data, columns=["cutoff", "target", "count"])
	filename = 'conv_SNRs/%s/conv_SNR_%s_%s.csv' % (window, window, no_realizations)
	df.to_csv(filename, sep=" ")

def save_master_avg_SNR(avg_SNR_master, window, no_realizations):
	df = pd.DataFrame(avg_SNR_master, columns=["avg_SNR"])
	filename = 'hierarchical_matches/sorted_%s/master_avg_SNR_%s.csv' % (window, no_realizations)
	df.to_csv(filename, sep=" ")

def collect_second_triggers(SNR_hierarchical, indices_avg_SNR, window, k):
#	sec_triggers=[]
#	for i in range(k):
#		ind = indices_avg_SNR[i]
#		sec_triggers = np.append(sec_triggers, SNR_hierarchical[ind*window : (ind+1)*window])

	sec_triggers = []
	#for i in range(0, int(k*window)):
	for i in range(k):
		ind = indices_avg_SNR[i]
		for j in range(window):
			sec_triggers.append(SNR_hierarchical[ind*window + j])

	return sec_triggers

def generate_snrs_from_triangular_dist(filename, low, mode, high, window, no_realizations):
	snrs = np.random.triangular(low, mode, high, 300)
	df = pd.DataFrame(snrs, columns=['random_cutoffs'])
	df.to_csv(filename)

def find_cutoff_indices(avg_SNR_master, min_cutoff_ind, window, no_realizations):
	high = avg_SNR_master[min_cutoff_ind]
	low = 2.0
	mode = 4
	filename = 'conv_SNRs/%s/random_cutoffs_%s_%s.csv' % (window, window, no_realizations)
	
	#generate_snrs_from_triangular_dist(filename, low, mode, high, window, no_realizations)
	#
	#df = pd.read_csv(filename)
	#temp_snrs = df['random_cutoffs']
	#temp_snrs = np.append(temp_snrs, high)
	#snrs = np.sort(temp_snrs)

	snrs = np.linspace(low, high, 300)
	
	indices = []
	for i in range(len(snrs)):
		ind = np.abs(avg_SNR_master - snrs[i]).argmin()
		indices.append(ind)
	print(np.max(snrs), "Index", np.max(indices), np.argmax(indices))
	return indices


fig = plt.figure()
ax = fig.add_subplot(111)
start_global = time.time()
window = 8 
no_realizations = 500 
min_cutoff_ind = 50 

write_non_hierarchical(no_realizations)
exit()

start = time.time()
filename = 'non_hierarchical_matches/ntriggers_non_hierarchical_%s.csv' %no_realizations
og_snr, og_counts = read_non_hierarchical(filename)
end = time.time()
print('Max SNR in the series', np.max(og_snr), "Non hierarchical Reading time", end-start)

start = time.time()
avg_SNR = []
SNR_hierarchical = []
for k in range(no_realizations):
	filename_avg = "hierarchical_matches/sorted_%s/avg_%s_%s" % (window, window, k)
	func.accumulate_triggers(avg_SNR, filename_avg)

	filename_second = "hierarchical_matches/sorted_%s/triggers_%s_%s"%(window, window, k)
	func.accumulate_triggers(SNR_hierarchical, filename_second)
end = time.time()
print("Hierarchical reading time", end-start)

start = time.time()
indices_avg_SNR = np.argsort(avg_SNR)
indices_avg_SNR = indices_avg_SNR[::-1]
temp = np.array(avg_SNR)
avg_SNR_master = np.sort(temp)[::-1]
#save_master_avg_SNR(avg_SNR_master, window, no_realizations)
end = time.time()
print("Max avg snr", np.max(avg_SNR_master))
print("Sorting time", end-start)

start = time.time()
cutoff_indices = find_cutoff_indices(avg_SNR_master, min_cutoff_ind, window, no_realizations)
end = time.time()
print("Cutoff indices time", end - start)

start = time.time()
sec_triggers = []
conv_cutoffs = []
conv_SNR = []
conv_counts = []
conv_tol = 0.01
for k in cutoff_indices:
#for k in range(min_cutoff_ind, min_cutoff_ind+100000,100):
	start = time.time()
	start2 = time.time()
	first_step_cutoff = avg_SNR_master[k]

	sec_triggers = collect_second_triggers(SNR_hierarchical, indices_avg_SNR, window, k)
	end2 = time.time()
	print(first_step_cutoff, "Sec triggers collection time", end2 - start2)

	temp_SNR, temp_triggers = func.hierarchical_distribution(ax, sec_triggers, og_snr, og_counts, first_step_cutoff, conv_tol, min_cutoff_ind)
	if (temp_SNR !=0 and temp_triggers !=0):
		conv_cutoffs.append(first_step_cutoff)
		conv_SNR.append(temp_SNR)
		conv_counts.append(k)
	end = time.time()
	print("One cutoff execution time", end-start)
	save_conv_SNRs(conv_cutoffs, conv_SNR, conv_counts, window, no_realizations)


end = time.time()
print("Computation time", end-start)

#
#	conv_cutoffs = []
#	conv_SNR = []
#	conv_counts = []
#	conv_tol = 0.01
#	min_cutoff_ind = 1 
#	start = time.time()
#	for i in range(min_cutoff_ind, len(avg_SNR), 50):
#		cutoff_index = i 
#		first_step_cutoff = avg_SNR[cutoff_index]
#		#print ('First step cutoff is %f' %first_step_cutoff)
#
#		temp_SNR, temp_triggers = func.hierarchical_distribution(ax, SNR_hierarchical[:window*cutoff_index], og_snr, og_counts, no_realizations, first_step_cutoff, conv_tol)
#		if (temp_SNR !=0 and temp_triggers !=0):
#			conv_cutoffs.append(first_step_cutoff)
#			conv_SNR.append(temp_SNR)
#			conv_counts.append(i)
#	end = time.time()
#	save_conv_SNRs(conv_cutoffs, conv_SNR, conv_counts, window, k)
#	print(k, "Max target SNR", np.max(conv_SNR), "max first_step_cutoff", np.max(conv_cutoffs), "max avg SNR", avg_SNR[0])
#
	
end_global = time.time()
print("Total time for ", window , "=", end_global - start_global)

#plt.yscale('log')
#plt.xlabel('SNR')
#plt.ylabel('Triggers per second')
#plt.legend(loc=8)
#plt.savefig("FAP_hierarchical.png", dpi=600)
#plt.show()

