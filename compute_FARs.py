import sys
import h5py
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
from guppy import hpy
h = hpy()
plt.rcParams.update({
    "text.usetex": True})

#def write_non_hierarchical(ax, no_realizations, bank, non_hierarchical_file, part):
def write_non_hierarchical(ax, no_realizations, bank, non_hierarchical_file):
	SNR = []
	og_snr = []
	og_counts = []
	#for i in range(part*no_realizations, (part+1)*no_realizations):
	for i in range(no_realizations):
		print(i)
		filename = "non_hierarchical_matches/%s/snr_1_%s" % (bank, i)
		func.accumulate_triggers(SNR, filename)
	
	og_snr, og_counts = func.non_hierarchical_distribution(ax, SNR)
	start = time.time()
	save_non_hierarchical(og_snr, og_counts, no_realizations, non_hierarchical_file)
	end = time.time()
	print('Max SNR in the series', np.max(og_snr), "Writing time", end-start)
	

def save_non_hierarchical(og_snr, og_counts, no_realizations, non_hierarchical_file):
	#snr = np.array(og_snr)[:, np.newaxis]
	#counts = np.array(og_counts)[:, np.newaxis]
	#data = np.concatenate((snr, counts), axis = 1)
	#df = pd.DataFrame(data, columns=["snr", "counts"])
	#df.to_csv(non_hierarchical_file)
	with h5py.File(non_hierarchical_file, 'w') as f:
	    f.create_dataset("snr", data=og_snr)	
	    f.create_dataset("counts", data=og_counts)
	f.close()

def read_non_hierarchical(filename):
	#yolo = 'non_hierarchical_matches/best_bank/ntriggers_non_hierarchical_16_5.csv'
	#df = pd.read_csv(yolo)
	#og_snr = df['snr']
	#og_counts = df['counts']
	#print ('csv', og_snr.shape)
	hf = h5py.File(filename, 'r')
	og_snr = hf['snr']
	og_counts = hf['counts']
	return og_snr, og_counts

def save_conv_SNRs(first_cutoff, conv_SNR, first_triggers, window, no_realizations, conv_filename):
	cutoff = np.array(first_cutoff)[:, np.newaxis]
	target = np.array(conv_SNR)[:, np.newaxis]
	count = np.array(first_triggers)[:, np.newaxis]
	data = np.concatenate((cutoff, target, count), axis = 1)
	df = pd.DataFrame(data, columns=["cutoff", "target", "count"])
	df.to_csv(conv_filename, sep=" ")

def save_master_avg_SNR(avg_SNR_master, window, no_realizations):
	df = pd.DataFrame(avg_SNR_master, columns=["avg_SNR"])
	filename = 'hierarchical_matches/sorted_%s/master_avg_SNR_%s.csv' % (window, no_realizations)
	df.to_csv(filename, sep=" ")

def collect_second_triggers(SNR_hierarchical, indices_avg_SNR, window, k):
	sec_triggers = []
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
	low = 1.5 
	mode = 4
	#filename = 'conv_SNRs/%s/random_cutoffs_%s_%s.csv' % (window, window, no_realizations)
	
	#generate_snrs_from_triangular_dist(filename, low, mode, high, window, no_realizations)
	#
	#df = pd.read_csv(filename)
	#temp_snrs = df['random_cutoffs']
	#temp_snrs = np.append(temp_snrs, high)
	#snrs = np.sort(temp_snrs)

	snrs = np.linspace(low, high, 200)
	
	indices = []
	for i in range(len(snrs)):
		ind = np.abs(avg_SNR_master - snrs[i]).argmin()
		indices.append(ind)
	return indices

fig = plt.figure()
ax = fig.add_subplot(111)

start_global = time.time()
window = 2 
no_realizations = 5000
part = 2 
min_cutoff_ind = 1 
print("\t \t \t Window =", window, "\t No. of files", no_realizations)

##bank = 'only_noise'
#print("Writing NON-hierarchical file")
##non_hierarchical_file = 'non_hierarchical_matches/ntriggers_non_hierarchical_%s.csv' % no_realizations
#bank = 'best_bank'
##non_hierarchical_file = 'non_hierarchical_matches/best_bank/ntriggers_non_hierarchical_%s_%s_%s.hdf5' % (window, no_realizations, part)
#non_hierarchical_file = 'non_hierarchical_matches/best_bank/ntriggers_non_hierarchical_%s.hdf5' % (no_realizations)
#write_non_hierarchical(ax, no_realizations, bank, non_hierarchical_file)
#exit()

start = time.time()
#filename = 'non_hierarchical_matches/ntriggers_non_hierarchical_%s.csv' %no_realizations
#filename = 'non_hierarchical_matches/best_bank/ntriggers_non_hierarchical_%s_%s.csv' %(window, no_realizations)
#filename = 'non_hierarchical_matches/best_bank/ntriggers_non_hierarchical_%s_%s_%s.hdf5' %(window, no_realizations, part)
filename = 'non_hierarchical_matches/best_bank/ntriggers_non_hierarchical_%s.hdf5' %(no_realizations)
print("Reading %s" %filename)
og_snr, og_counts = read_non_hierarchical(filename)
end = time.time()
print('Max SNR in the series', np.max(og_snr), "Non hierarchical Reading time", end-start, "Size NON-hierarchical", sys.getsizeof(og_snr))

start = time.time()
#avg_SNR, SNR_hierarchical = func.accumulate_hdfs(no_realizations, window, part)
avg_SNR, SNR_hierarchical = func.accumulate_hdfs(no_realizations, window)
end = time.time()
print("Hierarchical reading time", end-start, "Size of Avg_SNR", sys.getsizeof(avg_SNR), "Size of SNR_hierarchical", sys.getsizeof(SNR_hierarchical))



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
cutoff_indices = cutoff_indices[::-1]
end = time.time()
print("Cutoff indices time", end - start)

start = time.time()
conv_cutoffs = []
conv_SNR = []
conv_counts = []
conv_tol = 0.01
for k in cutoff_indices:
	start = time.time()
	start2 = time.time()
	first_step_cutoff = avg_SNR_master[k]

	sec_triggers = collect_second_triggers(SNR_hierarchical, indices_avg_SNR, window, k)
	end2 = time.time()
	print('Cutoff', first_step_cutoff, "Sec triggers collection time", end2 - start2, "Size of second triggers", sys.getsizeof(sec_triggers))

	temp_SNR, temp_triggers = func.hierarchical_distribution(sec_triggers, og_snr, og_counts, first_step_cutoff, conv_tol, min_cutoff_ind)
	if (temp_SNR !=0 and temp_triggers !=0):
		conv_cutoffs.append(first_step_cutoff)
		conv_SNR.append(temp_SNR)
		conv_counts.append(k)
	end = time.time()
	print("One cutoff execution time", end-start, "Target_snr", temp_SNR)
	print('==========================')

	#conv_filename = 'conv_SNRs/%s/conv_SNR_%s_%s.csv' % (window, window, no_realizations)
	conv_filename = 'conv_SNRs/best_bank/%s/conv_SNR_%s_%s.csv' % (window, window, no_realizations)
	save_conv_SNRs(conv_cutoffs, conv_SNR, conv_counts, window, no_realizations, conv_filename)
	del sec_triggers

end = time.time()
print("Computation time", end-start)
	
end_global = time.time()
print("Total time for ", window , "=", end_global - start_global)

print(h.heap())












#avg_SNR = []
#SNR_hierarchical = []
#for k in range(no_realizations):
#	#filename_avg = "hierarchical_matches/sorted_%s/avg_%s_%s" % (window, window, k)
#	filename_avg = "hierarchical_matches/best_bank/sorted_%s/avg_%s_%s" % (window, window, k)
#	func.accumulate_triggers(avg_SNR, filename_avg)
#
#	#filename_second = "hierarchical_matches/sorted_%s/triggers_%s_%s"%(window, window, k)
#	filename_second = "hierarchical_matches/best_bank/sorted_%s/triggers_%s_%s"%(window, window, k)
#	func.accumulate_triggers(SNR_hierarchical, filename_second)

