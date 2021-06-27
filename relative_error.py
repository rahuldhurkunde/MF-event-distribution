import numpy as np
import pycbc 
from pycbc.waveform import get_fd_waveform
import matplotlib.pyplot as plt
from pycbc import psd, noise, events, types
import scipy.stats as stats
from scipy import interpolate
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000
import time
import csv
import pandas as pd
import array



class SNR:
	def __init__(self, pycbc, our, index):
		self.pycbc = pycbc
		self.our = our
		self.index = index


def plot_relative_error(snr, our_snr,  trigger_list):
	error = []
	for i in range(len(trigger_list)):
		ind = trigger_list[i]
		er = np.abs(1 - np.abs(our_snr[ind]/snr[ind]))
		error.append(er)
	#print(max(error), np.argmax(error))
	#plt.xlabel("Time (sec)")
	#plt.ylabel("Relative error")
	#plt.title('Relative error for SNR > 4.0')
	#plt.plot(error, label = 'error')
	#plt.legend()
	#plt.savefig("Rel_err.png", dpi = 600)
	#plt.show()
	return max(error)

def generate_wf(fmin, fmax, deltaF):
	params = np.loadtxt('best_bank')
	mass1 = params[0][0]
	mass2 = params[0][4]
	spin1z = params[0][3]
	spin2z = params[0][7]
	hp, hc = get_fd_waveform(approximant="IMRPhenomPv2", mass1=mass1, mass2=mass2, spin1z=spin1z, spin2z = spin2z, f_lower = fmin, f_final = fmax, delta_f = deltaF)
	return hp, hc

def read_snr(i):
	filename = 'non_hierarchical_matches/best_bank/snr_1_%s' %(i)
	temp = np.loadtxt(filename)
	our_snr = []
	for i in range(len(temp)):
		our_snr.append(temp[i])
	return our_snr    

def read_noise(i, deltaT):
	filename = '/work/rahul.dhurkunde/CudaPCA/best_bank/Injections/noise_%s' %(i+1)
	noise = np.loadtxt(filename)
	noise = types.timeseries.TimeSeries(noise, delta_t = deltaT)
	noise_fft = noise.to_frequencyseries()
	return noise	

def compute_trigger_list_initial(snr, threshold):
	trigger_ind = []
	for i in range(len(snr)):
		if (snr[i] > threshold):
			trigger_ind.append(i)
	return trigger_ind

def compute_trigger_list_next(snr, trigger_ind, threshold):
	new_trigger = []
	for i in range(len(trigger_ind)):
		scan_index = trigger_ind[i]
		if (snr[scan_index] > threshold):
			new_trigger.append(scan_index)
	return new_trigger

deltaF = 1.0/128
deltaT = 1.0/2048
fmin = 15.0
fmax = 1024.0
hp,hc = generate_wf(fmin, fmax, deltaF)
psd = psd.analytical.aLIGOZeroDetHighPower(int(2048*128/2+1), delta_f = deltaF, low_freq_cutoff = fmin)

start = time.time()
cutoff = np.linspace(3.0, 5.4, 11)
#cutoff = [3.0,4.0]

no_of_realizations = 100
rel_err = []

#temp_max_err = np.empty((no_of_realizations, len(cutoff)), float)
#PCA = []
#PYCBC = []
#for i in range(no_of_realizations):
#	print(i)
#	our_snr = read_snr(i)
#	noise = read_noise(i, deltaT)
#	snr = pycbc.filter.matchedfilter.matched_filter(hp, noise, psd = psd, low_frequency_cutoff = fmin, high_frequency_cutoff = fmax)
#	snr = np.abs(snr)
#	PCA = np.append(PCA, our_snr)
#	PYCBC = np.append(PYCBC, snr)
#
#PCA = np.array(PCA)[:, np.newaxis]
#PYCBC = np.array(PYCBC)[:, np.newaxis]
#SNRs = np.concatenate((PCA, PYCBC), axis = 1)
#df = pd.DataFrame(SNRs, columns=["ours", "pycbc"])
#df.to_csv('rel_err_files/combined_SNRs.csv')

#start = time.time()
#df = pd.read_csv('SNRs.csv')
#
#snr = df['pycbc']
#our_snr = df['ours']
#end = time.time()
#print("Reading time", end-start)


#start = time.time()
#snr_list = []
#for i in range(len(snr)):
#	temp = SNR(snr[i], our_snr[i], index)
#	snr_list.append(temp)	
#end = time.time()
#print("Object initialization", end-start)

#trigger_ind = []
#for k in range(len(cutoff)):
#	if (k == 0):
#		trigger_ind = compute_trigger_list_initial(snr, cutoff[k])
#		max_error = plot_relative_error(snr, our_snr, trigger_ind)
#		print(k, len(trigger_ind), max_error)
#		rel_err.append(max_error)
#	else:
#		trigger_ind = compute_trigger_list_next(snr, trigger_ind, cutoff[k])
#		max_error = plot_relative_error(our_snr, snr, trigger_ind)
#		print(k, len(trigger_ind), max_error)
#		rel_err.append(max_error)
#
#print("Error_list iz", rel_err)
#np.savetxt("Rel_err", rel_err)

#exit()
err = np.loadtxt('Rel_err')
plt.plot(cutoff, err)
plt.xlabel("Threshold")
plt.ylabel("Relative error")
plt.grid()
plt.savefig("Rel_err.png", dpi=600)
plt.show()
end = time.time()
print("Cutoff time = ", end-start)



