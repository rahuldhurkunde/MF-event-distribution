import sys
import h5py
import time
import pandas as pd
import numpy as np

def read_txt(SNR, filename):
	a = np.loadtxt(filename)
	for i in range(len(a)):
		SNR.append(a[i])	
	return SNR 

def save_hdf(snr, filename):
	with h5py.File(filename, 'w') as f:
	    f.create_dataset("snr", data=snr)	
	f.close()
	
def read_hdf(filename):
	hf = h5py.File(filename, 'r')
	temp = hf.get('snr')
	snr = np.array(temp)
	return snr

def convert_to_hdf(no_realizations, window):
	for k in range(no_realizations):
		print(k)
		avg_SNR = []
		filename_avg = "sorted_%s/avg_%s_%s" % (window, window, k)
		snr = read_txt(avg_SNR, filename_avg)
		avg_hdf = "hdf_%s/avg_%s_%s.hdf5" %(window, window, k)
		save_hdf(snr, avg_hdf)

		SNR_hierarchical = []
		filename_second = "sorted_%s/triggers_%s_%s"%(window, window, k)
		snr_sec = read_txt(SNR_hierarchical, filename_second)
		triggers_hdf = "hdf_%s/triggers_%s_%s.hdf5" %(window, window, k)
		save_hdf(snr_sec, triggers_hdf)

def accumulate_hdfs(no_realizations, window):
	avg_SNR = []
	SNR_hierarchical = []
	for k in range(no_realizations):
		print(k)
		avg_hdf = "hdf_%s/avg_%s_%s.hdf5" %(window, window, k)
		avg = read_hdf(avg_hdf)
		avg_SNR = np.append(avg_SNR, avg)

		triggers_hdf = "hdf_%s/triggers_%s_%s.hdf5" %(window, window, k)
		triggers = read_hdf(triggers_hdf)
		SNR_hierarchical = np.append(SNR_hierarchical, triggers)
	return avg_SNR, SNR_hierarchical

window = 32 
no_realizations = 5000 

convert_to_hdf(no_realizations, window)
#avg_SNR, SNR_hierarchical = accumulate_hdfs(no_realizations, window)

