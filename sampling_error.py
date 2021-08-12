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
plt.rcParams["text.usetex"] =True

def generate_pycbc_template(params, psd, fmin, fmax, deltaF, tt):
    mass1 = params[tt][0]
    mass2 = params[tt][4]
    spin1z = params[tt][3]
    spin2z = params[tt][7]
    hp, hc = get_fd_waveform(approximant="IMRPhenomPv2", mass1=mass1, mass2=mass2, spin1z=spin1z, spin2z = spin2z, f_lower = fmin, f_final = fmax, delta_f = deltaF)
    sigma = pycbc.filter.matchedfilter.sigma(hp, psd = psd, low_frequency_cutoff=fmin, high_frequency_cutoff=fmax)
    hp *= 1.0/sigma
    #plt.plot(hp.sample_frequencies, hp, 'x', label='pycbc')
    return hp

def read_our_template(vec_file):
	x = np.loadtxt(vec_file)
	freq = []
	vec_re = []
	vec_im = []
	for i in range(len(x)):
		freq.append(x[i][0])
		vec_re.append(x[i][1])
		vec_im.append(x[i][2])
	return freq, vec_re, vec_im

def interpolate_vec(psd, freq, vec_re, vec_im):
	re = interpolate.interp1d(freq, vec_re)
	im = interpolate.interp1d(freq, vec_im)
	new_freq = [(1918+i)/128.0 for i in range(127360)]
    
	interp_re = np.zeros(len(psd))
	interp_im = np.zeros(len(psd))
	dft_freq = []
	k = 0
	for i in range(len(psd)):
		dft_freq.append(i/128.0)
		if (i >= 1918 and i < 127360):
			interp_re[i] = re(new_freq[k])
			interp_im[i] = im(new_freq[k]) 
			k += 1
		else:
			interp_re[i] = 0.0
			interp_im[i] = 0.0
	return dft_freq, interp_re, interp_im


def unwhiten_vectors(dft_freq, interp_re, interp_im, psd):
	for i in range(len(psd)):
		if ( psd[i] == 0 ):
			interp_re[i] = 0.0
			interp_im[i] = 0.0
		else:
			interp_re[i] = interp_re[i]*np.sqrt(psd[i])
			interp_im[i] = interp_im[i]*np.sqrt(psd[i])
	our_template = interp_re + 1j*interp_im
	return our_template

def plot_mismatch(params, filename):
	a = np.loadtxt(filename)
	mismatch=[]
	total_mass=[]
	for i in range(len(a)):
		mismatch.append(a[i])
		total_mass.append(params[i][0]+params[i][4])

	unique_mismatch, ind = np.unique(mismatch, return_index = True)
	unique_mass = []
	for i in range(len(ind)):
		unique_mass.append(total_mass[ind[i]])
	plt.plot( unique_mass, unique_mismatch, '.')
	plt.grid()
	plt.xlabel('$M$(total mass)')
	plt.ylabel('mismatch')
	plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
	plt.savefig('Sampling_mismatch.png', dpi=600)
	plt.show()

fmin = 15.0
fmax = 1024.0
deltaF = 1/128.0
psd = psd.analytical.aLIGOZeroDetHighPower(int(2048*128/2+1), delta_f = deltaF, low_freq_cutoff = fmin)
mismatch_array = []
params = np.loadtxt('best_bank')

for tt in range(0):
	hp = generate_pycbc_template(params, psd, fmin, fmax, deltaF, tt)

	vec_file = '../CudaPCA/best_bank/vectors/vec_%d' %tt
	freq, vec_re, vec_im = read_our_template(vec_file)
	dft_freq, interp_re, interp_im = interpolate_vec(psd, freq, vec_re, vec_im)
	our_template = unwhiten_vectors(dft_freq, interp_re, interp_im, psd)
	#plt.plot(dft_freq, interp_re, '+',label = 'before whitening')
	#plt.plot(dft_freq, our_template, '+', label = 'ours')
	#plt.legend()
	#plt.show()	

	our_template = pycbc.types.frequencyseries.FrequencySeries(our_template, delta_f = 1.0/128)

	overlap = pycbc.filter.matchedfilter.overlap(hp, our_template/2.0, psd = psd, low_frequency_cutoff = fmin, high_frequency_cutoff = fmax, normalized = True)
	mismatch = 1 - overlap
	mismatch_array.append(mismatch)
	print(tt, 'Mismatch is', mismatch)
	data = np.transpose(mismatch_array)
#np.savetxt('sampling_mismatch', data)    
filename = 'sampling_mismatch'
plot_mismatch(params, filename)
