import numpy as np
import pycbc 
from pycbc.waveform import get_fd_waveform
import matplotlib.pyplot as plt
from pycbc import psd, noise, events, types
import scipy.stats as stats
from scipy import interpolate
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000

def compute_match_btw_waveform(our_re, our_im, pycbc_re, pycbc_im, psd):
    product = np.complex(len(pycbc)*2-1, dtype=complex)
    print (len(product))
    #for i in range(len(pycbc_re)):
        #product[i] =  


def plot_relative_error(obs, exp, trigger_list):
	error = []
	for i in range(len(trigger_list)):
		ind = trigger_list[i]
		er = np.abs(1 - np.abs(our_snr[ind]/snr[ind]))
		#er = our_snr[i] - snr[i]
		#er = np.abs(exp_re[i] + 1j*exp_im[i] - obs_re[i] - 1j*obs_im[i])/np.abs(obs_re[i] + 1j*obs_im[i])
		error.append(er)
	print(max(error), np.argmax(error))
	plt.xlabel("Time (sec)")
	plt.ylabel("Relative error")
	plt.title('Relative error for SNR > 4.0')
	#plt.plot(exp.sample_times, error)
	plt.plot(error, label = 'error')
	#plt.plot(obs, '+' ,label='us')
	#plt.plot(exp, 'x', label='pycbc')
	plt.legend()
	plt.savefig("Rel_err.png", dpi = 600)
	plt.show()
    
def waveforms_relative_error(freq, hp_red_re, hp_red_im, hp, psd):
    error = []
    for i in range(len(hp_red_re)):
        #err = 1 -abs(hp_red_re[i].real/hp[i].real)
        err = abs(hp_red_re[i].real - hp[i].real)
        error.append(err)
    plt.plot(error)
    plt.show()
    #hp_re = []
    #hp_im = []
    #for i in range(len(hp)):
    #    hp_re.append(hp[i].real)
    #    hp_im.append(hp[i].imag)

    #new_hp_re = interpolate.interp1d(freq, hp_red_re)
    #new_hp_im = interpolate.interp1d(freq, hp_red_im)
    #new_freq = [(1917+i)/128 for i in range(127360)]
    #plt.plot(hp.sample_frequencies, hp_re, 'x', label='pycbc')
    #plt.plot(new_freq, new_hp_re(new_freq), '+', label='ours interpolated')
    #plt.show()
    #
    #error = []
    #for i in range(len(new_freq)):
    #    freq = (1917 + i)/128.0
    #    err = np.abs(1 - new_hp_re(freq)/hp_re[1920 + i])
    #    #err = np.abs(hp_re[1917+i] + 1j*hp_im[1917+i] - new_hp_re(freq) - 1j*new_hp_im(freq))/np.abs(hp_re[1917+i] + 1j*hp_im[1917+i])
    #    error.append(err)
    #plt.ylabel('relative error')
    ##plt.xlabel('freq index')
    #plt.plot(new_freq, error)
    #plt.xticks(np.arange(min(new_freq), max(new_freq)+1, 60.0))
    #plt.legend()
    #plt.show()
    #print("Values at i = 1707", new_freq[1707], 3627/128.0)
    #print("Waveform values - interpolated us", new_hp_re(28.3359375), "pycbc", hp_re[3627])
    #print("Error values - absolute relative error", error[1707], 'relative error', 1 - new_hp_re(28.3359375)/hp_re[3627])
    #final_interp_hp = np.zeros(len(hp))
    #final_interp_hp[:1917] = 0.0
    #final_interp_hp[1917:129277] = new_hp_re(new_freq)
    #final_interp_hp[129277:] = 0.0
    #final_interp_hp = pycbc.types.frequencyseries.FrequencySeries(final_interp_hp, delta_f = 1.0/128)
    #hp_re = pycbc.types.frequencyseries.FrequencySeries(hp_re, delta_f=1.0/128)
    #snr = pycbc.filter.matchedfilter.match(hp_re, final_interp_hp, psd = psd, low_frequency_cutoff = 15.0, high_frequency_cutoff = 1024.0)
    #snr = np.abs(snr)
    #print(np.max(snr))
    
def read_snr():
    temp = np.loadtxt('noise_match')
    our_snr = []
    for i in range(len(temp)):
        our_snr.append(temp[i])
    return our_snr    
    
def compute_trigger_list(snr, threshold):
	trigger_ind = []
	for i in range(len(snr)):
		if (snr[i] > threshold):
			trigger_ind.append(i)
	return trigger_ind


deltaF = 1.0/128
deltaT = 1.0/2048
fmin = 15.0
fmax = 1024.0

a = np.loadtxt("vec_0")
f = np.loadtxt("freq_reduced_128_new")
interp_template_re = []
interp_template_im = []
freq = []
for i in range(len(a)):
    interp_template_re.append(a[i][1])
    interp_template_im.append(a[i][2])
    freq.append(f[i])

b = np.loadtxt("noise")
noise = []
for i in range(len(b)):
    noise.append(b[i])
noise = types.timeseries.TimeSeries(noise, delta_t = deltaT)
noise_fft = noise.to_frequencyseries()

asd = []
c = np.loadtxt('interp_asd')
for i in range(len(a)):
    asd.append(c[i][1]**2)
psd_us = types.frequencyseries.FrequencySeries(asd, delta_f = deltaF)
psd = psd.analytical.aLIGOZeroDetHighPower(int(2048*128/2+1), delta_f = deltaF, low_freq_cutoff = fmin)
#plot_relative_error(psd_us,psd)
#exit()

mismatch_array = []
for tt in range(0):
    params = np.loadtxt('TB')
    mass1 = params[tt][0]
    mass2 = params[tt][4]
    spin1z = params[tt][3]
    spin2z = params[tt][7]
    hp, hc = get_fd_waveform(approximant="IMRPhenomPv2", mass1=mass1, mass2=mass2, spin1z=spin1z, spin2z = spin2z, f_lower = fmin, f_final = fmax, delta_f = deltaF)
    sigma = pycbc.filter.matchedfilter.sigma(hp, psd = psd, low_frequency_cutoff=fmin, high_frequency_cutoff=fmax)
    hp *= 2.0/sigma

    x = np.loadtxt('reconstructed_wfs/template_%d' %tt)
    our_template = []
    for i in range(len(hp)):
        our_template.append(x[i][1] + 1j*x[i][2])
    our_template = types.frequencyseries.FrequencySeries(our_template, delta_f = deltaF)
    for i in range(len(hp)):
        if ( psd[i] == 0.0):
            our_template[i] = 0.0
        else:    
            our_template[i] = our_template[i]*np.sqrt(psd[i])
    #plt.plot(our_template, '+',label='us')
    #plt.plot(hp, 'x', label='template')
    #plt.legend()
    #plt.show()
    overlap = pycbc.filter.matchedfilter.overlap(hp/2.0, our_template/2.0, psd = psd, low_frequency_cutoff = fmin, high_frequency_cutoff = fmax, normalized = False)
    mismatch = 1 - overlap
    mismatch_array.append(mismatch)
    print(tt, 'Mismatch is', mismatch)
#np.savetxt('mismatch', mismatch_array)    
mismatch_array = np.loadtxt('mismatch_values')

fig = plt.figure()
ax = fig.add_subplot(111)
#plt.yscale('log')
#locmaj = mpl.ticker.LogLocator(base=2,numticks=20)
#ax.yaxis.set_major_locator(locmaj) 
#plt.grid()
plt.ylabel('Mismatch')
plt.xlabel('Tolerance')
plt.hist(mismatch_array)
#plt.boxplot(mismatch_array)
plt.title('Mismatch for tolerance of 1e-6')
#plt.savefig('One_box_plot.png', dpi=600)
plt.show()
#waveforms_relative_error(freq, interp_template_re, interp_template_im, hp, psd)
#waveforms_relative_error(freq, our_template, our_template, hp, psd)
exit()

params = np.loadtxt('TB')
mass1 = params[0][0]
mass2 = params[0][4]
spin1z = params[0][3]
spin2z = params[0][7]
hp, hc = get_fd_waveform(approximant="IMRPhenomPv2", mass1=mass1, mass2=mass2, spin1z=spin1z, spin2z = spin2z, f_lower = fmin, f_final = fmax, delta_f = deltaF)

our_snr = read_snr()
snr = pycbc.filter.matchedfilter.matched_filter(hp, noise, psd = psd, low_frequency_cutoff = fmin, high_frequency_cutoff = fmax)
snr = np.abs(snr)

trigger_ind = compute_trigger_list(snr, 4.0) 
plot_relative_error(our_snr, snr, trigger_ind)
#plt.plot(hp.sample_frequencies, hp, 'x')
#plt.plot(freq, interp_template_re, '+')
#plt.show()

#hp_re =[]
#hp_im =[]
#index =[]
#print(hp[1990], hp[1990].imag)
#f = open('pycbc_wf', 'w')
#for i in range(len(hp)):
#    hp_re.append(hp[i].real)
#    hp_im.append(hp[i].imag)
#    index.append(i)
#    data = np.column_stack((i, hp_re[i], hp_im[i]))
#    np.savetxt(f, data)
#f.close()

