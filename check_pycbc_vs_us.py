import numpy as np
import pycbc 
from pycbc.waveform import get_fd_waveform
import matplotlib.pyplot as plt
from pycbc import psd, noise, events, types
import scipy.stats as stats
from scipy import interpolate
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000

def plot_relative_error(obs, exp):
    error = []
    for i in range(len(obs)):
        #er = obs[i] - exp[i]
        er = np.abs(1 - obs[i]/exp[i])
        error.append(er)
    print(max(error), np.argmax(error))
    #plt.plot(snr.sample_times, np.abs(snr), '.')
    #plt.plot(snr.sample_times, our_snr, '.')
    plt.xlabel("Time (sec)")
    plt.ylabel("Relative error")
    plt.title('Relative error of SNR btw pycbc and our scheme')
    #plt.plot(exp.sample_times, error)
    #plt.plot(error, label = 'error')
    plt.plot(obs, '+' ,label='us')
    plt.plot(exp, 'x', label='pycbc')
    plt.legend()
    plt.savefig("Abs_err.png", dpi = 600)
    plt.show()
    
def waveforms_relative_error(freq, hp_red_re, hp_red_im, hp, psd):
    hp_re = []
    hp_im = []
    for i in range(len(hp)):
        hp_re.append(hp[i].real)
        hp_im.append(hp[i].imag)

    new_hp_re = interpolate.interp1d(freq, hp_red_re)
    new_hp_im = interpolate.interp1d(freq, hp_red_im)
    new_freq = [(1917+i)/128 for i in range(127360)]
    #plt.plot(hp.sample_frequencies, hp_re, 'x', label='pycbc')
    #plt.plot(new_freq, new_hp_re(new_freq), '+', label='ours interpolated')
    #plt.show()
    
    error = []
    for i in range(len(new_freq)):
        freq = (1917 + i)/128.0
        #err = np.abs(1 - new_hp_re(freq)/hp_re[1920 + i])
        err = np.abs(hp_re[1917+i] + 1j*hp_im[1917+i] - new_hp_re(freq) - 1j*new_hp_im(freq))/np.abs(hp_re[1917+i] + 1j*hp_im[1917+i])
        error.append(err)
    plt.ylabel('relative error')
    #plt.xlabel('freq index')
    plt.plot(new_freq, error)
    plt.xticks(np.arange(min(new_freq), max(new_freq)+1, 60.0))
    plt.legend()
    plt.show()
    print("Values at i = 1707", new_freq[1707], 3627/128.0)
    print("Waveform values - interpolated us", new_hp_re(28.3359375), "pycbc", hp_re[3627])
    print("Error values - absolute relative error", error[1707], 'relative error', 1 - new_hp_re(28.3359375)/hp_re[3627])
    final_interp_hp = np.zeros(len(hp))
    final_interp_hp[:1917] = 0.0
    final_interp_hp[1917:129277] = new_hp_re(new_freq)
    final_interp_hp[129277:] = 0.0
    final_interp_hp = pycbc.types.frequencyseries.FrequencySeries(final_interp_hp, delta_f = 1.0/128)
    hp_re = pycbc.types.frequencyseries.FrequencySeries(hp_re, delta_f=1.0/128)
    snr = pycbc.filter.matchedfilter.match(hp_re, final_interp_hp, psd = psd, low_frequency_cutoff = 15.0, high_frequency_cutoff = 1024.0)
    snr = np.abs(snr)
    print(np.max(snr))
    

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
psd = psd.analytical.aLIGOZeroDetHighPower(len(noise_fft), delta_f = deltaF, low_freq_cutoff = fmin)
#plot_relative_error(psd_us,psd)
#exit()

mass1 = 6.4156613
mass2 = 2.1420077
spin1z = 0.70474228
spin2z = -0.77516922
hp, hc = get_fd_waveform(approximant="IMRPhenomPv2", mass1=mass1, mass2=mass2, spin1z=spin1z, spin2z = spin2z, f_lower = fmin, f_final = fmax, delta_f = deltaF)
#waveforms_relative_error(freq, interp_template_re, interp_template_im, hp, psd)
#exit()

snr = pycbc.filter.matchedfilter.matched_filter(hp, noise, psd = psd, low_frequency_cutoff = fmin, high_frequency_cutoff = fmax)
snr = np.abs(snr)
#snr = 2.0 * snr

temp = np.loadtxt('noise_match')
our_snr = []
for i in range(len(temp)):
    our_snr.append(temp[i])


plot_relative_error(our_snr, snr)
#plt.plot(hp.sample_frequencies, hp, 'x')
#plt.plot(freq, interp_template_re, '+')
#plt.show()

