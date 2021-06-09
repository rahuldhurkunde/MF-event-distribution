import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import scipy.stats as stats
import pycbc
from pycbc import events
import event_distr_functions as func

SNR = []
fig = plt.figure()
ax = fig.add_subplot(111)

filename = "non_hierarchical_matches/only_noise/snr_1_0" 
func.accumulate_triggers(SNR, filename)

no_realizations = 1
og_snr, og_counts = func.non_hierarchical_distribution(ax, SNR, no_realizations)

window = 2 
filename = "hierarchical_matches/sorted_%s/avg_%s" % (window, window)
avg_SNR = []
func.accumulate_triggers(avg_SNR, filename)

conv_SNRarray = []
for i in range(25, len(avg_SNR),50):
	cutoff_index = i 
	first_step_cutoff = avg_SNR[cutoff_index]
	#print ('First step cutoff is %f' %first_step_cutoff)
	conv_tol = 0.01
	SNR_hierarchical = []
	filename_second = "hierarchical_matches/sorted_%s/triggers_%s" % (window, window)
	func.accumulate_triggers_hierarchical(SNR_hierarchical, filename_second, cutoff_index, window)

	conv_SNR, conv_triggers = func.hierarchical_distribution(ax, SNR_hierarchical, og_snr, og_counts, no_realizations, first_step_cutoff, conv_tol)
	conv_SNRarray.append([first_step_cutoff, conv_SNR, i])

filen = "conv_SNR_%d.txt" % window
np.savetxt(filen, conv_SNRarray)
	
plt.yscale('log')
plt.xlabel('SNR')
plt.ylabel('Triggers per second')
plt.legend(loc=8)
#plt.savefig("FAP_hierarchical.png", dpi=600)
plt.show()

