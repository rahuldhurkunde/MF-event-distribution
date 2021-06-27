import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import scipy.stats as stats
import pycbc
from pycbc import events
import event_distr_functions as func

def collect_second_triggers(SNR_hierarchical, window, k):
	sec_triggers = []
	for i in range(int(k*window)):
		sec_triggers.append(SNR_hierarchical[i])

	return sec_triggers

def find_cutoff_indices(avg_SNR, min_cutoff_ind):
	high = avg_SNR[min_cutoff_ind]
	low = np.min(avg_SNR)
	snrs = np.linspace(low, high, 300)
	
	indices = []
	for i in range(len(snrs)):
		ind = np.abs(avg_SNR - snrs[i]).argmin()
		indices.append(ind)
	print(np.max(snrs), "Index", np.max(indices), np.argmax(indices))
	return indices

fig = plt.figure()
ax = fig.add_subplot(111)

no_realizations = 1
window = 32 
min_cutoff_ind = 50 

SNR = []
filename = "non_hierarchical_matches/only_noise/snr_1_0" 
func.accumulate_triggers(SNR, filename)
og_snr, og_counts = func.non_hierarchical_distribution(ax, SNR)


avg_SNR = []
filename_avg = "hierarchical_matches/sorted_%s/avg_%s_0" % (window, window)
func.accumulate_triggers(avg_SNR, filename_avg)
print("len avg", len(avg_SNR))

SNR_hierarchical = []
filename_second = "hierarchical_matches/sorted_%s/triggers_%s_0"%(window, window)
func.accumulate_triggers(SNR_hierarchical, filename_second)
print("len hierarchical", len(SNR_hierarchical))


cutoff_indices = find_cutoff_indices(avg_SNR, min_cutoff_ind)

conv_tol = 0.01
conv_SNR = []
conv_cutoffs = []
for k in range(min_cutoff_ind, min_cutoff_ind + 1):
#for k in cutoff_indices:
	first_step_cutoff = avg_SNR[k]
	print ('First step cutoff is %f' %first_step_cutoff)

	sec_triggers = collect_second_triggers(SNR_hierarchical, window, k)
	temp_SNR, temp_triggers = func.hierarchical_distribution(ax, sec_triggers, og_snr, og_counts, first_step_cutoff, conv_tol, min_cutoff_ind)
	if (temp_SNR !=0 and temp_triggers !=0):
		conv_cutoffs.append(first_step_cutoff)
		conv_SNR.append(temp_SNR)
	print("Target SNR is %f" %temp_SNR)

#filename = "conv_SNRs/txtfiles/conv_SNR_%s_1" %window
#np.savetxt(filename, np.column_stack([conv_cutoffs, conv_SNR]))
#plt.plot(conv_cutoffs, conv_SNR, '.')

#func.theoretical_distribution(ax, no_realizations)
plt.axvline(x=temp_SNR, color = 'red')
plt.yscale('log')
plt.xlabel('SNR')
plt.ylabel('Triggers')
plt.legend(loc=8)
plt.savefig("dist_32.png", dpi=600)
plt.show()

