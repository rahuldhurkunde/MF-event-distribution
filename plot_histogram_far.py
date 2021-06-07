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

no_realizations = 1
#cutoff = [2.0, 2.5, 3.0]
cutoff = [2.0]

filename = "non_hierarchical_matches/only_noise/snr_1_0" 
func.accumulate_triggers(SNR, filename)

og_snr, og_counts = func.non_hierarchical_distribution(ax, SNR, no_realizations)
conv_tol = 0.01
for k in cutoff:
	SNR_hierarchical = []
	print ('First step cutoff is %f' %k)
	for i in range(no_realizations):
		filename_second = "hierarchical_matches/%s/snr_8_%s" % (k, i)
		func.accumulate_triggers(SNR_hierarchical, filename_second)
		#print (i)
	func.hierarchical_distribution(ax, SNR_hierarchical, og_snr, og_counts, no_realizations, k, conv_tol)
	
#func.theoretical_distribution(ax, no_realizations)
plt.yscale('log')
plt.xlabel('SNR')
plt.ylabel('Triggers')
plt.legend(loc=8)
plt.savefig("FAP_hierarchical.png", dpi=600)
plt.show()
exit()

#fig, ax = plt.subplots()
#a = np.loadtxt('first_step_nlouder.txt')
#b = np.loadtxt('second_step_nlouder_data.txt')
#c = np.loadtxt('theoretical_nlouder.txt')

plt.yscale('log')
locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12)
ax.yaxis.set_major_locator(locmaj)
#ax.plot(a[0], a[1] /500, color = 'blue', label = 'full')
#ax.plot(b[0], b[1] * 128, color = 'limegreen', label = 'averaged')
#ax.plot(c[0], c[1] /500, color = 'grey', label = 'theoretical')
#plt.plot( 7.0, 1, 's', label = 'true event')
#plt.plot( 5.32153, 1 , 'o')
plt.title('Event distribution second reconstruction from SNR = 2.5')
plt.xlabel("SNR")
plt.ylabel("No. of events per sec")
plt.grid()
ax.legend(loc = 7)
plt.savefig("FAP_noise_and_theoretical_.png", dpi = 600)
plt.show()


