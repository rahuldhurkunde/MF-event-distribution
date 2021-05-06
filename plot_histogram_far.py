import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pycbc
from pycbc import events
import event_distr_functions as func

SNR = []
SNR_2nd_step = []

no_realizations = 1 
for i in range(no_realizations):
	filename = "only_noise_matches/noise_snr_1_%s" % i
	func.accumulate_triggers(SNR, filename)
	filename_2nd = "only_noise_matches/red_noise_snr_8_%s" % i
	func.accumulate_triggers(SNR_2nd_step, filename_2nd)
	print (i)

fig = plt.figure()
ax = fig.add_subplot(111)

func.first_step_distribution(ax, SNR, no_realizations)
func.second_step_distribution(ax, SNR_2nd_step, no_realizations)
func.theoretical_distribution(ax, no_realizations)

a = np.loadtxt('first_step_nlouder.txt')
b = np.loadtxt('second_step_nlouder.txt')
c = np.loadtxt('theoretical_nlouder.txt')
plt.legend(loc = 7)
plt.savefig("FAP_noise_and_red.png", dpi=600)
#plt.show()
#ax.set_yscale('log')
#plt.xlabel("SNR")
#plt.ylabel("False alarm rate")
#plt.legend(loc = 7)
#plt.grid()
#plt.savefig("FAP_only_noise.png", dpi=600)
#plt.show()

#a_n = np.loadtxt('only_noise_matches/match_og')
#b_n = np.loadtxt('only_noise_matches/match_16')
#SNR_n = np.zeros(len(a_n))
#SNR_2_n = np.zeros(len(b_n))


#x1 = np.random.normal(0, 1.0, len(SNR))
#x2 = np.random.normal(0, 1.0, len(SNR))
#new_stat = np.sqrt(x1**2 + x2**2)
#x_g = np.histogram(new_stat, bins = 100, density = True)[1]
#dec_g = np.ones(len(SNR))
#n_g = events.coinc.calculate_n_louder(new_stat, x_g, dec_g, skip_background=True)
#ax.plot(x_g, n_g)
##plt.hist(SNR, bins=100, density=True)
###plt.plot(x, stats.chi2.pdf(x,df=2), color='r', lw=2)
##plt.show()
##exit()

#for i in range(len(b)):
#	SNR_2[i] = b[i][1]
#	SNR_2_n[i] = b_n[i][1]
#n_2, x_2 = np.histogram(SNR_2, bins = 100, density=True)
#n_2_n, x_2_n = np.histogram(SNR_2_n, bins = 100, density=True)

#yolo = calculate_events_above_threshold(SNR, x)


#dec_2 = np.ones(len(SNR_2))
#n_louder_2 = events.coinc.calculate_n_louder(SNR_2, x_2, dec_2, skip_background=True)
#n_louder_2_n = events.coinc.calculate_n_louder(SNR_2_n, x_2_n, dec_2, skip_background=True)
#ax.plot(x_2, n_louder_2, color='limegreen', label = 'Noise Background')
#ax.plot(x_2_n, n_louder_2_n, '--', color='limegreen')
##ax.plot(x, yolo), '.')
#ax.plot(2.635,1, 's', color = 'limegreen')


