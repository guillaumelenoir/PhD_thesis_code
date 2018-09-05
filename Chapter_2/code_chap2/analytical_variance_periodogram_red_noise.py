import numpy as np
import matplotlib.pyplot as plt
import wavepal as wv
from gen_car1_wv33 import gen_car1
import copy

path_to_figure_folder="../figures/"

# Generate a red noise on the time grid of ODP1148
data=np.genfromtxt("data/ODP1148-BF-18O.txt")
t=data[:,0]
sig=np.zeros(1)
sig[0]=2.
tau=np.zeros(1)
tau[0]=20.
mydata=gen_car1(t,1./tau,sig)
mydata=mydata[:,0]
# D=6000.
out1=wv.Wavepal(t, mydata)
out1.check_data()
out1.choose_trend_degree(-1)
out1.trend_vectors()
out1.carma_params(nmcmc=100000)
out1.freq_analysis(n_moments=2,freqstep=0.0001,weighted_WOSA=False,WOSA_segments="max",freqmin=0.003,freqmax=0.022,D=6000.,mywindow=8)
# D=2000.
out2=copy.copy(out1)
out2.freq_analysis(n_moments=2,freqstep=0.0001,weighted_WOSA=False,WOSA_segments="max",freqmin=0.003,freqmax=0.022,D=2000.,mywindow=8)
# D=1000.
out3=copy.copy(out1)
out3.freq_analysis(n_moments=2,freqstep=0.0001,weighted_WOSA=False,WOSA_segments="max",freqmin=0.003,freqmax=0.022,D=1000.,mywindow=8)
# D=600.
out4=copy.copy(out1)
out4.freq_analysis(n_moments=2,freqstep=0.0001,weighted_WOSA=False,WOSA_segments="max",freqmin=0.003,freqmax=0.022,D=600.,mywindow=8)
# Units and labels
freq_label=" ("+"ka"+"${}^{-1}$)"
# figures
freq=out1.freq
var_per1=out1.variance_anal
var_per2=out2.variance_anal
var_per3=out3.variance_anal
var_per4=out4.variance_anal
Q1=out1.nsmooth_vec[0]		# N.B.: with the above option: WOSA_segments="max", Q (=nsmooth_vec) is constant all along the frequencies
Q2=out2.nsmooth_vec[0]
Q3=out3.nsmooth_vec[0]
Q4=out4.nsmooth_vec[0]
plt.semilogy(freq,var_per1,label="Q="+str(Q1))
plt.semilogy(freq,var_per2,label="Q="+str(Q2))
plt.semilogy(freq,var_per3,label="Q="+str(Q3))
plt.semilogy(freq,var_per4,label="Q="+str(Q4))
plt.xlabel("Frequency"+freq_label)
plt.ylabel("Variance")
plt.title("Analytical variance of the WOSA periodogram\n for the background noise")
plt.legend(fancybox=True,fontsize='small')
figname=path_to_figure_folder+"anal_variance_rn.pdf"
plt.savefig(figname)
plt.close()
