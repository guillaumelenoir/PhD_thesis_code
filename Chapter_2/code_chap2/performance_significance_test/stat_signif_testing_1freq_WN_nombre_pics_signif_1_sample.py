
# signal = white noise + 1 period

# computes the average number of false "peak" detection PER FREQUENCY_UNIT. It thus relies on a peak detection algorithm.
# Illustration of the periodogram of 1 sample

import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from scipy.stats import chi2 as chi2distr
from periodogram import periodogram
from percentile_n_moments import percentile_n_moments
from dt_av_central import dt_av_central
from dt_av import dt_av
import copy
from scipy.interpolate import interp1d
import random
import time
from tapered_periodogram import tapered_periodogram
from freq_analysis_prelims import freq_analysis_prelims
from LS_WOSA import LS_WOSA
import numpy.linalg as la
import peakutils
from periodogram_grid_fixed import periodogram_grid_fixed
from WOSA_periodogram_grid_fixed import WOSA_periodogram_grid_fixed
from tapered_periodogram_grid_fixed import tapered_periodogram_grid_fixed


# Utilisation de peakutils: threshold a entrer = (true_threshold-min_fun)/(max_fun-min_fun)


#params
proba=0.9
A=1.
t=np.arange(1000.)
nt_del=900
nfreq=2000	# Nombre de points de frequence sur l'intervalle [0.,1./2./dt]
beta=0.75	# for the WOSA
D_WOSA=200.		# for the WOSA
taper_number=8	# for the tapered per. and the WOSA
taper_number2=10	# for the tapered per. and the WOSA
coverage=90.	# for the WOSA
freq_signal=0.02

#### WOSA tapered LS periodogram - 2-moment approx for the signif. levels
nt=t.size
dt=t[1]-t[0]
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
zero_array=np.zeros((nt,2))
zero_array2=np.zeros((nt-nt_del,2))
fmin=1./D_WOSA
fmax=1./2./dt
nfreq1=int(np.floor((fmax-fmin)/(1./2./dt)*float(nfreq)))
freq1=np.linspace(fmin,fmax,nfreq1)
## dt s.t. 10% of the points of a regular grid
rand_ind=np.sort(random.sample(range(nt),nt-nt_del))
t2=t[rand_ind]
freq2_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t2,freq1,D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
x2=A*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(nt-nt_del)
sigma_sq2=fact_var2*np.var(x2)
myperiodogram2=np.zeros(freq2_new.size)
signif_level=np.zeros(freq2_new.size)
for k in trange(freq2_new.size):
	M2=LS_WOSA(t2,freq2_new[k],0,zero_array2,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
	myperiodogram2[k]=la.norm(np.dot(x2,M2))**2/float(nsmooth_vec[0])
	mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
	moment1=np.trace(mymat)
	moment2=la.norm(mymat,ord='fro')**2
	M=moment1**2/moment2
	g=moment2/moment1
	signif_level[k]=g*chi2distr.ppf(proba,M)
myperiodogram2_bis=copy.copy(myperiodogram2)
myperiodogram2/=(sigma_sq2*signif_level)   # => threshold = 1.
maxper=np.amax(myperiodogram2)
minper=np.amin(myperiodogram2)
indexes2=peakutils.indexes(myperiodogram2,thres=(1.-minper)/(maxper-minper),min_dist=1)
plt.plot(freq2_new,myperiodogram2_bis,'k',label="WOSA tapered periodogram")
plt.plot(freq2_new,sigma_sq2*signif_level,'b',label="90% conf. level")
plt.plot(freq2_new[indexes2],myperiodogram2_bis[indexes2],'ro')
plt.xlim(freq2_new[0],freq2_new[-1])
plt.legend(fancybox=True,fontsize='small')
plt.xlabel("Frequency",fontsize=12)
plt.ylabel("Power",fontsize=12)
plt.tick_params(labelsize=12)
plt.savefig("final_npeaks_1sample_WOSA.pdf")
plt.show()

#### Tapered LS periodogram
signif_level=chi2distr.ppf(proba,2)
## dt s.t. 10% of the points of a regular grid
myperiodogram2=np.zeros(freq2_new.size)
for k in trange(freq2_new.size):
	myperiodogram2[k]=tapered_periodogram(t2,x2,freq2_new[k],taper_number)
maxper=np.amax(myperiodogram2)
minper=np.amin(myperiodogram2)
indexes2=peakutils.indexes(myperiodogram2,thres=(sigma_sq2*signif_level-minper)/(maxper-minper),min_dist=1)
plt.plot(freq2_new,myperiodogram2,'k',label="Tapered periodogram")
plt.plot(freq2_new,sigma_sq2*signif_level*np.ones(freq2_new.size),'b',label="90% conf. level")
plt.plot(freq2_new[indexes2],myperiodogram2[indexes2],'ro')
plt.xlim(freq2_new[0],freq2_new[-1])
plt.legend(fancybox=True,fontsize='small')
plt.xlabel("Frequency",fontsize=12)
plt.ylabel("Power",fontsize=12)
plt.tick_params(labelsize=12)
plt.savefig("final_npeaks_1sample_tapered.pdf")
plt.show()

#### Tapered LS periodogram (other taper)
signif_level=chi2distr.ppf(proba,2)
## dt s.t. 10% of the points of a regular grid
myperiodogram2=np.zeros(freq2_new.size)
for k in trange(freq2_new.size):
	myperiodogram2[k]=tapered_periodogram(t2,x2,freq2_new[k],taper_number2)
maxper=np.amax(myperiodogram2)
minper=np.amin(myperiodogram2)
indexes2=peakutils.indexes(myperiodogram2,thres=(sigma_sq2*signif_level-minper)/(maxper-minper),min_dist=1)
plt.plot(freq2_new,myperiodogram2,'k',label="Tapered periodogram")
plt.plot(freq2_new,sigma_sq2*signif_level*np.ones(freq2_new.size),'b',label="90% conf. level")
plt.plot(freq2_new[indexes2],myperiodogram2[indexes2],'ro')
plt.xlim(freq2_new[0],freq2_new[-1])
plt.legend(fancybox=True,fontsize='small')
plt.xlabel("Frequency",fontsize=12)
plt.ylabel("Power",fontsize=12)
plt.tick_params(labelsize=12)
plt.savefig("final_npeaks_1sample_tapered2.pdf")
plt.show()

### LS periodogram
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
## dt s.t. 10% of the points of a regular grid
myperiodogram2=np.zeros(freq2_new.size)
for k in trange(freq2_new.size):
	myperiodogram2[k]=periodogram(t2,x2,freq2_new[k])
maxper=np.amax(myperiodogram2)
minper=np.amin(myperiodogram2)
indexes2=peakutils.indexes(myperiodogram2,thres=(sigma_sq2*signif_level-minper)/(maxper-minper),min_dist=1)
plt.plot(freq2_new,myperiodogram2,'k',label="LS periodogram")
plt.plot(freq2_new,sigma_sq2*signif_level*np.ones(freq2_new.size),'b',label="90% conf. level")
plt.plot(freq2_new[indexes2],myperiodogram2[indexes2],'ro')
plt.xlim(freq2_new[0],freq2_new[-1])
plt.legend(fancybox=True,fontsize='small')
plt.xlabel("Frequency",fontsize=12)
plt.ylabel("Power",fontsize=12)
plt.tick_params(labelsize=12)
plt.savefig("final_npeaks_1sample_classical.pdf")
plt.show()
