import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from tapering_window_wv33 import tapering_window

mypath="../figures/"

data=np.genfromtxt("data/ODP1148-BF-18O.txt")
t=data[:,0]
t=t-t[0]
D=t[-1]  # NB: t[0]=0. (easier for using tapering_window)
freqstep=0.001
freq=np.arange(1./t[-1],0.5-1./t[-1],freqstep)  # dt_GCD=1 for that data set
# Tapers
g1=tapering_window(t,D,1)    # square
g2=tapering_window(t,D,4)	 # sin^2
g3=tapering_window(t,D,10)	 # Gaussian

nf=freq.size
app_g_cos_sin=np.zeros((nf,3))
app_g_cos2=np.zeros((nf,3))
app_gsq_cos2=np.zeros((nf,3))
for k in trange(nf):
	mysin2=np.sin(4.*np.pi*freq[k]*t)
	mycos2=np.cos(4.*np.pi*freq[k]*t)
	# computes beta
	beta1=np.arctan(np.sum(g1**2*mysin2)/np.sum(g1**2*mycos2))/2.0
	beta2=np.arctan(np.sum(g2**2*mysin2)/np.sum(g2**2*mycos2))/2.0
	beta3=np.arctan(np.sum(g3**2*mysin2)/np.sum(g3**2*mycos2))/2.0
	# for the 3 panels
	# compares the numerator with the denominator: to show how the numerator is small compared to the denominator, in %
	app_g_cos_sin[k,0]=np.sum(g1*np.cos(2.*np.pi*freq[k]*t-beta1)*np.sin(2.*np.pi*freq[k]*t-beta1))/np.sum(g1)*2.*100.
	app_g_cos_sin[k,1]=np.sum(g2*np.cos(2.*np.pi*freq[k]*t-beta2)*np.sin(2.*np.pi*freq[k]*t-beta2))/np.sum(g2)*2.*100.
	app_g_cos_sin[k,2]=np.sum(g3*np.cos(2.*np.pi*freq[k]*t-beta3)*np.sin(2.*np.pi*freq[k]*t-beta3))/np.sum(g3)*2.*100.
	app_g_cos2[k,0]=np.sum(g1*np.cos(4.*np.pi*freq[k]*t-2.*beta1))/np.sum(g1)*100.
	app_g_cos2[k,1]=np.sum(g2*np.cos(4.*np.pi*freq[k]*t-2.*beta2))/np.sum(g2)*100.
	app_g_cos2[k,2]=np.sum(g3*np.cos(4.*np.pi*freq[k]*t-2.*beta3))/np.sum(g3)*100.
	app_gsq_cos2[k,0]=np.sum(g1**2*np.cos(4.*np.pi*freq[k]*t-2.*beta1))/np.sum(g1**2)*100.
	app_gsq_cos2[k,1]=np.sum(g2**2*np.cos(4.*np.pi*freq[k]*t-2.*beta2))/np.sum(g2**2)*100.
	app_gsq_cos2[k,2]=np.sum(g3**2*np.cos(4.*np.pi*freq[k]*t-2.*beta3))/np.sum(g3**2)*100.
# Figures
plt.subplot(311)
plt.plot(freq,np.absolute(app_g_cos2[:,0]))
plt.plot(freq,np.absolute(app_g_cos2[:,1]))
plt.plot(freq,np.absolute(app_g_cos2[:,2]))
plt.text(0.01,17.,'(a)')
plt.ylabel("Abs. error (%)")
plt.ylim([0.,20.])
plt.subplot(312)
plt.plot(freq,np.absolute(app_gsq_cos2[:,0]))
plt.plot(freq,np.absolute(app_gsq_cos2[:,1]))
plt.plot(freq,np.absolute(app_gsq_cos2[:,2]))
plt.text(0.01,17.,'(b)')
plt.ylabel("Abs. error (%)")
plt.ylim([0.,20.])
plt.subplot(313)
plt.plot(freq,np.absolute(app_g_cos_sin[:,0]))
plt.plot(freq,np.absolute(app_g_cos_sin[:,1]))
plt.plot(freq,np.absolute(app_g_cos_sin[:,2]))
plt.text(0.01,17.,'(c)')
plt.xlabel("freq. (kyr${}^{-1}$)")
plt.ylabel("Abs. error (%)")
plt.ylim([0.,20.])
plt.savefig(mypath+"approx_Lomb_generalized.pdf")
plt.close()

