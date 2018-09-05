import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from tqdm import trange
from scipy.stats import chi2 as chi2distr
from periodogram import periodogram
from WOSA_periodogram import WOSA_periodogram
from periodogram_grid_fixed import periodogram_grid_fixed
from WOSA_periodogram_grid_fixed import WOSA_periodogram_grid_fixed
from WOSA_matrix import WOSA_matrix
from percentile_n_moments import percentile_n_moments
import time
import peakutils
#import copy
#import wavepal as wv

#### CASE 1: Regularly sampled - signal = white noise
#n_samples=100000
#proba=0.95
## 1 sample of the time series
#t1=np.arange(100.)
#x1=np.random.standard_normal(t1.size)
#plt.plot(t1,x1)
#plt.show()
## Stats about the significance
#freq=1./10.
#signif_level=chi2distr.ppf(proba,2)
#count=0
#for k in trange(n_samples):
#	x1=np.random.standard_normal(t1.size)
#	sigma_sq=float(t1.size)/float(t1.size+1)*np.var(x1)
#	myperiodogram=periodogram(t1,x1,freq)
#	if myperiodogram>sigma_sq*signif_level:
#		count+=1
#print float(count)/float(n_samples)*100.

### CASE 2: Regularly sampled - signal = white noise + 1 period
n_samples=10000
proba=0.95
A=np.linspace(0.,3.,num=100)
# 1 sample of the time series
t2=np.linspace(0.,999.,num=1000)     # t2.size must be even (easier)
freq_signal=250./float(t2.size)
freq_signal_int=250
x2=A[11]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(t2.size)
plt.plot(t2,x2)
plt.show()
# Periodogram components
ref_factor=100
freq_fine_grid=np.arange(0.,float(t2.size*ref_factor/2))/float(ref_factor*t2.size)
n_freq_fine_grid=freq_fine_grid.size
freq_matrix=np.zeros((n_freq_fine_grid,t2.size,2))
for k in trange(n_freq_fine_grid):
	freq_matrix[k,:,:]=periodogram_grid_fixed(t2,freq_fine_grid[k])
freq_signal_int*=ref_factor
# Stats about the significance
#freq_probed1=freq_signal
signif_level=chi2distr.ppf(proba,2)
mypercentage_in=np.zeros(A.size)
mypercentage_out=np.zeros(A.size)
mypercentage_out_fourier=np.zeros(A.size)
for j in trange(A.size):
	count_in=0
	count_out=0
	count_out_fourier=0
	for k in range(n_samples):
		x2=A[j]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(t2.size)
		sigma_sq=float(t2.size)/float(t2.size+1)*np.var(x2)
		myperiodogram=np.dot(freq_matrix[freq_signal_int,:,0],x2)**2+np.dot(freq_matrix[freq_signal_int,:,1],x2)**2
		if myperiodogram>sigma_sq*signif_level:
			count_in+=1
		#freq_probed2=np.random.uniform(low=0.,high=0.5)
		#while (freq_probed2<freq_signal+1./float(t2.size)) and (freq_probed2>freq_signal-1./float(t2.size)):
		#	freq_probed2=np.random.uniform(low=0.,high=0.5)
		#myperiodogram=periodogram(t2,x2,freq_probed2)
		freq_probed2_int=np.random.randint(low=1,high=n_freq_fine_grid)
		freq_probed2=freq_fine_grid[freq_probed2_int]
		#while (freq_probed2<freq_signal+1./float(t2.size)) and (freq_probed2>freq_signal-1./float(t2.size)):
		#	freq_probed2_int=np.random.randint(low=1,high=n_freq_fine_grid)
		#	freq_probed2=freq_fine_grid[freq_probed2_int]
		myperiodogram=np.dot(freq_matrix[freq_probed2_int,:,0],x2)**2+np.dot(freq_matrix[freq_probed2_int,:,1],x2)**2
		if myperiodogram>sigma_sq*signif_level:
			count_out+=1
		#freq_probed3_int=np.random.randint(low=1,high=t2.size/2)
		#while freq_probed3_int==freq_signal_int:
		#	freq_probed3_int=np.random.randint(low=1,high=t2.size/2)
		#freq_probed3=float(freq_probed3_int)/float(t2.size)
		#myperiodogram=periodogram(t2,x2,freq_probed3)
		freq_probed3_int=np.random.randint(low=1,high=t2.size/2)*ref_factor
		while freq_probed3_int==freq_signal_int:
			freq_probed3_int=np.random.randint(low=1,high=t2.size/2)*ref_factor
		myperiodogram=np.dot(freq_matrix[freq_probed3_int,:,0],x2)**2+np.dot(freq_matrix[freq_probed3_int,:,1],x2)**2
		if myperiodogram>sigma_sq*signif_level:
			count_out_fourier+=1
		mypercentage_in[j]=float(count_in)/float(n_samples)*100.
		mypercentage_out[j]=float(count_out)/float(n_samples)*100.
		mypercentage_out_fourier[j]=float(count_out_fourier)/float(n_samples)*100.
plt.plot(A,mypercentage_in,label="in")
plt.plot(A,mypercentage_out,label="out")
plt.plot(A,mypercentage_out_fourier,label="out-Fourier")
plt.plot(A,5.*np.ones(A.size),label="5% detection line")
plt.legend()
plt.show()
mynumber_peaks=np.zeros(A.size)
n_samples_peaks=100
myperiodogram=np.zeros(n_freq_fine_grid)
y2=np.zeros((t2.size,2))
for j in trange(A.size):
	count_peaks=0
	for k in range(n_samples_peaks):
		#print k
		x2=A[j]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(t2.size)
		y2[:,0]=x2[:]
		y2[:,1]=x2[:]
		sigma_sq=float(t2.size)/float(t2.size+1)*np.var(x2)
		max_per=0
		#for l in range(n_freq_fine_grid):
			#time1=time.time()
		#	myperiodogram[l]=np.dot(freq_matrix[l,:,0],x2)**2+np.dot(freq_matrix[l,:,1],x2)**2
			#time1=time.time()-time1
			#print time1
			#time2=time.time()
		#	max_per=max(myperiodogram[l],max_per)
			#time2=time.time()-time2
			#print time2
			#time3=time.time()
		#	if (freq_fine_grid[l]<freq_signal+1./float(t2.size)) and (freq_fine_grid[l]>freq_signal-1./float(t2.size)):
		#		myperiodogram[l]=0.
			#time3=time.time()-time3
			#print time3
			#print ''
		myperiodogram=np.tensordot(freq_matrix,y2,axes=([1,2],[0,1]))**2
		max_per=np.amax(myperiodogram)
		indexes=peakutils.indexes(myperiodogram,thres=sigma_sq*signif_level/max_per,min_dist=1)
		count_peaks+=indexes.size
		#plt.plot(freq_fine_grid,myperiodogram)
		#plt.plot(freq_fine_grid,sigma_sq*signif_level*np.ones(n_freq_fine_grid))
		#plt.plot(freq_fine_grid[indexes],myperiodogram[indexes],'r.')
		#plt.show()
	mynumber_peaks[j]=float(count_peaks)/float(n_samples_peaks)
plt.plot(A,mynumber_peaks)
plt.legend()
plt.show()

#### CASE 3: Irregularly sampled - signal = white noise + 1 period
#n_samples=1000
#proba=0.95
#freq_signal=float(10)/100.           # A faire: freq_signal random, mais de Fourier, a chaque iteration
#A=np.linspace(0.,3.,num=50)
## 1 sample of the time series
#t2=np.arange(100.)
#for k in range(90):
#	t2=np.delete(t2,np.random.randint(low=1,high=100-k-1))
#print t2
#x2=np.random.standard_normal(t2.size)+A[11]*np.sin(2.*np.pi*freq_signal*t2)
#plt.plot(t2,x2)
#plt.show()
## Stats about the significance
#freq_probed1=freq_signal
#signif_level=chi2distr.ppf(proba,2)
#mypercentage_in=np.zeros(A.size)
#mypercentage_out=np.zeros(A.size)
#mypercentage_out_fourier=np.zeros(A.size)
#for j in trange(A.size):
#	count_in=0
#	count_out=0
#	count_out_fourier=0
#	for k in range(n_samples):
#		t2=np.arange(100.)
#		for k in range(90):
#			t2=np.delete(t2,np.random.randint(low=1,high=100-k-1))
#		x2=np.random.standard_normal(t2.size)+A[j]*np.sin(2.*np.pi*freq_signal*t2)
#		myperiodogram=periodogram(t2,x2,freq_probed1)
#		if myperiodogram>signif_level:
#			count_in+=1
#		freq_probed2=np.random.uniform(low=0.,high=0.5)
#		while (freq_probed2<freq_signal+1./float(t2.size)) and (freq_probed2>freq_signal-1./float(t2.size)):
#			freq_probed2=np.random.uniform(low=0.,high=0.5)
#		myperiodogram=periodogram(t2,x2,freq_probed2)
#		if myperiodogram>signif_level:
#			count_out+=1
#		freq_probed3=float(np.random.randint(low=1,high=100))/100.
#		while freq_probed3==freq_signal:
#			freq_probed3=float(np.random.randint(low=1,high=100))/100.
#		myperiodogram=periodogram(t2,x2,freq_probed3)
#		if myperiodogram>signif_level:
#			count_out_fourier+=1
#		mypercentage_in[j]=float(count_in)/float(n_samples)*100.
#		mypercentage_out[j]=float(count_out)/float(n_samples)*100.
#		mypercentage_out_fourier[j]=float(count_out_fourier)/float(n_samples)*100.
#plt.plot(A,mypercentage_in,label="in")
#plt.plot(A,mypercentage_out,label="out")
#plt.plot(A,mypercentage_out_fourier,label="out-Fourier")
#plt.legend()
#plt.show()

### CASE 4: Regularly sampled - signal = white noise + 1 period - WOSA periodogram
n_samples=10000
t2=np.linspace(0.,999.,num=1000)     # t2.size must be even (easier)
beta=0.75
D=100  # number of data points for the WOSA segments ; t2.size/D doit etre entier et D doit etre entier.
freq_signal=25./float(D)     # 0.25
freq_signal_int=25
proba=0.95
n_moments=10
A=np.linspace(0.,3.,num=100)
# Stats about the significance
M2=WOSA_matrix(t2,freq_signal,beta,D)
Q=M2.shape[1]/2
myeig=la.eigvalsh(np.dot(np.transpose(M2),M2))
mytrace=np.zeros((1,n_moments))
for o in range(n_moments):
	mytrace[0,o]=np.sum(myeig**(o+1))
signif_level,check_percentiles=percentile_n_moments(mytrace,proba*np.ones(1),0*np.ones(1,dtype='int'),'gamma-polynomial',100000)
plt.plot(range(2,n_moments+1),check_percentiles[0,:,0])
plt.show()
## check that the signif. level is idem at all the frequencies
#myfreq=np.linspace(0.,0.5,num=1001)
#mytrace=np.zeros((myfreq.size,n_moments))
#for k in trange(myfreq.size):
#	M2=WOSA_matrix(t2,myfreq[k],beta,D)
#	myeig=la.eigvalsh(np.dot(np.transpose(M2),M2))
#	for o in range(n_moments):
#		mytrace[k,o]=np.sum(myeig**(o+1))
#signif_level,_=percentile_n_moments(mytrace,proba*np.ones(1),0*np.ones(1,dtype='int'),'gamma-polynomial',100000)
#plt.plot(myfreq,signif_level)
#plt.show()
# Periodogram components
ref_factor=100
freq_fine_grid=np.arange(0.,float(D*ref_factor/2))/float(ref_factor*D)
n_freq_fine_grid=freq_fine_grid.size
freq_matrix=np.zeros((n_freq_fine_grid,t2.size,2*Q))
for k in trange(n_freq_fine_grid):
	freq_matrix[k,:,:]=WOSA_periodogram_grid_fixed(t2,freq_fine_grid[k],beta,D)
freq_signal_int*=ref_factor
mypercentage_in=np.zeros(A.size)
mypercentage_out=np.zeros(A.size)
mypercentage_out_fourier=np.zeros(A.size)
# 1 sample of the time series
x2=A[5]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(t2.size)
sigma_sq=float(t2.size)/float(t2.size+1)*np.var(x2)
plt.plot(t2,x2)
plt.show()
myperiodogram_no_smoothing=np.zeros(n_freq_fine_grid)
myperiodogram_with_smoothing=np.zeros(n_freq_fine_grid)
for k in range(n_freq_fine_grid):
	myperiodogram_no_smoothing[k]=periodogram(t2,x2,freq_fine_grid[k])
	for l in range(2*Q):
		myperiodogram_with_smoothing[k]+=np.dot(freq_matrix[k,:,l],x2)**2
	myperiodogram_with_smoothing[k]/=float(Q)
plt.plot(freq_fine_grid,myperiodogram_no_smoothing)
plt.plot(freq_fine_grid,sigma_sq*chi2distr.ppf(proba,2)*np.ones(n_freq_fine_grid))
plt.show()
plt.plot(freq_fine_grid,myperiodogram_with_smoothing)
plt.plot(freq_fine_grid,sigma_sq*signif_level[0]*np.ones(n_freq_fine_grid))
plt.show()
for j in trange(A.size):
	count_in=0
	count_out=0
	count_out_fourier=0
	for k in range(n_samples):
		x2=A[j]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(t2.size)
		sigma_sq=float(t2.size)/float(t2.size+1)*np.var(x2)
		myperiodogram=0.
		for l in range(2*Q):
			myperiodogram+=np.dot(freq_matrix[freq_signal_int,:,l],x2)**2
		myperiodogram/=float(Q)
		if myperiodogram>sigma_sq*signif_level:
			count_in+=1
		#freq_probed2=np.random.uniform(low=0.,high=0.5)
		#while (freq_probed2<freq_signal+1./float(D)) and (freq_probed2>freq_signal-1./float(D)):
		#	freq_probed2=np.random.uniform(low=0.,high=0.5)
		#myperiodogram=WOSA_periodogram(t2,x2,freq_probed2,beta,D)
		freq_probed2_int=np.random.randint(low=1,high=n_freq_fine_grid)
		freq_probed2=freq_fine_grid[freq_probed2_int]
		while (freq_probed2<freq_signal+1./float(D)) and (freq_probed2>freq_signal-1./float(D)):
			freq_probed2_int=np.random.randint(low=1,high=n_freq_fine_grid)
			freq_probed2=freq_fine_grid[freq_probed2_int]
		myperiodogram=0.
		for l in range(2*Q):
			myperiodogram+=np.dot(freq_matrix[freq_probed2_int,:,l],x2)**2
		myperiodogram/=float(Q)
		if myperiodogram>sigma_sq*signif_level:
			count_out+=1
		#freq_probed3_int=np.random.randint(low=1,high=D/2)
		#while freq_probed3_int==freq_signal_int:
		#	freq_probed3_int=np.random.randint(low=1,high=D/2)
		#freq_probed3=float(freq_probed3_int)/float(D)
		#myperiodogram=WOSA_periodogram(t2,x2,freq_probed3,beta,D)
		freq_probed3_int=np.random.randint(low=1,high=D/2)*ref_factor
		while freq_probed3_int==freq_signal_int:
			freq_probed3_int=np.random.randint(low=1,high=D/2)*ref_factor
		myperiodogram=0.
		for l in range(2*Q):
			myperiodogram+=np.dot(freq_matrix[freq_probed3_int,:,l],x2)**2
		myperiodogram/=float(Q)
		if myperiodogram>sigma_sq*signif_level:
			count_out_fourier+=1
		mypercentage_in[j]=float(count_in)/float(n_samples)*100.
		mypercentage_out[j]=float(count_out)/float(n_samples)*100.
		mypercentage_out_fourier[j]=float(count_out_fourier)/float(n_samples)*100.
plt.plot(A,mypercentage_in,label="in")
plt.plot(A,mypercentage_out,label="out")
plt.plot(A,mypercentage_out_fourier,label="out-Fourier")
plt.plot(A,5.*np.ones(A.size),label="5% detection line")
plt.legend()
plt.show()
