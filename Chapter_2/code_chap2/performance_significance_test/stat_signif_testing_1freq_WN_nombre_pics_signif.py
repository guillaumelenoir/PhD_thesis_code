
# signal = white noise + 1 period

# computes the average number of false "peak" detection PER FREQUENCY_UNIT. It thus relies on a peak detection algorithm.

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
n_samples=1000				#1000
proba=0.9
A=np.linspace(0.,2.,num=21)  #21
t=np.arange(1000.)
nt_del=900
nfreq=2000	# Nombre de points de frequence sur l'intervalle [0.,1./2./dt]
beta=0.75	# for the WOSA
D_WOSA=200.		# for the WOSA
taper_number=8	# for the tapered per. and the WOSA
taper_number2=10
coverage=90.	# for the WOSA

### LS periodogram
dt=t[1]-t[0]
nA=A.size
nt=t.size
signif_level=chi2distr.ppf(proba,2)
npeaks11=np.zeros(A.size)
npeaks12=np.zeros(A.size)
npeaks13=np.zeros(A.size)
fact_var=float(nt)/float(nt+1)
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
fmin=1./(t[-1]-t[0])
fmax=1./2./dt
nfreq1=int(np.floor((fmax-fmin)/(1./2./dt)*float(nfreq)))
freq1=np.linspace(fmin,fmax,nfreq1)
freq1_matrix0=np.zeros((nfreq1,t.size))
freq1_matrix1=np.zeros((nfreq1,t.size))
for k in trange(nfreq1):
	mymat=periodogram_grid_fixed(t,freq1[k])
	freq1_matrix0[k,:]=mymat[:,0]
	freq1_matrix1[k,:]=mymat[:,1]
for j in trange(nA):
	count=0
	count2=0.
	count3=0.
	for k in trange(n_samples):
		## dt=cste
		freq_signal=np.random.uniform(low=fmin,high=fmax)
		x=A[j]*np.sin(2.*np.pi*freq_signal*t)+np.random.standard_normal(nt)
		sigma_sq=fact_var*np.var(x)
		myperiodogram=np.tensordot(freq1_matrix0,x,axes=([1],[0]))**2+np.tensordot(freq1_matrix1,x,axes=([1],[0]))**2 # use of tensordot is faster than computing at each frequency independently
		maxper=np.amax(myperiodogram)
		minper=np.amin(myperiodogram)
		indexes=peakutils.indexes(myperiodogram,thres=(sigma_sq*signif_level-minper)/(maxper-minper),min_dist=1)
		count+=indexes.size
		## dt s.t. 10% of the points of a regular grid
		rand_ind=np.sort(random.sample(range(nt),nt-nt_del))
		t2=t[rand_ind]
		D=t2[-1]-t2[0]
		fmin2=1./D
		fmax2=1./2./max(dt_av_central(t2,D,1),dt_av(t2,D,1))
		freq_signal=np.random.uniform(low=fmin2,high=fmax2)
		x2=A[j]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(nt-nt_del)
		nfreq2=int(np.floor((fmax2-fmin2)/(1./2./dt)*float(nfreq)))
		freq2=np.linspace(fmin2,fmax2,nfreq2)
		sigma_sq2=fact_var2*np.var(x2)
		myperiodogram2=np.zeros(nfreq2)
		for k in range(nfreq2):
			myperiodogram2[k]=periodogram(t2,x2,freq2[k])
		maxper=np.amax(myperiodogram2)
		minper=np.amin(myperiodogram2)
		indexes2=peakutils.indexes(myperiodogram2,thres=(sigma_sq2*signif_level-minper)/(maxper-minper),min_dist=1)
		count2+=float(indexes2.size)/(fmax2-fmin2)
#		## linear interpolation
#		nt_interp=int((t2[-1]-t2[0])/dt)+1
#		t_interp=np.linspace(t2[0],t2[-1],nt_interp)
#		f=interp1d(t2,x2)
#		x3=f(t_interp)
#		fmin3=1./D
#		nfreq3=int(np.floor((fmax-fmin3)/(1./2./dt)*float(nfreq)))
#		freq3=np.linspace(fmin3,fmax,nfreq3)
#		fact_var3=float(nt_interp)/float(nt_interp+1)
#		sigma_sq3=fact_var3*np.var(x3)
#		myperiodogram3=np.zeros(nfreq3)
#		for k in range(nfreq3):
#			myperiodogram3[k]=periodogram(t_interp,x3,freq3[k])
#		maxper=np.amax(myperiodogram3)
#		minper=np.amin(myperiodogram3)
#		indexes3=peakutils.indexes(myperiodogram3,thres=(sigma_sq3*signif_level-minper)/(maxper-minper),min_dist=1)
#		count3+=float(indexes3.size)/(fmax-fmin3)

#		plt.plot(freq1,myperiodogram,'b')
#		plt.plot(freq1,sigma_sq*signif_level*np.ones(freq1.size),'b')
#		plt.plot(freq1[indexes],myperiodogram[indexes],'kx')
#		plt.plot(freq2,myperiodogram2,'g')
#		plt.plot(freq2,sigma_sq2*signif_level*np.ones(freq2.size),'g')
#		plt.plot(freq2[indexes2],myperiodogram2[indexes2],'kx')
#		plt.show()
#		time.sleep(5)

	npeaks11[j]=float(count)/float(n_samples)/(fmax-fmin)
	npeaks12[j]=count2/float(n_samples)
#	npeaks13[j]=count3/float(n_samples)
npeaks11/=1./2./dt     # Nombre de pics sur l'intervalle [0,1./2./dt]
npeaks12/=1./2./dt
#plt.plot(A**2/2.,npeaks11,label="Reg. sampled")
#plt.plot(A**2/2.,npeaks12,label="Irreg. sampled")     # It is large because the time series carries much less data points
#plt.plot(A**2/2.,npeaks13,label="Interp. linear")
#plt.legend()
#plt.show()

#### Tapered LS periodogram
dt=t[1]-t[0]
nA=A.size
nt=t.size
signif_level=chi2distr.ppf(proba,2)
npeaks21=np.zeros(A.size)
npeaks22=np.zeros(A.size)
#npeaks23=np.zeros(A.size)
npeaks24=np.zeros(A.size)
npeaks25=np.zeros(A.size)
fact_var=float(nt)/float(nt+1)
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
fmin=1./(t[-1]-t[0])
fmax=1./2./dt
nfreq1=int(np.floor((fmax-fmin)/(1./2./dt)*float(nfreq)))
freq1=np.linspace(fmin,fmax,nfreq1)
freq1_matrix0=np.zeros((nfreq1,t.size))
freq1_matrix1=np.zeros((nfreq1,t.size))
for k in trange(nfreq1):
	mymat=tapered_periodogram_grid_fixed(t,freq1[k],taper_number)
	freq1_matrix0[k,:]=mymat[:,0]
	freq1_matrix1[k,:]=mymat[:,1]
freq1_matrix0_bis=np.zeros((nfreq1,t.size))
freq1_matrix1_bis=np.zeros((nfreq1,t.size))
for k in trange(nfreq1):
	mymat=tapered_periodogram_grid_fixed(t,freq1[k],taper_number2)
	freq1_matrix0_bis[k,:]=mymat[:,0]
	freq1_matrix1_bis[k,:]=mymat[:,1]
for j in trange(nA):
	count=0
	count2=0.
	#count3=0.
	count4=0
	count5=0.
	for k in trange(n_samples):
		## dt=cste
		freq_signal=np.random.uniform(low=fmin,high=fmax)
		x=A[j]*np.sin(2.*np.pi*freq_signal*t)+np.random.standard_normal(nt)
		sigma_sq=fact_var*np.var(x)
		myperiodogram=np.tensordot(freq1_matrix0,x,axes=([1],[0]))**2+np.tensordot(freq1_matrix1,x,axes=([1],[0]))**2 # use of tensordot is faster than computing at each frequency independently
		maxper=np.amax(myperiodogram)
		minper=np.amin(myperiodogram)
		indexes=peakutils.indexes(myperiodogram,thres=(sigma_sq*signif_level-minper)/(maxper-minper),min_dist=1)
		count+=indexes.size
		#
		myperiodogram4=np.tensordot(freq1_matrix0_bis,x,axes=([1],[0]))**2+np.tensordot(freq1_matrix1_bis,x,axes=([1],[0]))**2
		maxper=np.amax(myperiodogram4)
		minper=np.amin(myperiodogram4)
		indexes=peakutils.indexes(myperiodogram4,thres=(sigma_sq*signif_level-minper)/(maxper-minper),min_dist=1)
		count4+=indexes.size
		## dt s.t. 10% of the points of a regular grid
		rand_ind=np.sort(random.sample(range(nt),nt-nt_del))
		t2=t[rand_ind]
		D=t2[-1]-t2[0]
		fmin2=1./D
		fmax2=1./2./max(dt_av_central(t2,D,taper_number),dt_av(t2,D,taper_number))
		freq_signal=np.random.uniform(low=fmin2,high=fmax2)
		x2=A[j]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(nt-nt_del)
		nfreq2=int(np.floor((fmax2-fmin2)/(1./2./dt)*float(nfreq)))
		freq2=np.linspace(fmin2,fmax2,nfreq2)
		sigma_sq2=fact_var2*np.var(x2)
		myperiodogram2=np.zeros(nfreq2)
		for k in range(nfreq2):
			myperiodogram2[k]=tapered_periodogram(t2,x2,freq2[k],taper_number)
		maxper=np.amax(myperiodogram2)
		minper=np.amin(myperiodogram2)
		indexes2=peakutils.indexes(myperiodogram2,thres=(sigma_sq2*signif_level-minper)/(maxper-minper),min_dist=1)
		count2+=float(indexes2.size)/(fmax2-fmin2)
		#
		fmin5=1./D
		fmax5=1./2./max(dt_av_central(t2,D,taper_number2),dt_av(t2,D,taper_number2))
		freq_signal=np.random.uniform(low=fmin5,high=fmax5)
		x5=A[j]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(nt-nt_del)
		nfreq5=int(np.floor((fmax5-fmin5)/(1./2./dt)*float(nfreq)))
		freq5=np.linspace(fmin5,fmax5,nfreq5)
		sigma_sq5=fact_var2*np.var(x5)
		myperiodogram5=np.zeros(nfreq5)
		for k in range(nfreq5):
			myperiodogram5[k]=tapered_periodogram(t2,x5,freq5[k],taper_number2)
		maxper=np.amax(myperiodogram5)
		minper=np.amin(myperiodogram5)
		indexes5=peakutils.indexes(myperiodogram5,thres=(sigma_sq5*signif_level-minper)/(maxper-minper),min_dist=1)
		count5+=float(indexes5.size)/(fmax5-fmin5)
#		## linear interpolation
#		nt_interp=int((t2[-1]-t2[0])/dt)+1
#		t_interp=np.linspace(t2[0],t2[-1],nt_interp)
#		f=interp1d(t2,x2)
#		x3=f(t_interp)
#		fmin3=1./D
#		nfreq3=int(np.floor((fmax-fmin3)/(1./2./dt)*float(nfreq)))
#		freq3=np.linspace(fmin3,fmax,nfreq3)
#		fact_var3=float(nt_interp)/float(nt_interp+1)
#		sigma_sq3=fact_var3*np.var(x3)
#		myperiodogram3=np.zeros(nfreq3)
#		for k in range(nfreq3):
#			myperiodogram3[k]=tapered_periodogram(t_interp,x3,freq3[k],taper_number)
#		maxper=np.amax(myperiodogram3)
#		minper=np.amin(myperiodogram3)
#		indexes3=peakutils.indexes(myperiodogram3,thres=(sigma_sq3*signif_level-minper)/(maxper-minper),min_dist=1)
#		count3+=float(indexes3.size)/(fmax-fmin3)
	npeaks21[j]=float(count)/float(n_samples)/(fmax-fmin)
	npeaks24[j]=float(count4)/float(n_samples)/(fmax-fmin)
	npeaks22[j]=count2/float(n_samples)
	npeaks25[j]=count5/float(n_samples)
#	npeaks23[j]=count3/float(n_samples)
npeaks21/=1./2./dt     # Nombre de pics sur l'intervalle [0,1./2./dt]
npeaks24/=1./2./dt
npeaks22/=1./2./dt
npeaks25/=1./2./dt
#plt.plot(A**2/2.,npeaks11,label="Reg. sampled")
#plt.plot(A**2/2.,npeaks12,label="Irreg. sampled")     # It is large because the time series carries much less data points
#plt.plot(A**2/2.,npeaks13,label="Interp. linear")
#plt.legend()
#plt.show()
plt.plot(A**2/2.,npeaks11,'k',label="RS")
plt.plot(A**2/2.,npeaks21,'k-.',label="RS - Tapered (1)")
plt.plot(A**2/2.,npeaks24,'k--',label="RS - Tapered (2)")
plt.plot(A**2/2.,npeaks12,'b',label="IS")
plt.plot(A**2/2.,npeaks22,'b-.',label="IS - Tapered (1)")
plt.plot(A**2/2.,npeaks25,'b--',label="IS - Tapered (2)")
plt.xlim(0.,2.)
plt.legend(fancybox=True,fontsize='small')
plt.xlabel("SNR",fontsize=12)
plt.ylabel("# detected peaks",fontsize=14)
plt.title("False peak detection",fontsize=14)
plt.tick_params(labelsize=12)
plt.savefig("final_npeaks_classic_vs_tapered.pdf")
plt.close()

#### WOSA tapered LS periodogram - 2-moment approx for the signif. levels
dt=t[1]-t[0]
nA=A.size
nt=t.size
signif_level=chi2distr.ppf(proba,2)
npeaks31=np.zeros(A.size)
npeaks32=np.zeros(A.size)
npeaks33=np.zeros(A.size)
fact_var=float(nt)/float(nt+1)
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
zero_array=np.zeros((nt,2))
zero_array2=np.zeros((nt-nt_del,2))
fmin=1./D_WOSA
fmax=1./2./dt
nfreq1=int(np.floor((fmax-fmin)/(1./2./dt)*float(nfreq)))
freq1=np.linspace(fmin,fmax,nfreq1)    # This will be used as a prior frequency grid in all the following cases
for j in trange(nA):
	count=0.
	count2=0.
	count3=0.
	for k in trange(n_samples):
		## dt=cste
		freq1_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t,freq1,D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
		freq_signal=np.random.uniform(low=freq1_new[0],high=freq1_new[-1])
		x=A[j]*np.sin(2.*np.pi*freq_signal*t)+np.random.standard_normal(nt)
		sigma_sq=fact_var*np.var(x)
		myperiodogram=np.zeros(freq1_new.size)
		signif_level=np.zeros(freq1_new.size)
		for k in range(freq1_new.size):
			M2=LS_WOSA(t,freq1_new[k],0,zero_array,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
			myperiodogram[k]=la.norm(np.dot(x,M2))**2/float(nsmooth_vec[0])
			mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
			moment1=np.trace(mymat)
			moment2=la.norm(mymat,ord='fro')**2
			M=moment1**2/moment2
			g=moment2/moment1
			signif_level[k]=g*chi2distr.ppf(proba,M)
		
#		myperiodogram_bis=copy.copy(myperiodogram)
		
		myperiodogram/=(sigma_sq*signif_level)   # => threshold = 1.
		maxper=np.amax(myperiodogram)
		minper=np.amin(myperiodogram)
		indexes=peakutils.indexes(myperiodogram,thres=(1.-minper)/(maxper-minper),min_dist=1)
		count+=float(indexes.size)/(freq1_new[-1]-freq1_new[0])
		
#		plt.plot(freq1_new,myperiodogram_bis,'b')
#		plt.plot(freq1_new,sigma_sq*signif_level,'b')
#		plt.plot(freq1_new[indexes],myperiodogram_bis[indexes],'kx')
#		plt.show()
#		time.sleep(5)
		
		## dt s.t. 10% of the points of a regular grid
		rand_ind=np.sort(random.sample(range(nt),nt-nt_del))
		t2=t[rand_ind]
		freq2_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t2,freq1,D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
		freq_signal=np.random.uniform(low=freq2_new[0],high=freq2_new[-1])
		x2=A[j]*np.sin(2.*np.pi*freq_signal*t2)+np.random.standard_normal(nt-nt_del)
		sigma_sq2=fact_var2*np.var(x2)
		myperiodogram2=np.zeros(freq2_new.size)
		signif_level=np.zeros(freq2_new.size)
		for k in range(freq2_new.size):
			M2=LS_WOSA(t2,freq2_new[k],0,zero_array2,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
			myperiodogram2[k]=la.norm(np.dot(x2,M2))**2/float(nsmooth_vec[0])
			mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
			moment1=np.trace(mymat)
			moment2=la.norm(mymat,ord='fro')**2
			M=moment1**2/moment2
			g=moment2/moment1
			signif_level[k]=g*chi2distr.ppf(proba,M)
		
#		myperiodogram2_bis=copy.copy(myperiodogram2)
		
		myperiodogram2/=(sigma_sq2*signif_level)   # => threshold = 1.
		maxper=np.amax(myperiodogram2)
		minper=np.amin(myperiodogram2)
		indexes2=peakutils.indexes(myperiodogram2,thres=(1.-minper)/(maxper-minper),min_dist=1)
		count2+=float(indexes2.size)/(freq2_new[-1]-freq2_new[0])

#		plt.plot(freq2_new,myperiodogram2_bis,'b')
#		plt.plot(freq2_new,sigma_sq2*signif_level,'b')
#		plt.plot(freq2_new[indexes2],myperiodogram2_bis[indexes2],'kx')
#		plt.show()
#		time.sleep(5)

#		## linear interpolation
#		nt_interp=int((t2[-1]-t2[0])/dt)+1
#		t_interp=np.linspace(t2[0],t2[-1],nt_interp)
#		f=interp1d(t2,x2)
#		x3=f(t_interp)
#		fact_var3=float(nt_interp)/float(nt_interp+1)
#		sigma_sq3=fact_var3*np.var(x3)
#		freq3_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t_interp,freq1,D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
#		freq_signal=np.random.uniform(low=freq3_new[0],high=freq3_new[-1])
#		myperiodogram3=np.zeros(freq3_new.size)
#		signif_level=np.zeros(freq3_new.size)
#		zero_array3=np.zeros((nt_interp,2))
#		for k in range(freq3_new.size):
#			M2=LS_WOSA(t_interp,freq3_new[k],0,zero_array3,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#			myperiodogram3[k]=la.norm(np.dot(x3,M2))**2/float(nsmooth_vec[0])
#			mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#			moment1=np.trace(mymat)
#			moment2=la.norm(mymat,ord='fro')**2
#			M=moment1**2/moment2
#			g=moment2/moment1
#			signif_level[k]=g*chi2distr.ppf(proba,M)
#		myperiodogram3/=(sigma_sq3*signif_level)   # => threshold = 1.
#		maxper=np.amax(myperiodogram3)
#		minper=np.amin(myperiodogram3)
#		indexes3=peakutils.indexes(myperiodogram3,thres=(1.-minper)/(maxper-minper),min_dist=1)
#		count3+=float(indexes3.size)/(freq3_new[-1]-freq3_new[0])
	npeaks31[j]=float(count)/float(n_samples)
	npeaks32[j]=count2/float(n_samples)
#	npeaks33[j]=count3/float(n_samples)
npeaks31/=1./2./dt     # Nombre de pics sur l'intervalle [0,1./2./dt]
npeaks32/=1./2./dt
#plt.plot(A**2/2.,npeaks31,label="Reg. sampled")
#plt.plot(A**2/2.,npeaks32,label="Irreg. sampled")     # It is large because the time series carries much less data points
#plt.plot(A**2/2.,npeaks33,label="Interp. linear")
#plt.legend()
#plt.show()
plt.plot(A**2/2.,npeaks11,'k',label="RS")
plt.plot(A**2/2.,npeaks21,'k-.',label="RS - Tapered")
plt.plot(A**2/2.,npeaks31,'k--',label="RS - Tapered+WOSA")
plt.plot(A**2/2.,npeaks12,'b',label="IS")
plt.plot(A**2/2.,npeaks22,'b-.',label="IS - Tapered")
plt.plot(A**2/2.,npeaks32,'b--',label="IS - Tapered+WOSA")
plt.xlim(0.,2.)
plt.legend(fancybox=True,fontsize='small')
plt.xlabel("SNR",fontsize=12)
plt.ylabel("# detected peaks",fontsize=14)
plt.title("False peak detection",fontsize=14)
plt.tick_params(labelsize=12)
plt.savefig("final_npeaks_classic_vs_tapered_vs_WOSA.pdf")
plt.close()
