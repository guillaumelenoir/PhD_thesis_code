
# signal = white noise + 1 period

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

#params
n_samples=1000
proba=0.9
A=np.linspace(0.,np.sqrt(2),num=21)
t=np.arange(1000.)
nt_del=900
beta=0.75	# for the WOSA
D_WOSA=200.		# for the WOSA
taper_number=8	# for the tapered per. and the WOSA
coverage=90.	# for the WOSA

### LS periodogram
dt=t[1]-t[0]
nA=A.size
nt=t.size
signif_level=chi2distr.ppf(proba,2)
mypercentage11=np.zeros(A.size)
mypercentage12=np.zeros(A.size)
mypercentage13=np.zeros(A.size)
#mypercentage11_out=np.zeros(A.size)
#mypercentage12_out=np.zeros(A.size)
#mypercentage13_out=np.zeros(A.size)
fact_var=float(nt)/float(nt+1)
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
fmin=1./(t[-1]-t[0])
fmax=1./2./dt
for j in trange(nA):
	count=0
	count2=0
	count3=0
#	count_out=0
#	count2_out=0
#	count3_out=0
	for k in range(n_samples):
		## dt=cste
		freq=np.random.uniform(low=fmin,high=fmax)
		x=A[j]*np.sin(2.*np.pi*freq*t)+np.random.standard_normal(nt)
		sigma_sq=fact_var*np.var(x)
		myperiodogram=periodogram(t,x,freq)
		if myperiodogram>sigma_sq*signif_level:
			count+=1
		#
#		freq_out=np.random.uniform(low=fmin,high=fmax)
#		myperiodogram_out=periodogram(t,x,freq_out)
#		if myperiodogram_out>sigma_sq*signif_level:
#			count_out+=1
		## dt s.t. 10% of the points of a regular grid
		rand_ind=np.sort(random.sample(range(nt),nt-nt_del))
		t2=t[rand_ind]
		D=t2[-1]-t2[0]
		freq2=np.random.uniform(low=1./D,high=1./2./max(dt_av_central(t2,D,1),dt_av(t2,D,1)))
		x2=A[j]*np.sin(2.*np.pi*freq2*t2)+np.random.standard_normal(nt-nt_del)
		sigma_sq2=fact_var2*np.var(x2)
		myperiodogram2=periodogram(t2,x2,freq2)
		if myperiodogram2>sigma_sq2*signif_level:
			count2+=1
		#
#		freq2_out=np.random.uniform(low=1./D,high=1./2./max(dt_av_central(t2,D,1),dt_av(t2,D,1)))
#		myperiodogram2_out=periodogram(t2,x2,freq2_out)
#		if myperiodogram2_out>sigma_sq2*signif_level:
#			count2_out+=1
		## linear interpolation
		nt_interp=int((t2[-1]-t2[0])/dt)+1
		t_interp=np.linspace(t2[0],t2[-1],nt_interp)
		f=interp1d(t2,x2)
		x3=f(t_interp)
		fact_var3=float(nt_interp)/float(nt_interp+1)
		sigma_sq3=fact_var3*np.var(x3)
		myperiodogram3=periodogram(t_interp,x3,freq2)
		if myperiodogram3>sigma_sq3*signif_level:
			count3+=1
		#
#		freq3_out=np.random.uniform(low=1./float(nt_interp)/dt,high=1./2./dt)
#		myperiodogram3_out=periodogram(t_interp,x3,freq3_out)
#		if myperiodogram3_out>sigma_sq3*signif_level:
#			count3_out+=1
		# spline (cubic) interpolation       -      NUMERICALLY UNSTABLE (may diverge)
		#time1=time.time()
		#f=interp1d(t2,x2,kind='cubic')
		#x4=f(t_interp)
		#sigma_sq4=fact_var3*np.var(x4)
		#myperiodogram4=periodogram(t_interp,x4,freq2)
		#if myperiodogram4>sigma_sq4*signif_level:
		#	count4+=1
	mypercentage11[j]=float(count)/float(n_samples)*100.
	mypercentage12[j]=float(count2)/float(n_samples)*100.
	mypercentage13[j]=float(count3)/float(n_samples)*100.
#	mypercentage11_out[j]=float(count_out)/float(n_samples)*100.
#	mypercentage12_out[j]=float(count2_out)/float(n_samples)*100.
#	mypercentage13_out[j]=float(count3_out)/float(n_samples)*100.
#plt.plot(A**2/2.,mypercentage11,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage12,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage13,label="Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
#plt.plot(A**2/2.,mypercentage11_out,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage12_out,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage13_out,label="Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.legend()
#plt.show()

### Tapered LS periodogram
dt=t[1]-t[0]
nA=A.size
nt=t.size
signif_level=chi2distr.ppf(proba,2)
mypercentage21=np.zeros(A.size)
mypercentage22=np.zeros(A.size)
mypercentage23=np.zeros(A.size)
#mypercentage21_out=np.zeros(A.size)
#mypercentage22_out=np.zeros(A.size)
#mypercentage23_out=np.zeros(A.size)
fact_var=float(nt)/float(nt+1)
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
fmin=1./(t[-1]-t[0])
fmax=1./2./dt
for j in trange(nA):
	count=0
	count2=0
	count3=0
#	count_out=0
#	count2_out=0
#	count3_out=0
	for k in range(n_samples):
		## dt=cste
		freq=np.random.uniform(low=fmin,high=fmax)
		x=A[j]*np.sin(2.*np.pi*freq*t)+np.random.standard_normal(nt)
		sigma_sq=fact_var*np.var(x)
		myperiodogram=tapered_periodogram(t,x,freq,taper_number)
		if myperiodogram>sigma_sq*signif_level:
			count+=1
		#
#		freq_out=np.random.uniform(low=fmin,high=fmax)
#		myperiodogram_out=tapered_periodogram(t,x,freq_out,taper_number)
#		if myperiodogram_out>sigma_sq*signif_level:
#			count_out+=1
		## dt s.t. 10% of the points of a regular grid
		rand_ind=np.sort(random.sample(range(nt),nt-nt_del))
		t2=t[rand_ind]
		D=t2[-1]-t2[0]
		freq2=np.random.uniform(low=1./D,high=1./2./max(dt_av_central(t2,D,taper_number),dt_av(t2,D,taper_number)))
		x2=A[j]*np.sin(2.*np.pi*freq2*t2)+np.random.standard_normal(nt-nt_del)
		sigma_sq2=fact_var2*np.var(x2)
		myperiodogram2=tapered_periodogram(t2,x2,freq2,taper_number)
		if myperiodogram2>sigma_sq2*signif_level:
			count2+=1
		#
#		freq2_out=np.random.uniform(low=1./D,high=1./2./max(dt_av_central(t2,D,taper_number),dt_av(t2,D,taper_number)))
#		myperiodogram2_out=tapered_periodogram(t2,x2,freq2_out,taper_number)
#		if myperiodogram2_out>sigma_sq2*signif_level:
#			count2_out+=1
		## linear interpolation
		nt_interp=int((t2[-1]-t2[0])/dt)+1
		t_interp=np.linspace(t2[0],t2[-1],nt_interp)
		f=interp1d(t2,x2)
		x3=f(t_interp)
		fact_var3=float(nt_interp)/float(nt_interp+1)
		sigma_sq3=fact_var3*np.var(x3)
		myperiodogram3=tapered_periodogram(t_interp,x3,freq2,taper_number)
		if myperiodogram3>sigma_sq3*signif_level:
			count3+=1
		#
#		freq3_out=np.random.uniform(low=1./float(nt_interp)/dt,high=fmax)
#		myperiodogram3_out=tapered_periodogram(t_interp,x3,freq3_out,taper_number)
#		if myperiodogram3_out>sigma_sq3*signif_level:
#			count3_out+=1
	mypercentage21[j]=float(count)/float(n_samples)*100.
	mypercentage22[j]=float(count2)/float(n_samples)*100.
	mypercentage23[j]=float(count3)/float(n_samples)*100.
#	mypercentage21_out[j]=float(count_out)/float(n_samples)*100.
#	mypercentage22_out[j]=float(count2_out)/float(n_samples)*100.
#	mypercentage23_out[j]=float(count3_out)/float(n_samples)*100.
#plt.plot(A**2/2.,mypercentage21,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage22,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage23,label="Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
#plt.plot(A**2/2.,mypercentage21_out,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage22_out,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage23_out,label="Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.legend()
#plt.show()

### WOSA tapered LS periodogram - 2-moment approx for the signif. levels
dt=t[1]-t[0]
nA=A.size
nt=t.size
mypercentage31=np.zeros(A.size)
mypercentage32=np.zeros(A.size)
mypercentage33=np.zeros(A.size)
#mypercentage31_out=np.zeros(A.size)
#mypercentage32_out=np.zeros(A.size)
#mypercentage33_out=np.zeros(A.size)
#mypercentage31_ampl=np.zeros(A.size)
#mypercentage32_ampl=np.zeros(A.size)
#mypercentage33_ampl=np.zeros(A.size)
#mypercentage31_ampl_out=np.zeros(A.size)
#mypercentage32_ampl_out=np.zeros(A.size)
#mypercentage33_ampl_out=np.zeros(A.size)
fact_var=float(nt)/float(nt+1)
fact_var2=float(nt-nt_del)/float(nt-nt_del+1)
zero_array=np.zeros((nt,2))
zero_array2=np.zeros((nt-nt_del,2))
fmin=1./D_WOSA
fmax=1./2./dt
for j in trange(nA):
	count=0
	count2=0
	count3=0
#	count_out=0
#	count2_out=0
#	count3_out=0
#	count_ampl=0
#	count2_ampl=0
#	count3_ampl=0
#	count_ampl_out=0
#	count2_ampl_out=0
#	count3_ampl_out=0
	for k in trange(n_samples):
		## dt=cste
		freq_new=np.zeros(0)
		while freq_new.size==0:
			freq=np.random.uniform(low=fmin,high=fmax)
			freq_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t,freq*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
		M2=LS_WOSA(t,freq,0,zero_array,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
		x=A[j]*np.sin(2.*np.pi*freq*t)+np.random.standard_normal(nt)
		sigma_sq=fact_var*np.var(x)
		myperiodogram=la.norm(np.dot(x,M2))**2/float(nsmooth_vec[0])
		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
		moment1=np.trace(mymat)
		moment2=la.norm(mymat,ord='fro')**2
		M=moment1**2/moment2
		g=moment2/moment1
		signif_level=g*chi2distr.ppf(proba,M)
		if myperiodogram>sigma_sq*signif_level:
			count+=1
		#
#		freq_new=np.zeros(0)
#		while freq_new.size==0:
#			freq_out=np.random.uniform(low=fmin,high=fmax)
#			freq_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t,freq_out*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
#		M2=LS_WOSA(t,freq_out,0,zero_array,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram_out=la.norm(np.dot(x,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram_out>sigma_sq*signif_level:
#			count_out+=1
#		#
#		_,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t,freq*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,True)
#		M2=LS_WOSA(t,freq,0,zero_array,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram_ampl=la.norm(np.dot(x,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram_ampl>sigma_sq*signif_level:
#			count_ampl+=1
#		#
#		_,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t,freq_out*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,True)
#		M2=LS_WOSA(t,freq_out,0,zero_array,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram_ampl_out=la.norm(np.dot(x,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram_ampl_out>sigma_sq*signif_level:
#			count_ampl_out+=1
		## dt s.t. 10% of the points of a regular grid
		rand_ind=np.sort(random.sample(range(nt),nt-nt_del))
		t2=t[rand_ind]
		freq_new=np.zeros(0)
		while freq_new.size==0:
			freq2=np.random.uniform(low=fmin,high=fmax)
			freq_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t2,freq2*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
		M2=LS_WOSA(t2,freq2,0,zero_array2,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
		x2=A[j]*np.sin(2.*np.pi*freq2*t2)+np.random.standard_normal(nt-nt_del)
		sigma_sq2=fact_var2*np.var(x2)
		myperiodogram2=la.norm(np.dot(x2,M2))**2/float(nsmooth_vec[0])
		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
		moment1=np.trace(mymat)
		moment2=la.norm(mymat,ord='fro')**2
		M=moment1**2/moment2
		g=moment2/moment1
		signif_level=g*chi2distr.ppf(proba,M)
		if myperiodogram2>sigma_sq2*signif_level:
			count2+=1
		#
#		freq_new=np.zeros(0)
#		while freq_new.size==0:
#			freq2_out=np.random.uniform(low=fmin,high=fmax)
#			freq_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t2,freq2_out*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
#		M2=LS_WOSA(t2,freq2_out,0,zero_array2,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram2_out=la.norm(np.dot(x2,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram2_out>sigma_sq2*signif_level:
#			count2_out+=1
#		#
#		_,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t2,freq2*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,True)
#		M2=LS_WOSA(t2,freq2,0,zero_array2,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram2_ampl=la.norm(np.dot(x2,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram2_ampl>sigma_sq2*signif_level:
#			count2_ampl+=1
#		#
#		_,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t2,freq2_out*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,True)
#		M2=LS_WOSA(t2,freq2_out,0,zero_array2,D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram2_ampl_out=la.norm(np.dot(x2,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram2_ampl_out>sigma_sq2*signif_level:
#			count2_ampl_out+=1
		## linear interpolation
		nt_interp=int((t2[-1]-t2[0])/dt)+1
		t_interp=np.linspace(t2[0],t2[-1],nt_interp)
		f=interp1d(t2,x2)
		x3=f(t_interp)
		fact_var3=float(nt_interp)/float(nt_interp+1)
		sigma_sq3=fact_var3*np.var(x3)
		_,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t_interp,freq2*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
		M2=LS_WOSA(t_interp,freq2,0,np.zeros((nt_interp,2)),D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
		myperiodogram3=la.norm(np.dot(x3,M2))**2/float(nsmooth_vec[0])
		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
		moment1=np.trace(mymat)
		moment2=la.norm(mymat,ord='fro')**2
		M=moment1**2/moment2
		g=moment2/moment1
		signif_level=g*chi2distr.ppf(proba,M)
		if myperiodogram3>sigma_sq3*signif_level:
			count3+=1
		#
#		freq_new=np.zeros(0)
#		while freq_new.size==0:
#			freq3_out=np.random.uniform(low=fmin,high=fmax)
#			freq_new,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t_interp,freq3_out*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,False)
#		M2=LS_WOSA(t_interp,freq3_out,0,np.zeros((nt_interp,2)),D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram3_out=la.norm(np.dot(x3,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram3_out>sigma_sq3*signif_level:
#			count3_out+=1
#		#
#		_,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t_interp,freq2*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,True)
#		M2=LS_WOSA(t_interp,freq2,0,np.zeros((nt_interp,2)),D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram3_ampl=la.norm(np.dot(x3,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram3_ampl>sigma_sq3*signif_level:
#			count3_ampl+=1
#		#
#		_,tau,myind_time,myind_freq,myind_Q,D_new,nsmooth_vec,nsmooth,weight_WOSA=freq_analysis_prelims(t_interp,freq3_out*np.ones(1),D_WOSA,beta,taper_number,coverage,True,True,-1,None,True)
#		M2=LS_WOSA(t_interp,freq3_out,0,np.zeros((nt_interp,2)),D_new,tau,nsmooth,nsmooth_vec[0],myind_time,myind_freq,taper_number,-1,weight_WOSA)
#		myperiodogram3_ampl_out=la.norm(np.dot(x3,M2))**2/float(nsmooth_vec[0])
#		mymat=np.dot(np.transpose(M2),M2)/float(nsmooth_vec[0])
#		moment1=np.trace(mymat)
#		moment2=la.norm(mymat,ord='fro')**2
#		M=moment1**2/moment2
#		g=moment2/moment1
#		signif_level=g*chi2distr.ppf(proba,M)
#		if myperiodogram3_ampl_out>sigma_sq3*signif_level:
#			count3_ampl_out+=1
	mypercentage31[j]=float(count)/float(n_samples)*100.
	mypercentage32[j]=float(count2)/float(n_samples)*100.
	mypercentage33[j]=float(count3)/float(n_samples)*100.
#	mypercentage31_out[j]=float(count_out)/float(n_samples)*100.
#	mypercentage32_out[j]=float(count2_out)/float(n_samples)*100.
#	mypercentage33_out[j]=float(count3_out)/float(n_samples)*100.
#	mypercentage31_ampl[j]=float(count_ampl)/float(n_samples)*100.
#	mypercentage32_ampl[j]=float(count2_ampl)/float(n_samples)*100.
#	mypercentage33_ampl[j]=float(count3_ampl)/float(n_samples)*100.
#	mypercentage31_ampl_out[j]=float(count_ampl_out)/float(n_samples)*100.
#	mypercentage32_ampl_out[j]=float(count2_ampl_out)/float(n_samples)*100.
#	mypercentage33_ampl_out[j]=float(count3_ampl_out)/float(n_samples)*100.
#plt.plot(A**2/2.,mypercentage31,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage32,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33,label="Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
#plt.plot(A**2/2.,mypercentage31_out,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage32_out,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33_out,label="Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
#plt.plot(A**2/2.,mypercentage11,label="Classic. Reg. sampled")
#plt.plot(A**2/2.,mypercentage12,label="Classic. Irreg. sampled")
#plt.plot(A**2/2.,mypercentage13,label="Classic. Interp. linear")
#plt.plot(A**2/2.,mypercentage21,label="Tapered Reg. sampled")
#plt.plot(A**2/2.,mypercentage22,label="Tapered Irreg. sampled")
#plt.plot(A**2/2.,mypercentage23,label="Tapered Interp. linear")
#plt.plot(A**2/2.,mypercentage31,label="WOSA Reg. sampled")
#plt.plot(A**2/2.,mypercentage32,label="WOSA Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33,label="WOSA Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
#plt.plot(A**2/2.,mypercentage11_out,label="Classic. Reg. sampled")
#plt.plot(A**2/2.,mypercentage12_out,label="Classic. Irreg. sampled")
#plt.plot(A**2/2.,mypercentage13_out,label="Classic. Interp. linear")
#plt.plot(A**2/2.,mypercentage21_out,label="Tapered Reg. sampled")
#plt.plot(A**2/2.,mypercentage22_out,label="Tapered Irreg. sampled")
#plt.plot(A**2/2.,mypercentage23_out,label="Tapered Interp. linear")
#plt.plot(A**2/2.,mypercentage31_out,label="WOSA Reg. sampled")
#plt.plot(A**2/2.,mypercentage32_out,label="WOSA Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33_out,label="WOSA Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
#plt.plot(A**2/2.,mypercentage31,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage32,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33,label="Interp. linear")
#plt.plot(A**2/2.,mypercentage31_ampl,label="Weighted Reg. sampled")
#plt.plot(A**2/2.,mypercentage32_ampl,label="Weighted Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33_ampl,label="Weighted Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
#plt.plot(A**2/2.,mypercentage31_out,label="Reg. sampled")
#plt.plot(A**2/2.,mypercentage32_out,label="Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33_out,label="Interp. linear")
#plt.plot(A**2/2.,mypercentage31_ampl_out,label="Weighted Reg. sampled")
#plt.plot(A**2/2.,mypercentage32_ampl_out,label="Weighted Irreg. sampled")
#plt.plot(A**2/2.,mypercentage33_ampl_out,label="Weighted Interp. linear")
#plt.plot(A**2/2.,10.*np.ones(A.size),label="10% detection line")
#plt.ylim(0.,100.)
#plt.legend()
#plt.show()
plt.plot(A**2/2.,mypercentage11,'k',label="RS")
plt.plot(A**2/2.,mypercentage21,'k-.',label="RS - Tapered")
plt.plot(A**2/2.,mypercentage31,'k--',label="RS - Tapered+WOSA")
plt.plot(A**2/2.,mypercentage12,'b',label="IS")
plt.plot(A**2/2.,mypercentage22,'b-.',label="IS - Tapered")
plt.plot(A**2/2.,mypercentage32,'b--',label="IS - Tapered+WOSA")
plt.plot(A**2/2.,mypercentage13,'r',label="IN")
plt.plot(A**2/2.,mypercentage23,'r-.',label="IN - Tapered")
plt.plot(A**2/2.,mypercentage33,'r--',label="IN - Tapered+WOSA")
plt.plot(A**2/2.,10.*np.ones(A.size),'m',label="10% rejection")
plt.xlim(0.,1.)
plt.ylim(0.,100.)
plt.legend(fancybox=True,fontsize='small',loc='lower right')
plt.xlabel("SNR",fontsize=12)
plt.ylabel("% rejection of H_0",fontsize=14)
plt.title("Power of the test",fontsize=14)
plt.tick_params(labelsize=12)
plt.savefig("final_power_classic_vs_tapered_vs_WOSA.pdf")


