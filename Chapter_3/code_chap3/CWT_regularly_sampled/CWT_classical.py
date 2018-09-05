import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2distr

def CWT_Morlet_regularly_sampled(t,y,deltaj=0.05,w0=5.5,Tmin=None,Tmax=None,percentile_cwt=None,alpha_rn=None,sigma_wn=None,t_axis_label="",mydata_axis_label="",t_units=None,mydata_units=None,time_string=None,period_string=None,power_string=None,dashed_periods=None,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,color_cl=None,linewidth_cl=2,cmap="jet",nlevels=50,reverse_xaxis=False,reverse_yaxis=False,decimals=3,plot_coi="fill",linewidth_coi=1.0,plot_perlim2="fill",linewidth_perlim2=1.0):
	
	""" Scalogram of the CWT-Morlet with NO SMOOTHING, in the case of a regularly sampled time series.
		Significance testing against a red noise is possible, if its coefficients are provided (alpha_rn and sigma_wn)
		Normalization c(s)=1/s, such that the |CWT|**2 represents the local squared amplitude.
		The support of exp(-x**2/2) is approximated to be of length 6, and centered around 0.
		
		Required inputs:
		- t [1-dim numpy array of floats]: the times of the time series, distinct and in ascending order.
		- mydata [1-dim numpy array of floats - size=time.size]: the data at the times given by 't'.
		Optional Inputs:
		- deltaj=0.05: parameter controlling the density of periods/scales. Scale vector is defined as scale_min*2.0**(float(j)*deltaj) for j=0, 1, 2, ... A smaller deltaj implies denser scales and a better precision. We have the following approximate relation between periods and scales: period=2*np.pi*scale.
		- w0=5.5: the usual parameter for the Morlet wavelet controlling the time-frequency resolution. Minimal allowed value is w0=5.5. Increasing w0 offers a better scale/period resolution but a worse time resolution. There is always a trade-off due to the so-called 'Heisenberg uncertainty principle'.
		- Tmin=None: minimal period. Default value is permin=2.*dt, where dt is the time step.
		- Tmax=None: maximal period. Default value is permax=np.pi*(t[-1]-t[0])/w0/3.
		- percentile_cwt=None: The x^th percentiles for the confidence levels. Must be a 1-dim numpy array. Default is the 95^th percentile (i.e. the 95% confidence level):percentile_cwt=np.zeros(1); percentile_cwt[0]=95.0
		- alpha_rn=None and sigma_wn=None: AR-1 coefficients (sigma_wn is the white noise st. dev.). Must be provided in the case we want to perform significance testing against a red noise. The easiest way to get it is to go through wavepal first until "carma_params" method. Example with trendless data:
			import wavepal as wv
			x=wv.Wavepal(t,y)
			x.check_data()
			x.choose_trend_degree(-1)
			x.trend_vectors()
			x.carma_params(signif_level_type="a")
			sig_rn=x.sigwn_unique
			alpha_rn=x.alpha_unique
		- t_axis_label="": label for the time axis in figures
		- mydata_axis_label="": label for the data axis in figures
		- t_units=None: units of 't' (string type)
		- mydata_units=None: units of 'mydata' (string type).
		- time_string=None: list of floats containing the location of the ticks for the time axis.
		- period_string=None: list of floats containing the location of the ticks for the period axis.
		- dashed_periods=None: list of floats containing the periods for which a dashed line is drawn.
		- fontsize_title=14: fontsize for the figure title.
		- fontsize_axes=12: fontsize for the figure axes.
		- fontsize_ticks=12: fontsize for the figure ticks.
		- color_cl=None: colors for the confidence levels/contours. Ex.: ['m','g','r'] will draw in magenta, green and red the 3 levels specified in 'percentile_cwt'. Default is a magenta contour for all the confidence levels.
		- linewidth_cl=2: linewidth of the confidence contours in the scalogram
		- cmap="jet": colormap for the scalogram. Other choices on http://matplotlib.org/users/colormaps.html
		- nlevels=50: number of automatically-chosen color levels.
		- reverse_xaxis=False: Reverse the horizontal axis if True
		- reverse_yaxis=False: Reverse the vertical axis if True
		- decimals=3: Numbers of decimals for the colorbar scale ticks.
		- plot_coi="fill": plotting-type for the cone of influence: "fill" or "line".
		- linewidth_coi=1.0: linewidth for the coi (if plot_coi="line").
		- plot_perlim2="fill": plotting-type for the refinement of the SNEZ: "fill" or "line".
		- linewidth_perlim2=1.0: linewidth for perlim2 (if plot_perlim2="line").
		-----------------------------
		(C) 2017 G. Lenoir"""
	
	# Time step
	dt=t[1]-t[0]

	# Evaluate the some input arguments
	if Tmin is None:
		Tmin=2.*dt
	if Tmax is None:
		Tmax=np.pi*(t[-1]-t[0])/3./w0
	if percentile_cwt is None:
		percentile_cwt=np.zeros(1)
		percentile_cwt[0]=95.0
	if t_units is None:
		t_label=""
	else:
		t_label=" ("+t_units+")"
	if mydata_units is None:
		power_label=""
	else:
		power_label=" ("+mydata_units+"${}^2$)"
	if color_cl is None:
		color_cl=['m']*percentile_cwt.size

	# CWT_Morlet
	s0=Tmin*w0/(2.*np.pi)
	smax=Tmax*w0/(2.*np.pi)
	n=y.size
	ffty=np.fft.fft(y)
	J=int(np.floor(np.log2(smax/s0)/deltaj))+1
	s=np.zeros(J)
	for j in range(J):
		s[j]=s0*2.0**(float(j)*deltaj)
	wavelet=np.zeros(n)
	CWT=np.zeros((n,J),dtype='complex128')
	for j in range(J):
		for k in range(n):
			wavelet[k]=np.exp((-((2.*np.pi*float(k)/(float(n)*dt))*s[j]-w0)**2)/2.)
		conv=ffty*wavelet
		CWT[:,j]=np.fft.ifft(conv)
	CWT*=2.
	coi1=3.*s+t[0]
	coi2=t[-1]-3.*s
	periodcwt=2.*np.pi*s/w0
	WPS=np.abs(CWT)**2

	# Analytical confidence levels red noise
	if ((alpha_rn is not None) and (sigma_wn is not None)):
		WPS_cl_anal=np.zeros((n,J,percentile_cwt.size))
		spectrum_rn=(sigma_wn**2)/alpha_rn/s*dt/np.sqrt(np.pi)*(1.0-np.exp(-alpha_rn*dt)**2)/(1.0+np.exp(-alpha_rn*dt)**2-2.0*np.exp(-alpha_rn*dt)*np.cos(2.*np.pi/periodcwt*dt))
		for l in range(percentile_cwt.size):
			for k in range(n):
				WPS_cl_anal[k,:,l]=spectrum_rn/2.0*chi2distr.ppf(percentile_cwt[l]/100.,2)

	# Refining the SNEZ
	perlim2=2.*dt*(1.++3./np.sqrt(2.)/w0)*np.ones(n)

	# Max and min scalogram - does not take into account the coi and the ref. of the SNEZ
	j_perlim2=max(int(np.ceil(np.log2(w0*perlim2[0]/2./np.pi/s0)/deltaj)),0)
	ind=int(np.ceil(3.*s[j_perlim2]/dt))
	maxscal=np.max(WPS[ind:(n-ind),j_perlim2])
	minscal=np.min(WPS[ind:(n-ind),j_perlim2])
	#WPS[:ind,0]=0.
	#WPS[(n-ind):,0]=0.
	for j in range(j_perlim2+1,J):
		ind=int(np.ceil(3.*s[j]/dt))
		maxscal=max(maxscal,np.max(WPS[ind:(n-ind),j]))
		minscal=min(minscal,np.min(WPS[ind:(n-ind),j]))
		#WPS[:ind,k]=0.
		#WPS[(n-ind):,k]=0.

	# Figure scalogram
	if period_string is None:
		myperiod=np.ceil(np.log2(periodcwt[0]))
		myperiod=2.**myperiod
		period_string=[myperiod]
		while True:
			myperiod*=2.
			if myperiod>periodcwt[-1]:
				break
			period_string.append(myperiod)
	mycontourf=plt.contourf(t,np.log2(periodcwt),np.transpose(WPS),nlevels,cmap=cmap,vmin=minscal,vmax=maxscal)
	if (alpha_rn is not None) and (sigma_wn is not None):
		for k in range(percentile_cwt.size):
			cv_anal=WPS/WPS_cl_anal[:,:,k]
			plt.contour(t,np.log2(periodcwt),np.transpose(cv_anal),levels=[1.],colors=color_cl[k],linewidths=linewidth_cl)
	if plot_coi=="fill":
		plt.fill_betweenx(np.log2(periodcwt),t[0],coi1,edgecolors=None,facecolor='black',alpha=0.5)
		plt.fill_betweenx(np.log2(periodcwt),coi2,t[-1],edgecolors=None,facecolor='black',alpha=0.5)
	elif plot_coi=="line":
		plt.plot(coi1,np.log2(periodcwt),'k',linewidth=linewidth_coi)
		plt.plot(coi2,np.log2(periodcwt),'k')
	if plot_perlim2=="fill":
		plt.fill_between(t,np.log2(periodcwt[0])*np.ones(n),np.log2(perlim2),edgecolors=None,facecolor='black',alpha=0.5)
	elif plot_perlim2=="line":
			plt.plot(t,np.log2(perlim2),'k',linewidth=linewidth_perlim2)
	ax = plt.gca()
	ax.tick_params(length=5, width=1, color='w')
	if dashed_periods is not None:
		for k in range(len(dashed_periods)):
			plt.plot(t,np.log2(dashed_periods[k])*np.ones(n),'w--')
	plt.xlabel(t_axis_label+t_label,fontsize=fontsize_axes)
	plt.ylabel("Period"+t_label,fontsize=fontsize_axes)
	mytitle="Wavelet Scalogram"
	if percentile_cwt.size>0:
		mytitle+=" & "
		for k in range(percentile_cwt.size-2):
			mytitle+=str(percentile_cwt[k])
			mytitle+="%, "
		if percentile_cwt.size>=2:
			mytitle+=str(percentile_cwt[-2])
			mytitle+="% and "
			mytitle+=str(percentile_cwt[-1])
			mytitle+="% "
		else:
			mytitle+=str(percentile_cwt[-1])
			mytitle+="% "
		mytitle+="Confidence levels"
	if time_string is not None:
		plt.xticks(time_string, time_string)
	plt.xticks(fontsize=fontsize_ticks)
	plt.yticks(np.log2(period_string), period_string, fontsize=fontsize_ticks)
	plt.xlim([t[0],t[-1]])
	plt.ylim([np.log2(periodcwt[0]), np.log2(periodcwt[-1])])
	plt.suptitle(mytitle, fontsize=fontsize_title)
	if reverse_xaxis is True:
		plt.gca().invert_xaxis()
	if reverse_yaxis is True:
		plt.gca().invert_yaxis()
	# Colorbar and its rescaling: min and max of the levels of color are defined over the non-shaded regions
	cbar=plt.colorbar(mycontourf)
	cbar.ax.set_ylabel("Power"+power_label,fontsize=fontsize_axes)
	cbar.ax.tick_params(labelsize=fontsize_ticks)
	minminscal=np.amin(np.amin(WPS))
	maxmaxscal=np.amax(np.amax(WPS))
	cbar.set_clim(minminscal,maxmaxscal)
	my_color_ticks=np.linspace(minminscal,maxmaxscal,10)
	cbar.set_ticks(my_color_ticks)
	my_color_ticklabels=np.linspace(minscal,maxscal,10)
	my_power_color_ticklabels=float(10**(int(np.floor(np.log10(my_color_ticklabels[-1])))))
	my_color_ticklabels=[el/my_power_color_ticklabels for el in my_color_ticklabels]
	my_color_ticklabels=np.around(my_color_ticklabels,decimals=decimals)
	mystring=' '.join(str(e) for e in my_color_ticklabels).split(' ')
	if ('1' in mystring[0]) or ('2' in mystring[0]) or ('3' in mystring[0]) or ('4' in mystring[0]) or ('5' in mystring[0]) or ('6' in mystring[0]) or ('7' in mystring[0]) or ('8' in mystring[0]) or ('9' in mystring[0]):
		mystring[0]="<="+mystring[0]
	mystring[-1]=">="+mystring[-1]
	cbar.set_ticklabels(mystring)
	if int(np.log10(my_power_color_ticklabels))!=0:
		cbar.ax.set_title('1e'+str(int(np.log10(my_power_color_ticklabels))),fontsize=fontsize_ticks)
	return plt
