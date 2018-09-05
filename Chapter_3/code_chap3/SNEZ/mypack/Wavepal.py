import numpy as np

class Wavepal:



	def __init__(self,t,mydata,t_axis_label="",mydata_axis_label="",t_units=None,mydata_units=None):
		
		""" Constructor of Wavepal class. It Initializes all the variables of Wavepal class with an access outside Wavepal (the user can access them). The list of available variables is given in this function.
			Required Inputs:
			- t [1-dim numpy array of floats]: the times of the time series, distinct and in ascending order.
			- mydata [1-dim numpy array of floats - size=time.size]: the data at the times given by 't'.
			Optional Inputs:
			- t_axis_label="": label for the time axis in figures
			- mydata_axis_label="": label for the data axis in figures
			- t_units=None: units of 't' (string type)
			- mydata_units=None: units of 'mydata' (string type).
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		# check inputs
		try:
			assert np.issubsctype(t,float)
		except:
			print "Error at input 't': must be a numpy array of 'float' type"
			return
		try:
			assert np.issubsctype(mydata,float)
		except:
			print "Error at input 'mydata': must be a numpy array of 'float' type"
			return
		try:
			assert type(t_axis_label) is str
		except AssertionError:
			print "Error at input 't_axis_label': must be of 'str' type"
			return
		try:
			assert type(mydata_axis_label) is str
		except AssertionError:
			print "Error at input 'mydata_axis_label': must be of 'str' type"
			return
		try:
			assert (t_units is None) or (type(t_units) is str)
		except AssertionError:
			print "Error at input 't_units': must be None or of 'str' type"
			return
		try:
			assert (mydata_units is None) or (type(mydata_units) is str)
		except AssertionError:
			print "Error at input 'mydata_units': must be None or of 'str' type"
			return
		# Initializes all the variables the user may have access when using Wavepal
		self.t=t
		self.mydata=mydata
		self.t_axis_label=t_axis_label
		self.mydata_axis_label=mydata_axis_label
		self.t_units=t_units
		self.mydata_units=mydata_units
		if t_units is None:
			self.t_label=""
			self.freq_label=""
		else:
			self.t_label=" ("+t_units+")"
			self.freq_label=" ("+t_units+"${}^{-1}$)"
		if mydata_units is None:
			self.mydata_label=""
			self.power_label=""
			self.varpow_label=""
		else:
			self.mydata_label=" ("+mydata_units+")"
			self.power_label=" ("+mydata_units+"${}^2$)"
			self.varpow_label=" ("+mydata_units+"${}^4$)"
		# Defined in function 'check_data'
		self.nt=None
		self.run_check_data=False
		# Defined in function 'plot_timestep'
		self.dt=None
		# Defined in function 'choose_trend_degree'
		self.pol_degree=None
		self.trend=None
		self.run_choose_trend_degree=False
		# Defined in function 'trend_vectors'
		self.myprojvec=None
		self.Vmat=None
		self.run_trend_vectors=False
		# Defined in function 'timefreq_analysis'
		self.theta=None
		self.scale=None
		self.scalelim1_smooth=None
		self.shannonnyquistexclusionzone=None
		self.scalogram=None
		self.run_timefreq_analysis=False



	def check_data(self):
		
		""" check_data checks the data, and modifies them if needed. 
				-> check for NaN or Inf values. Returns an error if NaN or Inf found.
				-> Put the time axis in ascending order and check that all the times are distinct. Correction if needed.
				-> Check that there are at least 50 data points. Returns an error if this condition is not fulfilled.
			Inputs:
			/
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		from dt_min import dt_min
		from distinct_ages import distinct_ages

		# Check that the data are not Nan or Inf
		for k in range(self.t.size):
			if np.isnan(self.t[k]) or np.isnan(self.mydata[k]) or np.isinf(self.t[k]) or np.isinf(self.mydata[k]):
				print "Error: Nan or Inf value at time/age ", self.t[k]
				return
		# Put the time axis in ascending order and check that all the times are distinct - Correction if needed
		tind=np.argsort(self.t)
		self.t=self.t[tind]
		self.mydata=self.mydata[tind]
		try:
			assert dt_min(self.t)>0.0
		except AssertionError:
			print "WARNING: The times/ages of the time series must be distinct"
			print "The program will automatically select the first occurrence of the time/age (and skip the others) where they are the same (times/ages are previously set in ascending order)"
			self.t, self.mydata=distinct_ages(self.t,self.mydata)
		# Check that we have at least 50 data points
		self.nt=self.t.size
		try:
			assert self.nt>=50 # useful for the ACF of residual noise and its confidence levels (in CARMA pack)
		except AssertionError:
			print "Error: Not enough data points - please provide at least 50 data points"
			return
		self.run_check_data=True

	
	
	def plot_timestep(self,hist=True,nbins=10,log_yaxis=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,reverse_xaxis=False):
	
		""" plot_timestep computes and generates the figure of the time step in function of time.
			Optional Inputs:
			- hist=True: if True, draws the histogram of the distribution of the time steps.
			- nbins=10: number of bins for the histogram.
			- log_yaxis=False: If True, the vertical axis is in log scale.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- reverse_xaxis=False: Reverse the horizontal axis if True
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		import matplotlib.pyplot as plt
		
		# check inputs
		try:
			assert type(hist) is bool
		except AssertionError:
			print "Error at input 'hist': must be True or False"
			return
		try:
			assert (type(nbins) is int) and (nbins>0)
		except AssertionError:
			print "Error at input 'nbins': must be of type 'int' and >0"
			return
		try:
			assert type(log_yaxis) is bool
		except AssertionError:
			print "Error at input 'log_yaxis': must be True or False"
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print "Error at input 'reverse_xaxis': must be True or False"
			return
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print "Error: Must have run function 'check_data'"
			return
		self.dt=np.zeros(self.nt-1)
		for k in range(1,self.nt):
			self.dt[k-1]=self.t[k]-self.t[k-1]
		if log_yaxis is True:
			plt.semilogy(self.t[1:],self.dt,"k.",zorder=1)
		else:
			plt.plot(self.t[1:],self.dt,"k.",zorder=1)
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.xlim(self.t[1], self.t[-1])
		plt.ylabel(self.t_axis_label+" step"+self.t_label,fontsize=fontsize_axes)
		plt.title(self.t_axis_label+" step in function of "+self.t_axis_label,fontsize=fontsize_title)
		if hist is True:
			pdf, bin_edges = np.histogram(self.dt, bins=nbins)
			bin_edges = bin_edges[0:pdf.size]
			# Stretch the PDF so that it is readable on the residual plot when plotted horizontally
			pdf = pdf / float(pdf.max()) * 0.4 * (self.t[-1]-self.t[1])
			# Add the histogram to the plot
			plt.barh(bin_edges, pdf, height=bin_edges[1] - bin_edges[0],left=self.t[1],alpha=0.5,zorder=2)
		plt.tick_params(labelsize=fontsize_ticks)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		return plt
	
	

	def choose_trend_degree(self,pol_degree):

		""" choose_trend_degree records the user choice for the degree of the polynomial trend.
			Required Inputs:
			- pol_degree [int]: the degree of the polynomial trend. pol_degree=-1 means no trend (trend is imposed to zero). pol_degree=0 means trend=constant=data average. pol_degree=1 means linear detrending (trend=a*t+b). etc.
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		from detrending import detrending
		
		# check inputs
		try:
			assert (type(pol_degree) is int) and pol_degree>=-1
		except AssertionError:
			print "Error at input 'pol_degree': must be an integer >= -1"
			return
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print "Error: Must have run function 'check_data'"
			return
		self.pol_degree=pol_degree
		self.trend=detrending(self.t,self.mydata,self.pol_degree)
		self.run_choose_trend_degree=True



	def trend_vectors(self):
		
		""" trend_vectors computes some arrays concerning the trend of the time series.
			Inputs:
			/
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		import numpy.linalg as la
		import copy
		
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print "Error: Must have run function 'check_data'"
			return
		try:
			assert self.run_choose_trend_degree is True
		except AssertionError:
			print "Error: Must have run function 'choose_trend_degree'"
			return

		tbis=self.t/self.t[-1]    # for numerical stability with the powers of the time, for the polynomial trend
		myprojvec=np.zeros((self.nt,self.pol_degree+3))
		Vmat=np.zeros((self.nt,self.pol_degree+3))
		for o in range(self.pol_degree+1):
			myprojvec_o=tbis**o
			Vmat[:,o]=copy.copy(myprojvec_o)
			for k in range(o):	# Gram-Schmidt
				h=myprojvec[:,k]
				myprojvec_o-=np.dot(h,myprojvec_o)*h
			myprojvec[:,o]=myprojvec_o/la.norm(myprojvec_o)
		self.myprojvec=myprojvec
		self.Vmat=Vmat
		self.run_trend_vectors=True



	def timefreq_analysis(self,theta=None,permin=None,permax=None,deltaj=0.05,w0=5.5,gauss_spread=3.0,eps=0.000001,dt_GCD=None,shannonnyquistexclusionzone=True,weighted_CWT=True,smoothing_coeff=0.,smoothing_type="fixed",shannonnyquistexclusionzone_type="G"):
		
		import numpy.linalg as la
		from dt_min import dt_min
		from tqdm import trange
		from CWT import CWT
		from timefreq_analysis_prelims import timefreq_analysis_prelims
		
		# check inputs
		try:
			assert (theta is None) or np.issubsctype(theta,float)
		except:
			print "Error at input 'theta': must be None or a numpy array of 'float' type"
			return
		try:
			assert (permin is None) or ((type(permin) is float) and (permin>0.))
		except AssertionError:
			print "Error at input 'permin': must be of 'float' type and >0."
			return
		try:
			assert (permax is None) or ((type(permax) is float) and (permax>0.))
		except AssertionError:
			print "Error at input 'permax': must be of 'float' type and >0."
			return
		if (permin is not None) and (permax is not None):
			try:
				assert permax>=permin
			except AssertionError:
				print "Error at input 'permin' and 'permax': must have permax>=permin"
				return
		try:
			assert (type(deltaj) is float) and (deltaj>0.)
		except AssertionError:
			print "Error at input 'deltaj': must be of 'float' type and >0."
			return
		try:
			assert (type(w0) is float) and (w0>=5.5)
		except AssertionError:
			print "Error at input 'w0': must be of 'float' type and >=5.5"
			return
		try:
			assert (type(gauss_spread) is float) and (gauss_spread>0.)
		except AssertionError:
			print "Error at input 'gauss_spread': must be of 'float' type and >0."
			return
		try:
			assert (type(eps) is float) and (eps>=0.)
		except AssertionError:
			print "Error at input 'eps': must be of 'float' type and >=0."
			return
		try:
			assert (dt_GCD is None) or ((type(dt_GCD) is float) and (dt_GCD>0.))
		except AssertionError:
			print "Error at input 'dt_GCD': must be None or of 'float' type and >0."
			return
		try:
			assert type(shannonnyquistexclusionzone) is bool
		except AssertionError:
			print "Error at input 'shannonnyquistexclusionzone': Must be True or False"
			return
		try:
			assert type(weighted_CWT) is bool
		except AssertionError:
			print "Error at input 'weighted_CWT': Must be True or False"
			return
		try:
			assert (type(smoothing_coeff) is float) and (smoothing_coeff>=0.)
		except AssertionError:
			print "Error at input 'smoothing_coeff': must be of 'float' type and >=0."
			return
		try:
			assert (type(smoothing_type) is str) and ((smoothing_type.lower()=="fixed") or (smoothing_type.lower()=="variable"))
		except AssertionError:
			print "Error at input 'smoothing_type': Must be 'fixed' or 'variable'"
			return
		smoothing_type=smoothing_type.lower()
		try:
			assert (type(shannonnyquistexclusionzone_type) is str) and ((shannonnyquistexclusionzone_type.lower()=="g") or (shannonnyquistexclusionzone_type.lower()=="h") or (shannonnyquistexclusionzone_type.lower()=="gh"))
		except AssertionError:
			print "Error at input 'shannonnyquistexclusionzone_type': Must be 'g', 'h' or 'gh'"
			return
		shannonnyquistexclusionzone_type=shannonnyquistexclusionzone_type.lower()
		# Period -> Scale conversion for 'permin' and 'permax' if they are provided by the user
		if permin is None:
			scalemin=None
		else:
			if weighted_CWT is True:
				scalemin=permin/(2.*np.pi)
			elif weighted_CWT is False:
				scalemin=permin*(w0+np.sqrt(2+w0**2))/(4.*np.pi*w0)
		if permax is None:
			scalemax=None
		else:
			if weighted_CWT is True:
				scalemax=permax/(2.*np.pi)
			elif weighted_CWT is False:
				scalemax=permax*(w0+np.sqrt(2+w0**2))/(4.*np.pi*w0)
		# Set Default values for input arguments
		if dt_GCD is None:
			dt_GCD=dt_min(self.t)
		if scalemin is None:
			scalemin=dt_GCD/np.pi*(1.+eps)
		if scalemax is None:
			scalemax=(self.t[-1]-self.t[0])/2./w0/gauss_spread
		if theta is None:
			self.theta=np.linspace(self.t[0],self.t[-1],self.t.size)
		else:
			self.theta=theta
		# theta is put in ascending order and check whether theta range is included in self.t range
		theta_ind_sort=np.argsort(self.theta)
		self.theta=self.theta[theta_ind_sort]
		try:
			assert self.theta[0]>=self.t[0]
		except AssertionError:
			print "Error: theta[0] must be >= than the smallest time of the time series"
			return
		try:
			assert self.theta[-1]<=self.t[-1]
		except AssertionError:
			print "Error: theta[-1] must be <= than the biggest time of the time series"
			return
		# Adjust scalemin and scalemax and builds the scale vector
		scalemin=np.maximum(dt_GCD/np.pi*(1.+eps),scalemin)
		scalemax=np.minimum((self.t[-1]-self.t[0])/2./w0/gauss_spread,scalemax)
		J=int(np.floor(np.log2(scalemax/scalemin)/deltaj))
		scale=np.zeros(J+1)
		for k in range(J+1):
			scale[k]=scalemin*2.**(float(k)*deltaj)
		# Build the CWT components
		scale,_,_,_,_,_,weight_cwt,scalelim1_ind,scalelim1_smooth,scalelim1_ind_smooth,Qmax,n_outside_scalelim1=timefreq_analysis_prelims(self.t,self.theta,scale,w0,gauss_spread,eps,dt_GCD,shannonnyquistexclusionzone,weighted_CWT,smoothing_coeff,smoothing_type,shannonnyquistexclusionzone_type)
		J=scale.size
		Q=self.theta.size
		self.shannonnyquistexclusionzone=shannonnyquistexclusionzone
		self.scale=scale
		self.scalelim1_smooth=scalelim1_smooth
		# data scalogram - This is the main loop, over the scales
		self.scalogram=np.zeros((Q,J))
		mydata_transp=np.transpose(self.mydata)
		print "Main loop, over the time-frequency plane:"
		for l in trange(J):
			scale_l=scale[l]
			M2, mytheta, mytheta_ind=CWT(self.t,self.theta,l,scale_l,self.myprojvec,self.pol_degree,weight_cwt[:,l],scalelim1_ind,n_outside_scalelim1[l],w0)
			# Data scalogram - intermediate calculus
			scalogram_int=np.dot(mydata_transp,M2)
			for k in range(n_outside_scalelim1[l]):
				mytheta_k=mytheta[k]
				ind_left=np.argmin(np.absolute(mytheta-(mytheta_k-smoothing_coeff*w0*scale_l)))
				ind_right=np.argmin(np.absolute(mytheta-(mytheta_k+smoothing_coeff*w0*scale_l)))
				mylength=ind_right+1-ind_left
				# Data scalogram
				self.scalogram[mytheta_ind[k],l]=la.norm(scalogram_int[2*ind_left:2*(ind_right+1)])**2/float(mylength)
		# computes the min/max for the color scale
		minscal=float("inf")
		maxscal=-float("inf")
		for k in range(Q):
			minscal=min(minscal,np.amin(self.scalogram[k,scalelim1_ind_smooth[k]:]))
			maxscal=max(maxscal,np.amax(self.scalogram[k,scalelim1_ind_smooth[k]:]))
		self.minscal=minscal
		self.maxscal=maxscal

		self.run_timefreq_analysis=True



	def plot_scalogram(self,time_string=None,scale_string=None,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,cmap="jet",nlevels=50,reverse_xaxis=False,reverse_yaxis=False,alpha_SNEZ=0.5):
		
		""" plot_scalogram generates the figure of the scalogram and its confidence levels. It also generates the figure of the global scalogram and its confidence levels.
			Optional Inputs:
			- time_string=None: list of floats containing the location of the ticks for the time axis.
			- scale_string=None: list of floats containing the location of the ticks for the scale axis.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- cmap="jet": colormap for the scalogram. Other choices on http://matplotlib.org/users/colormaps.html
			- nlevels=50: number of automatically-chosen color levels.
			- reverse_xaxis=False: Reverse the horizontal axis if True
			- reverse_yaxis=False: Reverse the vertical axis if True
			- alpha_SNEZ=0.5: Transparency for the SNEZ. It must take a value between 0 (completely transparent) and 1 (completely opaque). Only used if shannonnyquistexclusionzone=False in the method 'timefreq_analysis'.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt

		# check inputs
		try:
			assert (time_string is None) or (type(time_string) is list)
		except AssertionError:
			print "Error at input 'time_string': must be None or of type 'list'"
			return
		if type(time_string) is list:
			for k in range(len(time_string)):
				try:
					assert type(time_string[k]) is float
				except AssertionError:
					print "Error at input 'time_string': must be a list containing floats"
					return
		try:
			assert (scale_string is None) or (type(scale_string) is list)
		except AssertionError:
			print "Error at input 'scale_string': must be None or of type 'list'"
			return
		if type(scale_string) is list:
			for k in range(len(scale_string)):
				try:
					assert type(scale_string[k]) is float
				except AssertionError:
					print "Error at input 'scale_string': must be a list containing floats"
					return
		try:
			assert type(cmap) is str
		except AssertionError:
			print "Error at input 'cmap': must be of type 'str'"
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print "Error at input 'reverse_xaxis': must be True or False"
			return
		try:
			assert type(reverse_yaxis) is bool
		except AssertionError:
			print "Error at input 'reverse_yaxis': must be True or False"
			return
		try:
			assert ((type(alpha_SNEZ) is float) or (type(alpha_SNEZ) is int)) and (alpha_SNEZ>=0 and alpha_SNEZ<=1)
		except AssertionError:
			print "Error at input 'alpha_SNEZ': must of type float or int and must take a value in [0,1]"
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print "Error: Must have run function 'timefreq_analysis'"
			return
		if scale_string is None:
			myscale=np.ceil(np.log2(self.scale[0]))
			myscale=2.**myscale
			scale_string=[myscale]
			while True:
				myscale*=2.
				if myscale>self.scale[-1]:
					break
				scale_string.append(myscale)
		mycontourf=plt.contourf(self.theta,np.log2(self.scale),np.transpose(self.scalogram),nlevels,vmin=self.minscal,vmax=self.maxscal,cmap=cmap)
		if self.shannonnyquistexclusionzone is True:
			plt.fill_between(self.theta,np.log2(self.scale[0])*np.ones(self.theta.size),np.log2(self.scalelim1_smooth),edgecolors=None,facecolor='black')
		elif self.shannonnyquistexclusionzone is False:
			plt.fill_between(self.theta,np.log2(self.scale[0])*np.ones(self.theta.size),np.log2(self.scalelim1_smooth),edgecolors=None,facecolor='black',alpha=alpha_SNEZ)
		ax = plt.gca()
		ax.tick_params(length=5, width=1, color='w')
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.ylabel("Scale"+self.t_label,fontsize=fontsize_axes)
		mytitle="Wavelet Scalogram"
		if time_string is not None:
			plt.xticks(time_string, time_string)
		plt.xticks(fontsize=fontsize_ticks)
		plt.yticks(np.log2(scale_string), scale_string, fontsize=fontsize_ticks)
		plt.xlim([self.theta[0], self.theta[-1]])
		plt.ylim([np.log2(self.scale[0]), np.log2(self.scale[-1])])
		plt.suptitle(mytitle, fontsize=fontsize_title)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		if reverse_yaxis is True:
			plt.gca().invert_yaxis()
		cbar=plt.colorbar(mycontourf)
		cbar.ax.set_ylabel("Power"+self.power_label,fontsize=fontsize_axes)
		cbar.ax.tick_params(labelsize=fontsize_ticks)
		return plt


