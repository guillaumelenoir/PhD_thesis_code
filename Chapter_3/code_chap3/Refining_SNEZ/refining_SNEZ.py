import numpy as np
import matplotlib.pyplot as plt
import wavepal_pack as wv	# modified version of wavepal, without the cone of influence for 'plot_scalogram'

mypath="../../figures/Refining_SNEZ/"

# Read the data
myt=np.arange(0.,204.,4.)
myt=np.append(myt,np.arange(203.,400.,3.))
myt=np.append(myt,np.arange(402.,602.,2.))
print myt.size

mydata=np.sin(2.*np.pi*myt/10.)

# Initialize the class called Wavepal (which is a class of the package wv)
x=wv.Wavepal(myt, mydata, "Time", "")

# Check the data set
x.check_data()

# plot the time series
plt.plot(x.t,x.mydata,'r')
plt.plot(x.t,x.mydata,'b.')
plt.xlabel("Time",fontsize=12)
plt.title("Time series",fontsize=14)
plt.savefig(mypath+"timeseries.pdf")
plt.close()

# plot time step
plt_timestep=x.plot_timestep(hist=False)
plt_timestep.savefig(mypath+"timestep.pdf")
plt_timestep.close()

# Choose the degree of the polynomial for the subsequent analyses
x.choose_trend_degree(-1)

# Compute some variables related to the trend
x.trend_vectors()

# Compute and plot the scalogram
theta=np.linspace(myt[0],myt[-1],300)
x.timefreq_analysis(theta=theta,activate_perlim2=True,computes_global_scalogram=False,deltaj=0.01,smoothing_coeff=0.,computes_amplitude=True)
plot_scalogram=x.plot_scalogram(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,with_global_scalogram=False)
plot_scalogram.suptitle("Wavelet scalogram - With refinement of the SNEZ",fontsize=14)
plot_scalogram.savefig(mypath+"Scalogram_with_perlim2.pdf")
plot_scalogram.close()
plot_scalogram=x.plot_cwtamplitude(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12)
plot_scalogram.suptitle("Wavelet amplitude - With refinement of the SNEZ",fontsize=14)
plot_scalogram.savefig(mypath+"Amplitude_with_perlim2.pdf")
plot_scalogram.close()
x.timefreq_analysis(theta=theta,activate_perlim2=False,computes_global_scalogram=False,deltaj=0.01,smoothing_coeff=0.,computes_amplitude=True)
plot_scalogram=x.plot_scalogram(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,with_global_scalogram=False)
plot_scalogram.suptitle("Wavelet scalogram - No refinement of the SNEZ",fontsize=14)
plot_scalogram.savefig(mypath+"Scalogram_no_perlim2.pdf")
plot_scalogram.close()
plot_scalogram=x.plot_cwtamplitude(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12)
plot_scalogram.suptitle("Wavelet amplitude - No refinement of the SNEZ",fontsize=14)
plot_scalogram.savefig(mypath+"Amplitude_no_perlim2.pdf")
plot_scalogram.close()
