import numpy as np
import mypack as wv    # Simplified version of Wavepal, with the two cases for the local time step defining the Shannon-Nyquist exclusion zone
import matplotlib.pyplot as plt

mypath="../../figures/SNEZ/"

# Read the data
myt=np.genfromtxt("data/Ti_8_15_ages.txt")
myt=myt[4000:7000]
myt-=myt[0]
mydata=np.sin(2.*np.pi*myt/0.010)

# Initialize the class called Wavepal (which is a class of the package wv)
x=wv.Wavepal(myt, mydata, "Time", "")

# Check the data set
x.check_data()
print x.t.size

# plot the time series
plt.plot(x.t,x.mydata,'b.')
plt.xlabel("Time",fontsize=12)
plt.title("Time series",fontsize=14)
plt.savefig(mypath+"timeseries.pdf")
plt.close()

# plot time step
plt_timestep=x.plot_timestep(hist=False,log_yaxis=True)
plt_timestep.title("Time step in function of time",fontsize=14)
plt_timestep.savefig(mypath+"timestep.pdf")
plt_timestep.close()

# Choose the degree of the polynomial for the subsequent analyses
x.choose_trend_degree(-1)

# Compute some variables related to the trend
x.trend_vectors()

# Compute and plot the scalogram
theta=np.linspace(myt[0],myt[-1],1000)
scale_string=[0.0001, 0.0002, 0.0004, 0.0008, 0.0016, 0.0032, 0.0064, 0.0128, 0.0256, 0.0512]
x.timefreq_analysis(theta=theta,shannonnyquistexclusionzone_type="G",shannonnyquistexclusionzone=False)
plot_scalogram=x.plot_scalogram(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,scale_string=scale_string,alpha_SNEZ=0)
plot_scalogram.suptitle("Wavelet scalogram - No SNEZ",fontsize=14)
plot_scalogram.savefig(mypath+"w0_5,5_noSNEZ.pdf")
plot_scalogram.close()
x.timefreq_analysis(theta=theta,shannonnyquistexclusionzone_type="G",shannonnyquistexclusionzone=True)
plot_scalogram=x.plot_scalogram(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,scale_string=scale_string)
plot_scalogram.suptitle("Wavelet scalogram - SNEZ with matrix G",fontsize=14)
plot_scalogram.savefig(mypath+"G_w0_5,5.pdf")
plot_scalogram.close()
x.timefreq_analysis(theta=theta,shannonnyquistexclusionzone_type="H",shannonnyquistexclusionzone=True)
plot_scalogram=x.plot_scalogram(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,scale_string=scale_string)
plot_scalogram.suptitle("Wavelet scalogram - SNEZ with matrix H",fontsize=14)
plot_scalogram.savefig(mypath+"H_w0_5,5.pdf")
plot_scalogram.close()
x.timefreq_analysis(theta=theta,shannonnyquistexclusionzone_type="GH",shannonnyquistexclusionzone=True)
plot_scalogram=x.plot_scalogram(fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,scale_string=scale_string)
plot_scalogram.suptitle("Wavelet scalogram - SNEZ with matrices G and H",fontsize=14)
plot_scalogram.savefig(mypath+"GH_w0_5,5.pdf")
plot_scalogram.close()
