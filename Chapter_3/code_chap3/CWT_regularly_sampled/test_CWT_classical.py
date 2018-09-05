import numpy as np
from CWT_classical import CWT_Morlet_regularly_sampled
import wavepal as wv
from gen_car1 import gen_car1
import matplotlib.pyplot as plt

# Test 1
t=np.linspace(1.,1000.,453)
y=np.sqrt(5.)*np.sin(2.*np.pi*t/20.)
myplot=CWT_Morlet_regularly_sampled(t,y)
myplot.show()

# Test 2
t=np.linspace(5.,5000.,1000)
y=np.sin(2.*np.pi*t/20.)+np.sin(2.*np.pi*t/400.)+np.sqrt(0.5)*np.sin(2.*np.pi*t/100.)
myplot=CWT_Morlet_regularly_sampled(t,y,Tmin=15.)
myplot.show()

## Test 3 - Pure noise
t=np.linspace(1.,5000.,383)
y=gen_car1(t,0.45*np.ones(1),1.7689*np.ones(1))
y=y[:,0]+np.sin(2.*np.pi*t/300.)
x=wv.Wavepal(t,y)
x.check_data()
x.choose_trend_degree(-1)
x.trend_vectors()
x.carma_params(signif_level_type="a",nmcmc_carma_max=100000)
sig_wn=x.sigwn_unique
alpha_rn=x.alpha_unique
x.timefreq_analysis(w0=8.7,smoothing_coeff=0.,shannonnyquistexclusionzone=True,weighted_CWT=True,computes_global_scalogram=False)
plot_scalogram=x.plot_scalogram(linewidth_cl=4,with_global_scalogram=False)
plot_scalogram.savefig("wavepal_scal.pdf")
plot_scalogram.show()
myplot=CWT_Morlet_regularly_sampled(t,y,alpha_rn=alpha_rn,sigma_wn=sig_wn,linewidth_cl=4,w0=8.7)
myplot.show()

# Test 4 - with wavepal for the background noise
data=np.genfromtxt("GISP2d18O_splineinterp_30yr.txt")
#myt=data[:,0]
#mydata=data[:,1]
myt=data[:500,0]
mydata=data[:500,1]
mydata-=np.mean(mydata)
plt.plot(myt,mydata)
plt.show()
x=wv.Wavepal(myt, mydata, "Age", t_units="a")
x.check_data()
x.choose_trend_degree(-1)
x.trend_vectors()
x.carma_params(signif_level_type="a",nmcmc_carma_max=100000)
sig_wn=x.sigwn_unique
alpha_rn=x.alpha_unique
x.timefreq_analysis(w0=15.,smoothing_coeff=0.,shannonnyquistexclusionzone=True,computes_global_scalogram=False)
plot_scalogram=x.plot_scalogram(linewidth_cl=4,with_global_scalogram=False)
plot_scalogram.savefig("wavepal_scal.pdf")
plot_scalogram.show()
myplot=CWT_Morlet_regularly_sampled(myt,mydata,w0=15.,alpha_rn=alpha_rn,sigma_wn=sig_wn,t_axis_label="Age",t_units="a",linewidth_cl=4)
myplot.show()
