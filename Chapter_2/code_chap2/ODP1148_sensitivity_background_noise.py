import numpy as np
import matplotlib.pyplot as plt
import copy
import wavepal as wv

mypath="../figures/ODP1148_example/"

# Read the data
data=np.genfromtxt("data/ODP1148-BF-18O.txt")
myt=data[:,0]
mydata=data[:,1]

# Initialize the class called Wavepal (which is a class of the package wv)
x=wv.Wavepal(myt, mydata, "Age", "$\delta{}^{18}O$", t_units="ka", mydata_units="permil")

# Check the data set
x.check_data()

# Choose the degree of the polynomial for the subsequent analyses
x.choose_trend_degree(7)

# Compute some variables related to the trend
x.trend_vectors()

# Sensitivity analysis for the background noise (MCMC CL)
# -------------------------------------------------------
nmcmc=10000
# White noise
x.carma_params(p=0,q=0,path_to_figure_folder=mypath+"sensitivity_order_carma/wn/",nmcmc=nmcmc,signif_level_type="n",make_carma_fig=False,nbins=15,dpi=400)
x.freq_analysis(freqstep=0.0001,D=600.,mywindow=8)
cl00=x.periodogram_cl_mcmc[:,0]
# Red noise
x.carma_params(p=1,q=0,path_to_figure_folder=mypath+"sensitivity_order_carma/rn/",nmcmc=nmcmc,signif_level_type="n",make_carma_fig=False,nbins=15,dpi=400)
x.freq_analysis(freqstep=0.0001,D=600.,mywindow=8)
cl10=x.periodogram_cl_mcmc[:,0]
# CARMA(2,0)
x.carma_params(p=2,q=0,path_to_figure_folder=mypath+"sensitivity_order_carma/carma20/",nmcmc=nmcmc,signif_level_type="n",make_carma_fig=False,nbins=15,dpi=400)
x.freq_analysis(freqstep=0.0001,D=600.,mywindow=8)
cl20=x.periodogram_cl_mcmc[:,0]
# CARMA(2,1)
x.carma_params(p=2,q=1,path_to_figure_folder=mypath+"sensitivity_order_carma/carma21/",nmcmc=nmcmc,signif_level_type="n",make_carma_fig=False,nbins=15,dpi=400)
x.freq_analysis(freqstep=0.0001,D=600.,mywindow=8)
cl21=x.periodogram_cl_mcmc[:,0]
# Figures
plt.plot(x.freq,x.periodogram,label="Data periodogram")
plt.plot(x.freq,cl00,label="(p,q)=(0,0) - MCMC CL at 95% (from stoch. param.)")
plt.plot(x.freq,cl10,label="(p,q)=(1,0) - MCMC CL at 95% (from stoch. param.)")
plt.plot(x.freq,cl20,label="(p,q)=(2,0) - MCMC CL at 95% (from stoch. param.)")
plt.plot(x.freq,cl21,label="(p,q)=(2,1) - MCMC CL at 95% (from stoch. param.)")
plt.xlabel("Frequency"+" (ka${}^{-1}$)")
plt.ylabel("Power"+" (permil${}^{2}$)")
plt.suptitle("WOSA periodogram and Confidence levels")
plt.legend(fancybox=True,fontsize=7,bbox_to_anchor=(1.1, 1.05))
plt.savefig(mypath+"/sensitivity_order_carma/peridogram.pdf")
plt.close()
