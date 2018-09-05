import numpy as np
import matplotlib
matplotlib.use('Agg')
import wavepal as wv
import matplotlib.pyplot as plt
from time import time
from os import path
from gen_car1 import gen_car1

here = path.abspath(path.dirname(__file__))+"/"

nt=[100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]

# 1. 95% confidence level - 10 000 samples for MCMC approach
time95_anal=np.zeros(len(nt))
time95_mcmc=np.zeros(len(nt))

percentile=np.zeros(1)
percentile[0]=95.

## 1.1 White noise
for k in range(len(nt)):
    print "wn 95%: iteration "+str(k+1)+"/"+str(len(nt))
    t=np.linspace(1.,float(nt[k]),nt[k])
    data=np.sin(2.*np.pi*t/10.)+2.*np.random.standard_normal(nt[k])
    x=wv.Wavepal(t,data)
    x.check_data()
    x.choose_trend_degree(-1)
    x.trend_vectors()
    # case 1: analytical
    time1=time()
    x.carma_params(p=0,q=0,signif_level_type="a")
    x.freq_analysis(percentile=percentile)
    time2=time()
    time95_anal[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_wn_95_anal_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
    # case 2: MCMC
    time1=time()
    x.carma_params(p=0,q=0,signif_level_type="n",nmcmc=10000)
    x.freq_analysis(percentile=percentile)
    time2=time()
    time95_mcmc[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_wn_95_mcmc_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
# Fig. of the CPU
plt.plot(nt,time95_anal,'b',label="Analytical - 95%")
plt.plot(nt,time95_anal,'b.')
plt.plot(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.plot(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_95_normal.pdf")
plt.close()
plt.semilogx(nt,time95_anal,'b',label="Analytical - 95%")
plt.semilogx(nt,time95_anal,'b.')
plt.semilogx(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.semilogx(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_95_semilogx.pdf")
plt.close()
plt.semilogy(nt,time95_anal,'b',label="Analytical - 95%")
plt.semilogy(nt,time95_anal,'b.')
plt.semilogy(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.semilogy(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_95_semilogy.pdf")
plt.close()
plt.loglog(nt,time95_anal,'b',label="Analytical - 95%")
plt.loglog(nt,time95_anal,'b.')
plt.loglog(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.loglog(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_95_loglog.pdf")
plt.close()

## 1.2 Red noise
for k in range(len(nt)):
    print "rn 95%: iteration "+str(k+1)+"/"+str(len(nt))
    t=np.linspace(1.,float(nt[k]),nt[k])
    alpha=1./10.*np.ones(1)
    sig=2.*np.ones(1)
    data=np.sin(2.*np.pi*t/10.)+gen_car1(t,alpha,sig)
    data=data[:,0]
    x=wv.Wavepal(t,data)
    x.check_data()
    x.choose_trend_degree(-1)
    x.trend_vectors()
    # case 1: analytical
    time1=time()
    x.carma_params(signif_level_type="a")
    x.freq_analysis(percentile=percentile)
    time2=time()
    time95_anal[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_rn_95_anal_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
    # case 2: MCMC
    time1=time()
    x.carma_params(signif_level_type="n",nmcmc=10000,nmcmc_carma_max=1000000)
    x.freq_analysis(percentile=percentile)
    time2=time()
    time95_mcmc[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_rn_95_mcmc_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
# Fig. of the CPU
plt.plot(nt,time95_anal,'b',label="Analytical - 95%")
plt.plot(nt,time95_anal,'b.')
plt.plot(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.plot(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_95_normal.pdf")
plt.close()
plt.semilogx(nt,time95_anal,'b',label="Analytical - 95%")
plt.semilogx(nt,time95_anal,'b.')
plt.semilogx(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.semilogx(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_95_semilogx.pdf")
plt.close()
plt.semilogy(nt,time95_anal,'b',label="Analytical - 95%")
plt.semilogy(nt,time95_anal,'b.')
plt.semilogy(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.semilogy(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_95_semilogy.pdf")
plt.close()
plt.loglog(nt,time95_anal,'b',label="Analytical - 95%")
plt.loglog(nt,time95_anal,'b.')
plt.loglog(nt,time95_mcmc,'g',label="MCMC - 95%")
plt.loglog(nt,time95_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_95_loglog.pdf")
plt.close()

# 2. 99% confidence level - 100 000 samples for MCMC approach
time99_anal=np.zeros(len(nt))
time99_mcmc=np.zeros(len(nt))

percentile=np.zeros(1)
percentile[0]=99.

## 2.1 White noise
for k in range(len(nt)):
    print "wn 99%: iteration "+str(k+1)+"/"+str(len(nt))
    t=np.linspace(1.,float(nt[k]),nt[k])
    data=np.sin(2.*np.pi*t/10.)+2.*np.random.standard_normal(nt[k])
    x=wv.Wavepal(t,data)
    x.check_data()
    x.choose_trend_degree(-1)
    x.trend_vectors()
    # case 1: analytical
    time1=time()
    x.carma_params(p=0,q=0,signif_level_type="a")
    x.freq_analysis(percentile=percentile)
    time2=time()
    time99_anal[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_wn_99_anal_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
    # case 2: MCMC
    time1=time()
    x.carma_params(p=0,q=0,signif_level_type="n",nmcmc=100000)
    x.freq_analysis(percentile=percentile)
    time2=time()
    time99_mcmc[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_wn_99_mcmc_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
# Fig. of the CPU
plt.plot(nt,time99_anal,'b',label="Analytical - 99%")
plt.plot(nt,time99_anal,'b.')
plt.plot(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.plot(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_99_normal.pdf")
plt.close()
plt.semilogx(nt,time99_anal,'b',label="Analytical - 99%")
plt.semilogx(nt,time99_anal,'b.')
plt.semilogx(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.semilogx(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_99_semilogx.pdf")
plt.close()
plt.semilogy(nt,time99_anal,'b',label="Analytical - 99%")
plt.semilogy(nt,time99_anal,'b.')
plt.semilogy(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.semilogy(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_99_semilogy.pdf")
plt.close()
plt.loglog(nt,time99_anal,'b',label="Analytical - 99%")
plt.loglog(nt,time99_anal,'b.')
plt.loglog(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.loglog(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_wn_99_loglog.pdf")
plt.close()

## 2.2 Red noise
for k in range(len(nt)):
    print "rn 99%: iteration "+str(k+1)+"/"+str(len(nt))
    t=np.linspace(1.,float(nt[k]),nt[k])
    alpha=1./10.*np.ones(1)
    sig=2.*np.ones(1)
    data=np.sin(2.*np.pi*t/10.)+gen_car1(t,alpha,sig)
    data=data[:,0]
    x=wv.Wavepal(t,data)
    x.check_data()
    x.choose_trend_degree(-1)
    x.trend_vectors()
    # case 1: analytical
    time1=time()
    x.carma_params(signif_level_type="a")
    x.freq_analysis(percentile=percentile)
    time2=time()
    time99_anal[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_rn_99_anal_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
    # case 2: MCMC
    time1=time()
    x.carma_params(signif_level_type="n",nmcmc=100000,nmcmc_carma_max=1000000)
    x.freq_analysis(percentile=percentile)
    time2=time()
    time99_mcmc[k]=time2-time1
    # Fig. of the results of the last iteration
    if k==(len(nt)-1):
        plot_periodogram=x.plot_periodogram(fontsize_legend=8)
        plot_periodogram.savefig(here+"periodogram_rn_99_mcmc_"+str(nt[-1])+".pdf")
        plot_periodogram.close()
# Fig. of the CPU
plt.plot(nt,time99_anal,'b',label="Analytical - 99%")
plt.plot(nt,time99_anal,'b.')
plt.plot(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.plot(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_99_normal.pdf")
plt.close()
plt.semilogx(nt,time99_anal,'b',label="Analytical - 99%")
plt.semilogx(nt,time99_anal,'b.')
plt.semilogx(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.semilogx(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_99_semilogx.pdf")
plt.close()
plt.semilogy(nt,time99_anal,'b',label="Analytical - 99%")
plt.semilogy(nt,time99_anal,'b.')
plt.semilogy(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.semilogy(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_99_semilogy.pdf")
plt.close()
plt.loglog(nt,time99_anal,'b',label="Analytical - 99%")
plt.loglog(nt,time99_anal,'b.')
plt.loglog(nt,time99_mcmc,'g',label="MCMC - 99%")
plt.loglog(nt,time99_mcmc,'g.')
plt.title("Computing times")
plt.xlabel("# data points")
plt.ylabel("Computing time (s)")
plt.legend(fancybox=True,fontsize='small')
plt.savefig(here+"CPU_rn_99_loglog.pdf")
plt.close()
