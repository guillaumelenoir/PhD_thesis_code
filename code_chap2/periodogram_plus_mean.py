# figures for periodogram + mean

import numpy as np
import matplotlib.pyplot as plt

mypath="../figures/"

# First figure
t=[3., 4., 6., 7., 18., 19., 32., 33., 34., 35., 40., 41., 51., 55., 56., 57., 60., 68., 69., 72., 77., 78., 79., 80., 81., 90., 100.]
t_irr=np.zeros(len(t))
t_irr[:]=t[:]
t_reg=np.arange(0.,100.01,0.01)
myfun=lambda x: np.sin(2.*np.pi*x/25.)+1.
#plt.plot(t_reg,myfun(t_reg),"b--")
#plt.plot(t_irr,myfun(t_irr),"ro",)
#plt.plot(t_reg,np.ones(t_reg.size),'b')
#plt.plot(t_reg,np.mean(myfun(t_irr))*np.ones(t_reg.size),'r')
#plt.xlabel("Time")
#plt.axis([0., 100., -0.2, 2.2])
#figname=mypath+"periodogram_and_mean_irregular_sampling.pdf"
#plt.savefig(figname)
#plt.close()

# Second figure
t_regu=np.arange(0.,102.,2.)
myfunc=lambda x: np.sin(2.*np.pi*x/80.)+1.
#plt.plot(t_regu,myfunc(t_regu),"b--")
#plt.plot(t_regu,myfunc(t_regu),"r.")
#plt.plot(t_regu,np.ones(t_regu.size),'b')
#plt.plot(t_regu,np.mean(myfunc(t_regu))*np.ones(t_regu.size),'r')
#plt.xlabel("Time")
#plt.axis([0., 100., -0.2, 2.2])
#figname=mypath+"periodogram_and_mean_regular_sampling.pdf"
#plt.savefig(figname)
#plt.close()

# Both together
mainframe, axarr = plt.subplots(2,1)
aa=axarr[0]
aa.plot(t_reg,myfun(t_reg),"b--")
aa.plot(t_irr,myfun(t_irr),"ro",)
aa.plot(t_reg,np.ones(t_reg.size),'b')
aa.plot(t_reg,np.mean(myfun(t_irr))*np.ones(t_reg.size),'r')
aa.axis([0., 100., -0.2, 2.2])
aa.text(2.,0.,'(a)')
aa=axarr[1]
aa.plot(t_regu,myfunc(t_regu),"b--")
aa.plot(t_regu,myfunc(t_regu),"r.")
aa.plot(t_regu,np.ones(t_regu.size),'b')
aa.plot(t_regu,np.mean(myfunc(t_regu))*np.ones(t_regu.size),'r')
aa.set_xlabel("Time")
aa.axis([0., 100., -0.2, 2.2])
aa.text(2.,0.,'(b)')
figname=mypath+"periodogram_and_mean.pdf"
plt.savefig(figname)
plt.close()



