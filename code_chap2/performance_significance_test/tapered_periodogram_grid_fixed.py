import numpy as np
from gram_schmidt import gram_schmidt
from tapering_window import tapering_window

def tapered_periodogram_grid_fixed(t,freq,mywindow):

	R=np.zeros((t.size,2))
	taper=tapering_window(t,t[-1]-t[0],mywindow)
	R[:,0]=taper*np.cos(2.*np.pi*freq*t)
	R[:,1]=taper*np.sin(2.*np.pi*freq*t)
	Q=gram_schmidt(R)
	#mycos2=Q[:,0]
	#mysin2=Q[:,1]
	#myperiodogram=np.dot(mycos2,x)**2+np.dot(mysin2,x)**2

	return Q