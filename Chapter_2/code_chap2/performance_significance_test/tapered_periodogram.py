import numpy as np
from gram_schmidt import gram_schmidt
from tapering_window import tapering_window

def tapered_periodogram(t,x,freq,mywindow):

	"""Returns the Lomb-Scargle periodogram of the data x(t) at frequency freq."""

	R=np.zeros((t.size,2))
	taper=tapering_window(t,t[-1]-t[0],mywindow)
	R[:,0]=taper*np.cos(2.*np.pi*freq*t)
	R[:,1]=taper*np.sin(2.*np.pi*freq*t)
	Q=gram_schmidt(R)
	myperiodogram=np.dot(Q[:,0],x)**2+np.dot(Q[:,1],x)**2

	return myperiodogram