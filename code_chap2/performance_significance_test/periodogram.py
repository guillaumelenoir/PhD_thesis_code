import numpy as np
from gram_schmidt import gram_schmidt

def periodogram(t,x,freq):

	"""Returns the Lomb-Scargle periodogram of the data x(t) at frequency freq."""

	R=np.zeros((t.size,2))
	R[:,0]=np.cos(2.*np.pi*freq*t)
	R[:,1]=np.sin(2.*np.pi*freq*t)
	Q=gram_schmidt(R)
	myperiodogram=np.dot(Q[:,0],x)**2+np.dot(Q[:,1],x)**2

	return myperiodogram