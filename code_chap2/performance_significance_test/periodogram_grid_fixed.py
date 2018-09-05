import numpy as np
from gram_schmidt import gram_schmidt

def periodogram_grid_fixed(t,freq):

	R=np.zeros((t.size,2))
	R[:,0]=np.cos(2.*np.pi*freq*t)
	R[:,1]=np.sin(2.*np.pi*freq*t)
	Q=gram_schmidt(R)
	#mycos2=Q[:,0]
	#mysin2=Q[:,1]
	#myperiodogram=np.dot(mycos2,x)**2+np.dot(mysin2,x)**2

	return Q