import numpy as np
from gram_schmidt import gram_schmidt

def WOSA_matrix(t,freq,beta,D):

	"""Returns the Lomb-Scargle+WOSA periodogram of the data x(t) at frequency freq."""

	Q=int(np.floor((t.size-D)/(1.-beta)/float(D)))+1
	weight=np.zeros(t.size)
	R=np.zeros((t.size,2))
	M2=np.zeros((t.size,2*Q))
	myfact=int((1.-beta)*float(D))
	for k in range(Q):
		weight[:]=0.
		weight[(k-1)*myfact:(k-1)*myfact+D]=1.
		R[:,0]=weight*np.cos(2.*np.pi*freq*t)
		R[:,1]=weight*np.sin(2.*np.pi*freq*t)
		Y=gram_schmidt(R)
		M2[:,2*k]=Y[:,0]
		M2[:,2*k+1]=Y[:,1]
	M2/=np.sqrt(float(Q))
	
	return M2