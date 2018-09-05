

# Obsolete - Not used


import numpy as np
from gram_schmidt import gram_schmidt
import time

def WOSA_periodogram(t,x,freq,beta,D):

	"""Returns the Lomb-Scargle+WOSA periodogram of the data x(t) at frequency freq."""

	Q=int(np.floor((t.size-D)/(1.-beta)/float(D)))+1
	#weight=np.zeros(t.size)
	R=np.zeros((t.size,2))
	myperiodogram=0.
	myfact=int((1.-beta)*float(D))
	for k in range(Q):
		#weight[:]=0.
		#weight[k*myfact:k*myfact+D]=1.
		R[:,0]=weight*np.cos(2.*np.pi*freq*t)
		R[:,1]=weight*np.sin(2.*np.pi*freq*t)
		Y=gram_schmidt(R)
		mycos2=Y[:,0]
		mysin2=Y[:,1]
		myperiodogram+=np.dot(mycos2,x)**2+np.dot(mysin2,x)**2
	myperiodogram/=float(Q)

	return myperiodogram