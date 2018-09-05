import numpy as np
from gram_schmidt import gram_schmidt
from tapering_window import tapering_window

# The time step must be constant

def WOSA_periodogram_grid_fixed(t,freq,beta,D,mywindow):

	# !!! D = nombre de points - c'est un nombre entier!
	D=int(D)
	Q=int(np.floor(float(t.size-D)/(1.-beta)/float(D)))+1
	weight=np.zeros(t.size)
	R=np.zeros((t.size,2))
	myfact=int((1.-beta)*float(D))
	W=np.zeros((t.size,2*Q))
	for k in range(Q):
		weight[:]=0.
		weight[k*myfact:k*myfact+D]=tapering_window(t[k*myfact:k*myfact+D],t[k*myfact+D-1]-t[k*myfact],mywindow)
		R[:,0]=weight*np.cos(2.*np.pi*freq*t)
		R[:,1]=weight*np.sin(2.*np.pi*freq*t)
		Y=gram_schmidt(R)
		W[:,2*k]=Y[:,0]
		W[:,2*k+1]=Y[:,1]
		#myperiodogram+=np.dot(mycos2,x)**2+np.dot(mysin2,x)**2
	#myperiodogram/=float(Q)

	return W