import numpy as np
from tapering_window import tapering_window

# computes the average WEIGHTED central time-step (the weight = the taper window)

def dt_av_central(t,D,mywindow):
	
	""" dt_av returns the average weighted central time step on a given WOSA segment. The weight is the taper window and is computed at t.
		Inputs:
		- t [1-dim numpy array of floats]: the times.
		- D [float]: the temporal length of the WOSA segment.
		- mywindow [int]: window choice for the windowing of the WOSA segment. See tapering_window.py for more details.
		Outputs:
		- dt [float]: the average weigthed central time step.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	n=t.size
	G=tapering_window(t,D,mywindow)
	dt_weighted=(t[1]-t[0])*G[0]
	for k in range(1,n-1):
		dt_weighted+=(t[k+1]-t[k-1])*G[k]/2.
	dt_weighted+=(t[n-1]-t[n-2])*G[n-1]
	dt_weighted/=np.sum(G)

	return dt_weighted