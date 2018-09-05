import numpy as np
from scipy.stats import chi2 as chi2distr
from scipy.optimize import brenth

def fct_bound_percentiles(lambdas,proba):

	trace=np.sum(lambdas)
	norm_frob_sq=np.sum(lambdas**2)
	norm_inf=np.amax(lambdas)
	# 2-moment approx
	M=trace**2/norm_frob_sq
	g=norm_frob_sq/trace
	percentile_2moment_centered=g*chi2distr.ppf(proba,M)-trace
	print g*chi2distr.ppf(proba,M)
	# bound on the true value
	myfun=lambda t: t*min(t**2/4./norm_frob_sq,t/2./norm_inf)-(np.log(2)-np.log(proba))/2.
	percentile,rootresults=brenth(myfun,percentile_2moment_centered/1000.,percentile_2moment_centered*1000.,maxiter=1000,full_output=True,disp=False)
	if rootresults.converged==False:
		print "Error in fct_bound_percentiles"
		return
	else:
		return percentile+trace