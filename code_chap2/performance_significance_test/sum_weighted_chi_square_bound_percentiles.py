import numpy as np
import matplotlib.pyplot as plt
from fct_bound_percentiles import fct_bound_percentiles

# params
lambdavec=np.zeros(6); lambdavec[0]=2.5; lambdavec[1]=2.5; lambdavec[2]=5.; lambdavec[3]=5.; lambdavec[4]=7.5; lambdavec[5]=7.5
proba=0.95

percentile=fct_bound_percentiles(lambdavec,proba)
print percentile