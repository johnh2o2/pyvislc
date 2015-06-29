import numpy as np
from lsp import *
from math import *
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Since we're computing LSP powers *relative* to 
# their NULL values, hifac should in principle be as high
# as you want
def lsp_raw(t,x, ofac=6, hifac=2, MACC =2):
	freqs, lsppowers, ndim, imax, prob = fasper(t,x,ofac,hifac, MACC=MACC)
	return freqs, lsppowers

variance = lambda x, mu : sqrt(np.mean(np.power(x - mu,2)))
rand = np.random.random

def filter_normalization(t,x,filts):
	unique_filts = np.unique(filts)
	avg_mags = {}
	avg_rms = {}
	filt_data = {}
	for filt in unique_filts: filt_data[filt] = x[np.where(filts == filt)]
	max_nobs = None
	most_observed_filt = None
	for filt in unique_filts:
		if max_nobs is None:
			most_observed_filt = filt 
			max_nobs = len(filt_data[filt])
		elif len(filt_data[filt]) > max_nobs:
			max_nobs = len(filt_data[filt])
			most_observed_filt = filt
	for filt in unique_filts:
		inds = np.where(filts == filt)
		T, X = t[inds], x[inds]
	
		avg_mags[filt] = np.mean(X)
		avg_rms[filt] = variance(X,avg_mags[filt])
	for filt in unique_filts:
		inds = np.where(filts == filt)
		x[inds]-=avg_mags[filt]
		x[inds]*=(avg_rms[most_observed_filt]/avg_rms[filt])
		x[inds]+=avg_mags[most_observed_filt]
	return t,x


def get_cdf(x, nbins=1000):
	ext_fac = 0.1
	min_val, max_val = (1 - ext_fac)*min(x), (1+ext_fac)*max(x)
	assert(max_val > min_val)
	np.sort(x)
	xvals = np.linspace(min_val, max_val,num=nbins)
	dx = xvals[1] - xvals[0]
	cdf = np.zeros(len(xvals))
	for i in range(len(x)):
		index = int((x[i] - min_val)/dx)
		cdf[index:] += 1
	return xvals, cdf/len(x)

def get_inv_cdf(xvals, cdf):
	inds = []

	for i in range(len(cdf)):
		if cdf[i] > 0.0 and cdf[i] < 1.0:
			inds.append(i)
	
	inds.extend(  max(np.where(cdf == 0.0)) )
	inds.extend(  min(np.where(cdf == 1.0)) )
	X,CDF = xvals[inds], cdf[inds]
	return interp1d(CDF,X)

def get_Ndraws(N,invcdf):
	ndraws = np.zeros(N)
	for i in range(N):
		ndraws[i] = invcdf(rand())

	return ndraws

def scaled_lsp(t, x, NMC=10):
	assert(len(x) == len(t))
	print "Getting raw lsp"
	freqs, lsp_observed = lsp_raw(t,x)
	print "Getting CDF of magnitudes"
	#mags, cdf = get_cdf(x)
	#plt.plot(mags,cdf)
	#plt.show()

	print "Getting inverse cdf"
	invcdf = get_inv_cdf(mags, cdf)

	lsp_null = np.zeros(len(freqs))
	lsps = np.zeros(shape=(len(freqs), NMC))
	print "Using MC to figure out the null LSP"
	for i in range(NMC):
		print "   MC ",i
		MCLC = get_Ndraws(len(x), invcdf)
		#plt.scatter(t, MCLC, facecolor='b',alpha=0.5)
		#plt.scatter(t,x, facecolor='k', alpha=0.3)
		#plt.show()
		print "      doing lsp"
		frqs,MCLSP = lsp_raw(t, MCLC)
		print "      done."
		lsps[:,i] = MCLSP[:]
	print "Setting values"
	lsp_null = np.array([ np.mean(lsps[i,:]) 				for i in range(len(freqs)) ])
	lsp_var	 = np.array([ variance(lsps[i,:],lsp_null[i]) 	for i in range(len(freqs)) ])
	# maybe it will eventually make sense to do a non-parametric 
	# characterization of the LSP power distribution for a given frequency

	# Replace any zeros in the variance array if there are any...

	min_nonzero = min([ v for v in lsp_var if v > 0 ] )
	for i in range(len(lsp_var)):
		if lsp_var[i] == 0: lsp_var[i] = min_nonzero

	#plt.plot(freqs, lsp_null, color='k')
	#plt.plot(freqs, lsp_observed, color='g')
	#plt.show()

	lsp_statistic = np.divide(lsp_observed - lsp_null, lsp_var)

	return freqs, lsp_statistic
