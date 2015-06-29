import os, sys 
import pandas as pd 
import numpy as np
import readhatlc as rlc
from lsp import *
from math import *
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import defaults

data_dir = "/Users/jah5/Documents/Fall2014_Gaspar/data"
time_col = 'BJD'
grpsize = 50
dt_log = np.dtype([('hatid', 'S15'), ('flag', 'S20')])
API_KEY = {
    'lcdirect' : 'ZjNmZjQ0NzY4MTQxNzQ0Zjk1OTdlNzY1MTAxOTY1YTQyNDNlMzZlZmE2MWE3M2E3YTY0OWE1MDM5ZDU5NmRjYQ'
}


variance = lambda x, mu : sqrt(np.mean(np.power(x - mu,2)))
rand = np.random.random


def data_kind(arr):
    if nums(arr[0]) is None:
        return "categorical"
    else:
        return "numerical"
def encode(arr):
    cats = [ c for c in set(arr)]
    codes = np.arange(0,len(cats))
    trans = {}
    for i in range(0,len(cats)):
        trans[cats[i]] = codes[i]
    return cats, np.array([ trans[c] for c in arr ])
def phase_fold(t,y,T,phase_offset=0, bin=defaults.settings['bin-phase'], 
    nbins=10, remove_outliers=defaults.settings['remove-outliers'], 
    outlier_limit=defaults.settings['outlier-sigma']):
    """ Given a period (T), phase-folds the (t,y) data and bins the results.
        Returns: phase_bin_centers, binned_magnitudes, rms_magnitudes
    """

    # Initial time
    t0 = t[0]

    # Returns phase
    def phase(time):
        DT = (time - t0)/T + phase_offset
        DT -= int(DT)
        return DT

    # Make a record array, then sort by phase
    pf = np.zeros(len(t), dtype=np.dtype([('phase', np.float_) , ('mag', np.float_)]))
    for i,tv in enumerate(t):
        ph = phase(tv)
        pf[i]['phase'] = ph 
        pf[i]['mag'] = y[i]
    pf = np.sort(pf,order='phase')
    
    # If we aren't binning, return the results
    if not bin: return pf['phase'],pf['mag'], None

    dt = 1./nbins

    phases = np.zeros(nbins)
    

    # Set central bin values
    for i in range(0,nbins):
        phases[i] = (i + 0.5)*dt
    outliers = None 
    while outliers != 0:
        binned_values = np.zeros(nbins)
        errs = np.zeros(nbins)
        nvals = np.zeros(nbins)
        # Now bin!
        for v in pf:
            index = int(v['phase']/dt)
            binned_values[index] += v['mag']
            nvals[index] += 1

        # Divide by number of measurements per bin
        for i in range(0,nbins):
            if nvals[i] != 0:
                binned_values[i]/=nvals[i]

        # Now calculate the RMS error in each bin
        for v in pf:
            index = int(v['phase']/dt)
            errs[index] += pow(v['mag'] - binned_values[index],2)

        # Normalize RMS error
        for i in range(0,nbins):
            if nvals[i] != 0:
                errs[i]=sqrt(errs[i]/nvals[i])
        if remove_outliers:
            # Remove outliers
            pf_new = []
            outliers = 0
            for v in pf:
                index = int(v['phase']/dt)
                if errs[index] == 0: continue
                if abs(v['mag'] - binned_values[index])/errs[index] < outlier_limit:
                    pf_new.append(v)
                else: outliers += 1
            pf = pf_new
        else: outliers = 0
    return phases, binned_values, errs
def fetch_lcs(hatids):
    i = 0
    lh = len(hatids)
    IDgrp = []
    def fetch(grp):
        idlist=""
        for i,ID in enumerate(grp):
            idlist = "%s%s"%(idlist,ID)
            if i < len(grp) - 1: idlist = "%s,"%(idlist)
        url = 'https://hatsurveys.org/lightcurves/lc/direct?hatid=%s&apikey=%s'%(idlist, API_KEY['lcdirect'])
        fetch_command = "curl -J -O '%s'"%(url)
        unzip_command = "unzip *zip; rm *zip"
        os.system(fetch_command)
        os.system(unzip_command)
        os.system("mv *hatlc*gz %s"%(data_dir))
    j=0
    while j < len(hatids):
        if len(IDgrp) == grpsize:
            fetch(IDgrp)
            IDgrp = []
        if not os.path.exists("%s/%s-hatlc.csv.gz"%(data_dir,hatids[j])): 
            IDgrp.append(hatids[j])
        j += 1
    if len(IDgrp) > 0: fetch(IDgrp)

def get_lombscarg(Tt, x, nfrs = defaults.settings['nfrqs']):
    t = Tt - Tt[0]
    
    # LSP code does NOT use angular frequencies#
    freqs, power, nout, jmax, prob = fasper(t,x, 6., 6.) 

    periods = np.power(freqs,-1)
    return periods[::-1], power[::-1], periods[jmax]
def lsp_raw(t,x, ofac=6, hifac=6, MACC =2):
    freqs, lsppowers, ndim, imax, prob = fasper(t,x,ofac,hifac, MACC=MACC)
    return freqs, lsppowers

#def window_lsp(t, ofac=6, hifac=6, MACC=2):
#    freqs, lsppowers, ndim, imax, prob = fasper(t,np.ones(len(t)),ofac,hifac, MACC=MACC)
#    return freqs, lsppowers
def nums(st):
    try:
        return float(st)
    except ValueError:
        return None
def prune_nans(t,x):
    tvals = []
    xvals = []
    for i in range(0,len(t)):
        if not np.isnan(x[i]):
            tvals.append(t[i])
            xvals.append(x[i])
    return np.array(tvals), np.array(xvals)
def is_peak(powers, i, imin, imax, per, peak_pers):
    p = powers[i]
    for I in range(imin,imax+1):
        if I == i: continue
        if powers[I] > p: return False
    for pper in peak_pers:
        if abs(per - pper)/pper < 0.05: return False
    return True
def find_n_peaks(periods, powers, n_peaks, dn=defaults.settings['peaks-dn']):
    peak_periods = []
    peak_powers = []
    ls0 = np.zeros(len(periods), dtype=np.dtype([('periods', np.float_), ('powers', np.float_)]))
    ls = np.zeros(len(periods), dtype=np.dtype([('periods', np.float_), ('powers', np.float_), ('indices', np.int_)]))
    for i in range(0,len(periods)):
        ls['periods'][i] = periods[i] 
        ls['powers'][i] = powers[i]
        ls['indices'][i] = i
        ls0['periods'][i] = periods[i]
        ls0['powers'][i]  = powers[i]

    ls = np.sort(ls, order='powers')[::-1]
    for L in ls:
        if len(peak_periods) == n_peaks: break
        I = L['indices']
        imin = max([ I - dn, 0 ] )
        imax = min([ I + dn, len(ls) - 1])
        if is_peak(ls0['powers'], I, imin, imax, L['periods'], peak_periods): 
            peak_periods.append(L['periods'])
            peak_powers.append(L['powers'])
    
    return np.array(peak_periods), np.array(peak_powers)

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
        if avg_rms[filt] > 0:
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
    lsp_null = np.array([ np.mean(lsps[i,:])                for i in range(len(freqs)) ])
    lsp_var  = np.array([ variance(lsps[i,:],lsp_null[i])   for i in range(len(freqs)) ])
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


