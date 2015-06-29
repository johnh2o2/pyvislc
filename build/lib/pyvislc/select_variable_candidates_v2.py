import pandas as pd
from math import *
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import readhatlc as rhlc
from utils import *

# Steps: 
# 1) select a subset of the full dataset
# 2) calculate outlier-subtracted rms for magnitude bins
# 3) determine which sources have high RMS
# 4) print list of HAT-ID's into a file.
overwrite = False
MAG_COL = "TF1"
MAGNITUDE = 'vmag'
power_cut = 100
rmsfile = "stats_rms.dat"
NLCS = 10000
OUTLIER_SIGMA = 5

NBINS = 15

num_peaks = 3
rms = {}
magnitudes = {}
lsp_peaks = {}


print "Loading list of HAT sources..."
sources = pd.read_csv("/Users/jah5/Documents/Fall2014_Gaspar/data/hat-field-219-object-list.csv")
sources = sources.iloc[:-1]


print "Now selecting %d random HAT sources"%(NLCS)
indices = np.arange(0,len(sources))
#np.random.shuffle(indices)
selected_indices = indices[0:NLCS]

print "Setting magnitudes"
hatids = [ hid for hid in sources.iloc[selected_indices]['#hat_id'] ]
for i in selected_indices:
    s = sources.iloc[i]
    hid = s['#hat_id']
    magnitudes[hid] = s[MAGNITUDE]

# FETCH lightcurves from the server
print "Fetching lightcurves ... "
fetch_lcs(hatids)

sdata = None
pk_data = None
if not overwrite:
    print "Calculating stats..."
    dt_stats = np.dtype([('id','S15'),('rms',np.float_)])
    
    if os.path.exists(rmsfile):
        print "Found stats_rms.dat, reading in data"
        f = open(rmsfile,'r')
        sdata = np.loadtxt(rmsfile,dtype=dt_stats)
        for i,hid in enumerate(sdata['id']):
            if hid not in hatids: continue
            rms[hid] = sdata[i]['rms']
        f.close()

    # Get LSP values
    print "Calculating LSP values..."
    pk_stats = np.dtype([('id','S15'),('period',np.float_),('power',np.float_)])

    
    for pkno in range(num_peaks):
        peaks_fname = "stats_peaks_%d.dat"%(pkno)
        if os.path.exists(peaks_fname):
            print "Found %s, reading in data"%(peaks_fname)
            f = open(peaks_fname,'r')
            pk_data = np.loadtxt(peaks_fname,dtype=pk_stats)
            #print pk_data
            for i,hid in enumerate(pk_data['id']):
                if hid not in hatids: continue
                if hid not in lsp_peaks: lsp_peaks[hid] = {}
                lsp_peaks[hid][pkno] = (pk_data['period'][i], pk_data['power'][i])
            f.close()



def bin(mag, rm, nbins=NBINS, remove_outliers=True, outlier_dist=OUTLIER_SIGMA):
    sources = np.zeros(len(rm), dtype=np.dtype([('mag', np.float_), ('rms', np.float_)]))
    for i in range(0,len(mag)):
        sources[i]['mag'] = mag[i]
        sources[i]['rms'] = rm[i]
    #print len(sources), len(mag)
    bin_edges = []
    means = []
    stddev = []
    sources = np.sort(sources,order='mag')
    eps = 10E-5
    i=0
    j=None
    bin_edges.append(sources[0]['mag'] - eps)
    while i < len(sources) - 1:
        
        j = min([i + len(sources)/nbins, len(sources) - 1])
        bin_edges.append(sources[j]['mag'] + eps)

        bin_sources = sources[i:j+1]
        outliers = None
        while outliers != 0:
            #print "bin %d, nsrc = %d"%( i, len(bin_sources))
            RMS = [ b['rms'] for b in bin_sources ] 
            mu = np.mean(RMS )
            err = sqrt(np.mean(np.power(RMS - mu,2)))
            bs_new = []
            outliers = 0
            for b in bin_sources:
                if abs(b['rms'] - mu)/err < outlier_dist: bs_new.append(b)
                else: outliers += 1
            bin_sources = bs_new
        means.append(mu)
        stddev.append(err)
        i+=j
    return bin_edges, means, stddev

def is_outlier(mag, rm, bin_edges, avgs, stdevs):
    for i in range(0, len(bin_edges) - 1):
        if mag > bin_edges[i] and mag < bin_edges[i+1]:
            if abs(rm - avgs[i])/stdevs[i] > OUTLIER_SIGMA: return True
            else: return False

def update_files():
    print "Updating files"
    print len([ID for ID in rms])
    f = open(rmsfile,'w')
    for ID in rms:
        f.write("%s %e\n"%(ID,rms[ID]))
    f.close()


    for pkno in range(num_peaks):
        peaks_fname = "stats_peaks_%d.dat"%(pkno)
        f = open(peaks_fname,'w')
        for hid in lsp_peaks:    
            f.write("%s %e %e\n"%(hid, lsp_peaks[hid][pkno][0], lsp_peaks[hid][pkno][1] ))
        f.close()

def not_ok(ID):
    rms[ID] = -1.0
    for i in range(num_peaks):
        if ID not in lsp_peaks:
            lsp_peaks[ID] = {}
        lsp_peaks[ID][i] = ( -1 , -1)
for ID in hatids:
    # If you already have the RMS value, skip
    if not sdata is None:
        if not pk_data is None:
            if ID in sdata['id'] and ID in lsp_peaks: continue

    print "Analyzing %s"%(ID)
    
    # Try to read lightcurve
    fname = "%s/%s-hatlc.csv.gz"%(data_dir,ID)
    try:
        lc = rhlc.read_hatlc(fname)
    except:
        print "Can't read ",ID
        not_ok(ID)
        continue

    if lc is None: 
        not_ok(ID)
        continue

    if len(lc['BJD']) < 2: 
        not_ok(ID)
        continue
    t, x, filts = np.array([ float(v) for v in lc['BJD'] ]),np.array([ float(v) for v in lc[MAG_COL] ]),np.array([ int(v) for v in lc['FLT'] ] )
    ok_inds = []
    for i in range(len(x)):
        if not isnan(x[i]): ok_inds.append(i)
    
    if len(ok_inds) < 2: 
        print "No data in ",ID
        not_ok(ID)
        continue
    
    else:
        t,x,filts = t[ok_inds], x[ok_inds], filts[ok_inds]
        t, x = filter_normalization(t, x, filts)
        mu = np.mean(x)
        rms[ID] = sqrt(np.mean(np.power(x - mu,2)))
        freqs, powers = lsp_raw(t, x)
        periods = np.power(freqs, -1)

        pk_pers, pk_powers = find_n_peaks(periods, powers, num_peaks)
        for i in range(num_peaks):
            if ID not in lsp_peaks:
                lsp_peaks[ID] = {}
            lsp_peaks[ID][i] = (pk_pers[i], pk_powers[i])
        update_files()




mags_cleaned = []
rms_cleaned = []
hatids_cleaned = []
first_power = []
second_power = []
third_power = []
first_period = []
#print rms
for ID in hatids:
    if rms[ID] > 0:
        hatids_cleaned.append(ID)
        mags_cleaned.append(magnitudes[ID])
        rms_cleaned.append(rms[ID])
        first_power.append(lsp_peaks[ID][0][1])
        second_power.append(lsp_peaks[ID][1][1])
        third_power.append(lsp_peaks[ID][2][1])
        first_period.append(lsp_peaks[ID][0][0])

bin_edges, means, stdev = bin(mags_cleaned, rms_cleaned)

f, (ax1, ax2, ax3, ax4) = plt.subplots(4,1)
alphaval = 0.3
ax1.plot(rms_cleaned, first_power, 'o', alpha=alphaval)
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_xlabel("RMS")
ax1.set_ylabel("$P_1$")

ax2.plot(mags_cleaned, first_power,'o', alpha=alphaval, label="$P_1$")
ax2.plot(mags_cleaned, second_power,'o', alpha=alphaval, label="$P_2$")
ax2.plot(mags_cleaned, third_power, 'o', alpha=alphaval, label="$P_3$")
ax2.set_yscale("log")
ax2.set_xlabel("mag")
ax2.set_ylabel("$P$")
ax2.legend(loc='best')

ax3.plot(first_period, first_power,'o', alpha=alphaval)
ax3.set_yscale("log")
ax3.set_xscale("log")
ax3.set_xlabel("Period [days]")
ax3.set_ylabel("$P$")

ax4.plot(mags_cleaned, rms_cleaned, 'o', alpha=alphaval)
ax4.set_yscale("log")
ax4.set_xscale("log")
ax4.set_xlabel("RMS")
ax4.set_ylabel("mag")

plt.show()

candidates = []
for ID in hatids_cleaned:
    if lsp_peaks[ID][0][1] > power_cut or is_outlier(magnitudes[ID], rms[ID], bin_edges, means, stdev): candidates.append(ID)

print "%d candidates"%(len(candidates))

f = open("candidates.dat",'w')
for ID in candidates:
    f.write("%s\n"%(ID))

f.close()

