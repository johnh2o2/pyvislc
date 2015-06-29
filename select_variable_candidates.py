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

MAG_COL = "TF1"
MAGNITUDE = 'vmag'
NLCS = 1000
OUTLIER_SIGMA = 3
NBINS = 10


rms = {}
magnitudes = {}

print "Loading list of HAT sources..."
sources = pd.read_csv("/Users/jah5/Documents/Fall2014_Gaspar/data/hat-field-219-object-list.csv")
sources = sources.iloc[:-1]
#sources = sources.sort('vmag')


print "Now selecting %d random HAT sources"%(NLCS)
indices = np.arange(0,len(sources))
#np.random.shuffle(indices)
selected_indices = indices[0:NLCS]

hatids = [ hid for hid in sources.iloc[selected_indices]['#hat_id'] ]
for i in selected_indices:
    s = sources.iloc[i]
    hid = s['#hat_id']
    magnitudes[hid] = s[MAGNITUDE]

# FETCH lightcurves from the server
#print "Fetching lightcurves ... "
#fetch_lcs(hatids)

# Calculate stats (rms, frequency fits, etc)
print "Calculating stats..."
dt_stats = np.dtype([('id','S15'),('rms',np.float_)])
rmsfile = "stats_rms.dat"
if os.path.exists(rmsfile):
    print "Found stats_rms.dat, reading in data"
    f = open(rmsfile,'r')
    sdata = np.loadtxt(rmsfile,dtype=dt_stats)
    for i,hid in enumerate(sdata['id']):
        if hid not in hatids: continue
        rms[hid] = sdata[i]['rms']
    f.close()
else: sdata = None

def bin(mag, rm, nbins=NBINS, remove_outliers=True, outlier_dist=OUTLIER_SIGMA):
    sources = np.zeros(len(rm), dtype=np.dtype([('mag', np.float_), ('rms', np.float_)]))
    for i in range(0,len(mag)):
        sources[i]['mag'] = mag[i]
        sources[i]['rms'] = rm[i]
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
            print "bin %d, nsrc = %d"%( i, len(bin_sources))
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

for ID in hatids:
    # If you already have the RMS value, skip
    if not sdata is None:
        if ID in sdata['id']: continue

    print "Analyzing %s"%(ID)
    
    # Try to read lightcurve
    fname = "%s/%s-hatlc.csv.gz"%(data_dir,ID)
    try:
        lc = rhlc.read_hatlc(fname)
    except:
        print "Can't read ",ID
        continue

    mags = np.array([ m for m in lc[MAG_COL] if not isnan(m) ])

    if len(mags) < 2: 
        print "No data in ",ID
        rms[ID] = -1.0
        continue
    
    else:
        mu = np.mean(mags)
        rms[ID] = sqrt(np.mean(np.power(mags - mu,2)))

f = open(rmsfile,'w')
for ID in hatids:
    f.write("%s %e\n"%(ID,rms[ID]))
f.close()

mags_cleaned = []
rms_cleaned = []
hatids_cleaned = []
for ID in hatids:
    if rms[ID] > 0:
        hatids_cleaned.append(ID)
        mags_cleaned.append(magnitudes[ID])
        rms_cleaned.append(rms[ID])

bin_edges, means, stdev = bin(mags_cleaned, rms_cleaned)

candidates = []
for ID in hatids_cleaned:
    if is_outlier(magnitudes[ID], rms[ID], bin_edges, means, stdev): candidates.append(ID)

print "%d candidates"%(len(candidates))

f = open("candidates.dat",'w')
for ID in candidates:
    f.write("%s\n"%(ID))

f.close()

