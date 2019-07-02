#!/usr/bin/env python

from edgeprof import edgeprof
import numpy as np
from astropy import units as u
from astropy.table import Table
import os

basedir = '/Volumes/Scratch2/tonywong/EDGE/'
dbdir   = '../../../'

# Make profiles from smoothed-mask mom0 at native and 7" resolution
dotypes = ['de20_smo', 'smo7_smo']
datdir  = [basedir+'comb_de/native/mom/', basedir+'comb_de/smo7/mom/']
radbin  = [2, 3]

rfpars = Table.read(dbdir+'dat_glob/derived/edge_rfpars.csv', format='ascii.ecsv')
glist=rfpars['Name'].tolist()
glist.sort()

for i, type in enumerate(dotypes):
    for gal in glist:
        idx = np.where(rfpars['Name']==gal)[0][0]
        if (rfpars['rfInc'][idx] < 85. and ~np.isnan(rfpars['rfPA'][idx])):
            print('')
            edgeprof(gname=gal, type=type, delr=radbin[i], datdir=datdir[i], 
                     dbdir=dbdir, blanksig=0, fmask=0.1)
            edgeprof(gname=gal, type=type, delr=radbin[i], datdir=datdir[i], 
                     dbdir=dbdir, blanksig=1, replace=True, fmask=0.1)
        else:
            print('\nGalaxy {} is excluded: INC={} PA={}'.format(
                gal, rfpars['rfInc'][idx], rfpars['rfPA'][idx]))

