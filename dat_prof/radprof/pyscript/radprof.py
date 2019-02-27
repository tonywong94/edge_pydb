#!/usr/bin/env python

from edgeprof import *
import numpy as np
#from astropy import units as u
from astropy.table import Table, Column
#from astropy.io import fits 
#from astropy.wcs import WCS
from matplotlib import pyplot as plt
import os


# Plot radial profiles for all galaxies
"""
Currently generating profiles from the smoothed masked mom-0 maps.
Blue line uses zeroes to replace blanks and plots annuli >20% non-blanked.
Green shading uses 1-sigma to replace blanks and <1 sigma values and plots
annuli >10% non-replaced.
"""
type='smo'
redo=True

dbdir = '/Volumes/Scratch2/tonywong/EDGE/edge-sql-base/global_values/external/'
db = Table.read(dbdir+'DETableFinal.csv', format='ascii.csv')
db['ledaRA'].unit = 'hourangle'
db['ledaDE'].unit = 'deg'
db['ledaD25'].unit = 'arcmin'
db['caDistMpc'].unit = 'Mpc'
db['coScaleMol'].unit = 'kpc'
db['coNormMol'].unit = 'solMass/pc2'
# db2 = Table.read(dbdir+'external/edge_califa.csv', format='ascii.ecsv')
# t2 = Table([db2['Name'], db2['caDistMpc']])
# db3 = Table.read(dbdir+'derived/edge_rfpars.csv', format='ascii.ecsv')
# db3.rename_column('rfName', 'Name')
# db = join(db1, t2, keys='Name')
# db = join(db, db3, keys='Name')
db.sort('ledaRA')
glist=db['Name'].tolist()
if redo==True:
    for gal in glist:
        idx = np.where(db['Name']==gal)[0][0]
        if (db['coInc'][idx] < 85. and ~np.isnan(db['coPA'][idx])):
            edgeprof(gname=gal, type=type, delr=2, db=db, blanksig=0, deproj=True,
                fmask=0.2)
            edgeprof(gname=gal, type=type, delr=2, db=db, blanksig=1, deproj=True,
                ndblank=True, fmask=0.1)

# Make the plot
nx=7
ny=5
pages = int(np.ceil(float(len(glist)) / (nx*ny)))
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

for rtype in ['arcsec', 'r25', 'kpc']:
    for itype in ['wtmean', 'sigmol']:
        if rtype == 'r25':
            rlabel = 'Radius / '+r'$R_{25}$'
        else:
            rlabel = 'Radius ['+rtype+']'
        if itype == 'sigmol':
            prefix = 'sig'
            ilabel = r'log $\Sigma_{mol}$ [$M_\odot$ pc$^{-2}$]'
        else:
            prefix = 'i'
            ilabel = r'log $I_{CO}$ [K km/s]'
        for num in range(0,pages):
            aa = nx*ny*num
            bb = nx*ny+aa
            gals = glist[aa:bb]
            nrows = len(gals)//nx
            figure = plt.figure(0)
            figure.set_size_inches(nx*4.5, nrows*4.)
            for i, gal in enumerate(gals):
                igal = num*nx*ny + i
                row,col = divmod(i,nx)
                ax = plt.subplot2grid((nrows,nx),(row,col))
                pltcoprof(gal, ax, rtype=rtype, itype=itype, imtype=type, ebar=False,
                    mod0=db['coNormMol'][igal], modh=db['coScaleMol'][igal])
                if col != 0:
                    ax.set_yticklabels([])
                if row != nrows-1:
                    ax.set_xticklabels([])
            figure.subplots_adjust(hspace=0.1)
            figure.subplots_adjust(wspace=0.1)
            if num == 3:
                figure.text(0.5, 0.03, rlabel, ha='center', fontsize=24)
            else:
                figure.text(0.5, 0.06, rlabel, ha='center', fontsize=24)
            figure.text(0.09, 0.5, ilabel, va='center', 
                rotation='vertical', fontsize=24)
            plt.savefig(prefix+'prof_'+type+'_'+rtype+'-'+str(num)+'.pdf', 
                bbox_inches='tight')
            plt.close()
        os.system('pdfunite '+prefix+'prof_'+type+'_'+rtype+"-*.pdf "
            +prefix+'prof_'+type+'_'+rtype+'.pdf')
        os.system("rm -f "+prefix+'prof_'+type+'_'+rtype+"-*.pdf")
