#!/usr/bin/env python

# Combine the CALIFA data into binary tables.

import glob
import os
from astropy.table import Table, Column, join, vstack
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import sys
sys.path.append('../edge_pydb')
from fitsextract import fitsextract, getlabels

# Get the orientation parameters from LEDA
globaldir = '../dat_glob/'
ort = Table.read(globaldir+'external/edge_leda.csv', format='ascii.ecsv')
ort.add_index('Name')

# Read the FITS data
codir = '../img_comom/fitsdata/'
cadir = 'fitsdata/'
prodtype = ['ELINES', 'SFH', 'SSP', 'indices', 'flux_elines']

for prod in prodtype:
    zsel, labels, units, nsel = getlabels(prod)
    filelist = [fn for fn in glob.glob(cadir+'*'+prod+'*.fits.gz')
            if not os.path.basename(fn).startswith('x')]
    print('\n',filelist)
    rglist = []
    smlist = []

    for file in filelist:
        base = os.path.basename(file)
        if prod not in ['indices', 'flux_elines']:
            gal = base.split('.')[0]
        elif prod == 'indices':
            gal = base.split('.')[2]
        else:
            gal = base.split('.')[1]

        print('\nWorking on galaxy {} product {} nsel={}'.format(
            gal, prod, nsel))

        # Generate output header using CO astrometry
        cohd = fits.getheader(codir+gal+'.co.smo7_dil.emom0max.fits.gz')
        cahd = fits.getheader(cadir+'x'+base, ignore_missing_end=True)
        outhd = cahd.copy()
        for key in ['NAXIS1', 'NAXIS2', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 
                    'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2']:
            outhd[key] = cohd[key]
        cdkeys = ['CD1_1','CD1_2','CD2_1','CD2_2','CD1_3','CD2_3','CD3_1','CD3_2','CD3_3']
        for key in cdkeys:
            if key in outhd.keys():
                del outhd[key]

        # First process the native resolution file since it has the astrometry
        hdu = fits.open(cadir+'x'+base, ignore_missing_end=True)[0]
        hdu.data[hdu.data==0] = np.nan
        newim,foot = reproject_interp(hdu, outhd, independent_celestial_slices=True)
        #fits.writeto(base.replace('fits','rgd.fits'), newim, outhd, overwrite=True)
        tab0 = fitsextract(newim, header=outhd, keepnan=True, stride=[3,3,1], 
            bunit=units, col_lbl=labels, zselect=zsel, ra_gc=15*ort.loc[gal]['ledaRA'],
			dec_gc=ort.loc[gal]['ledaDE'], pa=ort.loc[gal]['ledaPA'],
            inc=ort.loc[gal]['ledaIncl'], ortlabel='LEDA', first=True)
        gname = Column([np.string_(gal)]*len(tab0), name='Name', 
                       description='Galaxy Name')
        tab0.add_column(gname, index=0)
        rglist.append(tab0)

        # Then process the smoothed file
        hdu = fits.open(cadir+base, ignore_missing_end=True)[0]
        hdu.data[hdu.data==0] = np.nan
        hdu.header = cahd
        newim,foot = reproject_interp(hdu, outhd, independent_celestial_slices=True)
        #fits.writeto(base.replace('fits','sm.fits'), newim, outhd, overwrite=True)
        tab1 = fitsextract(newim, header=outhd, keepnan=True, stride=[3,3,1], 
            bunit=units, col_lbl=labels, zselect=zsel, ra_gc=15*ort.loc[gal]['ledaRA'],
			dec_gc=ort.loc[gal]['ledaDE'], pa=ort.loc[gal]['ledaPA'],
            inc=ort.loc[gal]['ledaIncl'], ortlabel='LEDA', first=True)
        gname = Column([np.string_(gal)]*len(tab1), name='Name', 
                       description='Galaxy Name')
        tab1.add_column(gname, index=0)
        smlist.append(tab1)

    if len(rglist) > 0:
        rg_merge = vstack(rglist)
    print(rg_merge.colnames)
    print('There are',len(rg_merge),'rows in native table')

    if len(smlist) > 0:
        sm_merge = vstack(smlist)
    print(sm_merge.colnames)
    print('There are',len(sm_merge),'rows in smoothed table')

    if (len(filelist) > 1):
        outname = 'edge'
    else:
        outname = gal
    if prod == prodtype[0]:
        rg_merge.write(outname+'.pipe3d.hdf5', path=prod+'_rg', overwrite=True, 
                serialize_meta=True, compression=True)
    else:
        rg_merge.write(outname+'.pipe3d.hdf5', path=prod+'_rg', append=True, 
                serialize_meta=True, compression=True)
    sm_merge.write(outname+'.pipe3d.hdf5', path=prod+'_sm', append=True, 
                serialize_meta=True, compression=True)

