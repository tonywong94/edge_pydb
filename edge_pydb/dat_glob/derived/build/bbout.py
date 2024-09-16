#!/usr/bin/env python

# Capture global fit parameters from Bbarolo ringlog2.txt files

import os
import glob
import numpy as np
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.io.fits import getheader
from datetime import datetime

dotypes = ['natv', 'smo7']

for type in dotypes:
    if type == 'natv':
        cubedir = '/Volumes/Scratch2/tonywong/EDGE/comb_de_10/native/'
        fitdir  = '/Volumes/Scratch2/tonywong/EDGE/bbarolo/natv_fitvd_dilmsk/output/'
        fitdir2 = '/Volumes/Data/tonywong/sharenb/bbarolo/smo7_fitvd_dilmsk/output/'
    elif type == 'smo7':
        cubedir = '/Volumes/Scratch2/tonywong/EDGE/comb_de_10/smo7/'
        fitdir  = '/Volumes/Scratch2/tonywong/EDGE/bbarolo/smo7_fitvd_dilmsk/output/'
        fitdir2 = '/Volumes/Data/tonywong/sharenb/bbarolo_7as/gal_nrad8_smolist/output/'

    # -- Original cubes are needed to provide astrometry
    cubelist = sorted(glob.glob(cubedir+'*.co.*norm.fits.gz'))
    print('Found {} cubes in {}'.format(len(cubelist),cubedir))
    gallist = [os.path.basename(file).split('.')[0] for file in cubelist]

    tab=Table(names=['Name','bbRactr','bbDectr','bbVsys','bbKinInc','bbKinPA','bbMask'],
              dtype=('S13','f8','f8','f8','f8','f8','S13'))
    tab['Name'].description = 'Galaxy name'
    tab['bbRactr'].unit = 'deg'
    tab['bbRactr'].description = 'R.A. J2000 of center used by Bbarolo'
    tab['bbDectr'].unit = 'deg'
    tab['bbDectr'].description = 'Dec. J2000 of center used by Bbarolo'
    tab['bbVsys'].unit = 'km / s'
    tab['bbVsys'].description = 'Systemic velocity determined by Bbarolo (radio-LSR)'
    tab['bbKinInc'].unit = 'deg'
    tab['bbKinInc'].description = 'Inclination determined by Bbarolo'
    tab['bbKinPA'].unit = 'deg'
    tab['bbKinPA'].description = 'Position angle E from N determined by Bbarolo'
    tab['bbMask'].description = 'Identifier for cube mask'

    for i, gal in enumerate(gallist):
        hdr = getheader(cubelist[i])
        hdr['specsys'] = 'LSRK'
        hdr['ctype3']  = 'VRAD'
        hdr['velref']  = 257
        w = WCS(hdr)
        ringlog = fitdir + gal[:8] + '/rings_final2.txt'
        try:
            dat=Table.read(ringlog, format='ascii', include_names=
                    ['INC(deg)','P.A.(deg)','XPOS(pix)','YPOS(pix)','VSYS(km/s)'])
            mask = 'de10_'+type+'_dil'
        except:
            print('Problem reading',ringlog)
            ringlog = fitdir2 + gal[:8] + '/ringlog2.txt'
            if os.path.isfile(ringlog):
                dat=Table.read(ringlog, format='ascii', include_names=
                    ['INC(deg)','P.A.(deg)','XPOS(pix)','YPOS(pix)','VSYS(km/s)'])
                mask = 'smooth'
            else:
                print('Problem reading',ringlog)
                mask = 'fail'
        if mask != 'fail':
            print(gal,'succeeded')
            lon, lat, v, stok = w.wcs_pix2world(
                    [[dat['XPOS(pix)'][0],dat['YPOS(pix)'][0],0,0]],0)[0]
            vsys = dat['VSYS(km/s)'][0]
            inc  = dat['INC(deg)'][0]
            pa   = dat['P.A.(deg)'][0]
            tab.add_row([gal, lon, lat, vsys, inc, pa, mask])

    for cname in tab.colnames:
        if cname == 'bbRactr' or cname == 'bbDectr':
            tab[cname].format='.4f'
        elif cname == 'bbVsys' or cname == 'bbKinInc' or cname == 'bbKinPA':
            tab[cname].format='.2f'

    tab.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    tab.meta['comments'] = ('Global fit parameters from Bbarolo on '
            +type+' CO data')
    tab.write('edge_bbpars_'+type+'.csv', format='ascii.ecsv', delimiter=',',
              overwrite=True)

