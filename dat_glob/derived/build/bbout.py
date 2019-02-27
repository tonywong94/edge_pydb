#!/usr/bin/env python

# Output global fit parameters from Bbarolo

import glob
import numpy as np
from astropy.table import Table, Column
from astropy.wcs import WCS

flist = glob.glob('genout/*/ringlog2.txt')

tab=Table(names=['bbRactr','bbDectr','bbVsys','bbInc','bbKinPA'])
col=Column(['UGC05498NED01'],name='Name',dtype='str',description='Galaxy Name')
tab.add_row()
tab.add_column(col,index=0)
tab['bbRactr'].unit = 'deg'
tab['bbRactr'].description = 'R.A. of center determined by Bbarolo'
tab['bbDectr'].unit = 'deg'
tab['bbDectr'].description = 'Dec. of center determined by Bbarolo'
tab['bbVsys'].unit = 'km / s'
tab['bbVsys'].description = 'Systemic velocity determined by Bbarolo (radio-LSR)'
tab['bbInc'].unit = 'deg'
tab['bbInc'].description = 'Inclination determined by Bbarolo'
tab['bbKinPA'].unit = 'arcsec'
tab['bbKinPA'].description = 'Position angle E from N determined by Bbarolo'

i = 0
for file in flist:
    gal = file.split('/')[1]
    if i>0:
        tab.add_row()
    tab['Name'][i] = gal
    dat=np.genfromtxt(file,usecols=[4,5,9,10,11],
        names=['inc','pa','xpos','ypos','vsys'],skip_header=1,max_rows=1)
    tab['bbVsys'][i] = dat['vsys']
    tab['bbInc'][i] = dat['inc']
    tab['bbKinPA'][i] = dat['pa']
    w = WCS('cubes/'+gal+'.co.cmnormsub.fits')
    lon,lat,v,stok = w.wcs_pix2world([[dat['xpos'],dat['ypos'],0,0]],0)[0]
    tab['bbRactr'][i] = lon
    tab['bbDectr'][i] = lat
    i += 1

for cname in tab.colnames:
    if cname == 'bbRactr' or cname == 'bbDectr':
        tab[cname].format='.4f'
    elif cname != 'Name':
        tab[cname].format='.2f'

tab['Name'][np.where(tab['Name']=='NGC4211N')] = 'NGC4211NED02'
tab['Name'][np.where(tab['Name']=='UGC05498')] = 'UGC05498NED01'

tab.write('edge_bbpars.csv', format='ascii.ecsv', delimiter=',')

