#!/usr/bin/env python

import glob
import os
from astropy.table import Table,Column,vstack,join
from datetime import datetime

files=glob.glob('EDGE_CO_vrot/*vrot.txt')

tablelist=[]
for file in files:
    txttab = Table.read(file, format='ascii', names=('radius','Vrot','e_Vrot'))
    gal = os.path.basename(file).split('_')[0]
    # --- Add the Galaxy Name
    nrows = len(txttab)
    gname = Column([gal]*nrows, name='Name', description='Galaxy Name')
    txttab.add_column(gname, index=0)
    # --- Add the beam smearing corrected data
    file2 = file.replace('CO_vrot','CO_vrot_bsc')
    bsctab = Table.read(file2, format='ascii', names=('radius','Vrot_bsc','e_Vrot_bsc'))
    galtab = join(txttab, bsctab, keys='radius')
    tablelist.append(galtab)

t_jamprof = vstack(tablelist)
t_jamprof['radius'].unit = 'arcsec'
t_jamprof['Vrot'].description = 'Rotation curve'
t_jamprof['Vrot'].unit = 'km / s'
t_jamprof['e_Vrot'].description = 'Error in rotation curve'
t_jamprof['e_Vrot'].unit = 'km / s'
t_jamprof['Vrot_bsc'].description = 'Beam smearing corrected rotation curve'
t_jamprof['Vrot_bsc'].unit = 'km / s'
t_jamprof['e_Vrot_bsc'].description = 'Error in corrected rotation curve'
t_jamprof['e_Vrot_bsc'].unit = 'km / s'
t_jamprof.meta['comments'] = ('CO Rotation Curves from 2018MNRAS.477..254L')
t_jamprof.meta['date'] = datetime.today().strftime('%Y-%m-%d')
t_jamprof.write('jam_rotcurves.csv', overwrite=True, delimiter=',', format='ascii.ecsv')
