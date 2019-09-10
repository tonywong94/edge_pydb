#!/usr/bin/env python

# Combine the 7" CO moment maps into binary tables.

from datetime import datetime
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack
from edge_pydb import EdgeTable
from edge_pydb.conversion import msd_co
from edge_pydb.fitsextract import fitsextract

msktyp = ['dil', 'smo']

# Get the orientation parameters from LEDA
ort = EdgeTable('edge_leda.csv', cols=['Name', 'ledaRA', 'ledaDE', 'ledaPA', 'ledaIncl'])
#ort = EdgeTable('edge_rfpars.csv', cols=['Name', 'rfPA', 'rfInc', 'rfKinRA', 'rfKinDecl'])
ort.add_index('Name')

for i, msk in enumerate(msktyp):
    filelist = sorted(glob.glob('fitsdata/*.co.smo7_'+msk+'.emom0max.fits.gz'))
    tablelist=[]

    for file in filelist:
        print('Reading',file)
    	# Read the emom0max image first (available for all galaxies)
        gal = os.path.basename(file).split('.')[0]
        tab0 = fitsextract(file, bunit='K km/s', col_lbl='emom0max',
        					keepnan=True, stride=[3,3,1],
        					ra_gc=15*ort.loc[gal]['ledaRA'],
						    dec_gc=ort.loc[gal]['ledaDE'],
                            pa=ort.loc[gal]['ledaPA'],
                            inc=ort.loc[gal]['ledaIncl'],
                            ortlabel='LEDA', first=True)
        gname = Column([np.string_(gal)]*len(tab0), name='Name', description='Galaxy Name')
        tab0.add_column(gname, index=0)
        print(tab0[20:50])
    	# Read the other images
        if msk == 'smo':
            dotypes = ['mom0', 'emom0']
            unit    = ['K km/s', 'K km/s']
        else:
            dotypes = ['mom0', 'emom0', 'mom1', 'emom1', 'mom2', 'emom2', 'snrpk']
            unit    = ['K km/s', 'K km/s', 'km/s', 'km/s', 'km/s', 'km/s', '']
        for j, type in enumerate(dotypes):
            getfile = 'fitsdata/'+gal+'.co.smo7_'+msk+'.'+type+'.fits.gz'
            if os.path.exists(getfile):
                print('Reading',getfile)
                addtb = fitsextract(getfile, bunit=unit[j], col_lbl=type, 
                                keepnan=True, stride=[3,3,1])
                jointb = join(tab0, addtb, keys=['ix','iy'])
                tab0 = jointb
            else:
                newcol = Column(data=[np.nan]*len(tab0), name=type, 
                                unit=unit[j], dtype='f4')
                tab0.add_column(newcol)
        # Add the H2 column density, not deprojected
        sigmol = msd_co(tab0['mom0'], name='sigmol')
        e_sigmol = msd_co(tab0['emom0'], name='e_sigmol')
        tab0.add_columns([sigmol, e_sigmol])
        tablelist.append(tab0)

    if len(tablelist) > 0:
        t_merge = vstack(tablelist)
        t_merge['emom0max'].description = 'error in mom0 assuming 200 km/s window'
        t_merge['mom0'].description = 'integrated intensity using {} mask'.format(msk)
        t_merge['emom0'].description = 'error in mom0 assuming {} mask'.format(msk)
        t_merge['sigmol'].description = 'apparent H2+He surf density not deprojected'
        t_merge['e_sigmol'].description = 'error in sigmol not deprojected'
        if msk == 'dil':
            t_merge['mom1'].description = 'intensity wgtd mean velocity using {} mask'.format(msk)
            t_merge['emom1'].description = 'error in mom1 assuming {} mask'.format(msk)
            t_merge['mom2'].description = 'intensity wgtd vel disp using {} mask'.format(msk)
            t_merge['emom2'].description = 'error in mom2 assuming {} mask'.format(msk)
            t_merge['snrpk'].description = 'peak signal to noise ratio'
        t_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
        print(t_merge[20:50])
    if (len(filelist) > 1):
        outname = 'edge'
    else:
        outname = gal
    if i == 0:
    	t_merge.write(outname+'.comom_smo7.hdf5', path=msk, overwrite=True, 
                serialize_meta=True, compression=True)
    else:
    	t_merge.write(outname+'.comom_smo7.hdf5', path=msk, append=True, 
                serialize_meta=True, compression=True)
