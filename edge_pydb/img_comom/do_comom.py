#!/usr/bin/env python

# Combine the 7" CO moment maps into binary tables.

import sys
sys.path.append('../edge_pydb')
from fitsextract import fitsextract
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack

msktyp = ['dil', 'smo']

# Get the orientation parameters from LEDA
globaldir = '../dat_glob/'
ort = Table.read(globaldir+'external/edge_leda.csv', format='ascii.ecsv')
ort.add_index('Name')

for i, msk in enumerate(msktyp):
    filelist = glob.glob('fitsdata/*.co.smo7_'+msk+'.emom0max.fits.gz')
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
                                keepnan=True, stride=[3,3,1],
                                ra_gc=15*ort.loc[gal]['ledaRA'],
						        dec_gc=ort.loc[gal]['ledaDE'],
                                pa=ort.loc[gal]['ledaPA'],
                                inc=ort.loc[gal]['ledaIncl'],
                                ortlabel='LEDA')
                jointb = join(tab0, addtb, keys=['ix','iy'])
                tab0 = jointb
            else:
                newcol = Column(data=[np.nan]*len(tab0), name=type, 
                                unit=unit[j], dtype='f4')
                tab0.add_column(newcol)
        tablelist.append(tab0)
    if len(tablelist) > 0:
        t_merge = vstack(tablelist)
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
