#!/usr/bin/env python

# Combine the 7" CO moment maps into binary tables.

import sys
sys.path.append('../edge_pydb')
from fitsextract import fitsextract
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack
#from astropy.wcs import WCS

msktyp = ['dil', 'smo']

for i, msk in enumerate(msktyp):
    filelist = glob.glob('fitsdata/*.co.smo7_'+msk+'.mom0.fits.gz')
    tablelist=[]
    for file in filelist:
    	# Read the mom0 image
        gal = os.path.basename(file).split('.')[0]
        tab0 = fitsextract(file, bunit='K km/s', col_name='mom0',
        					keepnan=True, stride=[3,3,1])
        gname = Column([np.string_(gal)]*len(tab0), name='Name', description='Galaxy Name')
        tab0.add_column(gname, index=0)
        print(tab0[20:50])
    	# Read the other images
        if msk == 'smo':
            dotypes = ['emom0', 'emom0max']
            unit    = ['K km/s', 'K km/s']
        else:
            dotypes = ['emom0', 'emom0max', 'mom1', 'emom1', 'mom2', 'emom2', 'snrpk']
            unit    = ['K km/s', 'K km/s', 'km/s', 'km/s', 'km/s', 'km/s', '']
        for j, type in enumerate(dotypes):
            addtb = fitsextract('fitsdata/'+gal+'.co.smo7_'+msk+'.'+type+'.fits.gz', 
                                bunit=unit[j], col_name=type, keepnan=True, 
                                stride=[3,3,1])
            #addtb['imgdata'].name = type
            jointb = join(tab0, addtb)
            tab0 = jointb
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
