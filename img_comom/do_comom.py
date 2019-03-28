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

for msk in msktyp:
    filelist = glob.glob('fitsdata/*.co.smo7_'+msk+'.mom0.fits.gz')
    tablelist=[]
    for file in filelist:
        gal = os.path.basename(file).split('.')[0]
        tab0 = fitsextract(file, keepnan=True, stride=[3,3,1])
        tab0['imgdata'].name = 'mom0'
        gname = Column([np.string_(gal)]*len(tab0), name='Name', description='Galaxy Name')
        tab0.add_column(gname, index=0)
        print(tab0[20:50])
        if msk == 'smo':
            dotypes = ['emom0', 'emom0max']
        else:
            dotypes = ['emom0', 'emom0max', 'mom1', 'emom1', 'mom2', 'emom2', 'snrpk']
        for type in dotypes:
            addtb = fitsextract('fitsdata/'+gal+'.co.smo7_'+msk+'.'+type+'.fits.gz', 
                                keepnan=True, stride=[3,3,1])
            addtb['imgdata'].name = type
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
    t_merge.write(outname+'.comom.smo7_'+msk+'.hdf5', path='data', overwrite=True, 
                serialize_meta=True, compression=True)
