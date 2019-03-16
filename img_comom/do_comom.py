#!/usr/bin/env python

# Combine the comom files into a binary table.

import sys
sys.path.append('../edge_pydb')
from fitsextract import *
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack
#from astropy.wcs import WCS

filelist = glob.glob('fitsdata/*.co.smo7_dil.mom0.fits.gz')
tablelist=[]
for file in filelist:
    gal = os.path.basename(file).split('.')[0]
    tab0 = fitsextract(file, keepnan=True, stride=[3,3,1])
    tab0['imgdata'].name = 'mom0'
    gname = Column([np.string_(gal)]*len(tab0), name='Name', description='Galaxy Name')
    tab0.add_column(gname, index=0)
    print(tab0[20:50])
    for type in (['emom0', 'emom0max', 'mom1', 'emom1', 'mom2', 'emom2', 'snrpk']):
        addtb = fitsextract('fitsdata/'+gal+'.co.smo7_dil.'+type+'.fits.gz', 
                            keepnan=True, stride=[3,3,1])
        addtb['imgdata'].name = type
        jointb = join(tab0, addtb)
        tab0 = jointb
    tablelist.append(tab0)
if len(tablelist) > 0:
    t_merge = vstack(tablelist)
print(t_merge[20:50])
t_merge.write('edge.comom.smo7_dil.hdf5', path='data', overwrite=True, 
            serialize_meta=True, compression=True)
