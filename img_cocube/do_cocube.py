#!/usr/bin/env python

# Combine the 7" CO cubes into binary tables.

import sys
sys.path.append('../edge_pydb')
from fitsextract import fitsextract
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack

filelist = glob.glob('fitsdata/*.co.smo7msk.K.fits.gz')
tablelist=[]
for file in filelist:
    # Read the FOV masked, gain-corrected cube
    gal = os.path.basename(file).split('.')[0]
    tab0 = fitsextract(file, bunit='K', col_lbl='co_data',
                        keepnan=True, stride=[3,3,1])
    gname = Column([np.string_(gal)]*len(tab0), name='Name', description='Galaxy Name')
    tab0.add_column(gname, index=0)
    print(tab0[20:50])
    # Read the other images: noise, dilated mask cube, smooth mask cube
    files  = ['smo7_str.ecube', 'smo7_dil.mask', 'smo7_smo.mask']
    labels = ['co_rms', 'co_dilmsk', 'co_smomsk']
    unit   = ['K', '', '']
    for j, file in enumerate(files):
        getfile = 'fitsdata/'+gal+'.co.'+file+'.fits.gz'
        if os.path.exists(getfile):
            addtb = fitsextract(getfile, bunit=unit[j], col_lbl=labels[j], 
                            keepnan=True, stride=[3,3,1])
            jointb = join(tab0, addtb)
            tab0 = jointb
        else:
            newcol = Column(data=[np.nan]*len(tab0), name=labels[j], 
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
t_merge.write(outname+'.cocube_smo7.hdf5', path='data', overwrite=True, 
        serialize_meta=True, compression=True)
