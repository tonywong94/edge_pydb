#!/usr/bin/env python

# Comparison of moment map methods in edge_pydb.

import numpy as np
from edge_pydb import EdgeTable
import matplotlib.pyplot as plt
from edge_pydb.plotting import gridplot

momcols  = ['Name','ix','iy','mom0_12','mom0_13']

strmsk = EdgeTable('edge_carma_allpix.2d_smo7.hdf5', path='comom_str', cols=momcols)
strmsk['mom0_12'].name = '(a) no masking'
dilmsk = EdgeTable('edge_carma_allpix.2d_smo7.hdf5', path='comom_dil', cols=momcols)
dilmsk['mom0_12'].name = '(b) dilated mask'
smomsk = EdgeTable('edge_carma_allpix.2d_smo7.hdf5', path='comom_smo', cols=momcols)
smomsk['mom0_12'].name = '(c) smoothed mask'
print('Working on full EDGE database')

# consolidate the tables and calculate sSFR
strmsk.join(dilmsk, keys=['Name', 'ix', 'iy'])
strmsk.join(smomsk, keys=['Name', 'ix', 'iy'])
print(strmsk.colnames)

plt.rcParams.update({'font.size': 20})
gridplot(edgetab=strmsk, gallist='NGC4047', 
         columnlist=['(a) no masking','(b) dilated mask','(c) smoothed mask'],
         nx=3, ny=1, plotstyle='image', clipedge=True, pad=16, vshow=True, 
         pdfname='ngc4047_mom.pdf')

