#!/usr/bin/env python

# Merge the radial profiles for all galaxies with a given masking/resolution.
# de20_smo: native resolution, smoothed mask
# smo7_smo: smoothed to 7", smoothed mask

import os
from astropy.table import Table,Column,vstack
import numpy as np
import glob

for msk in ['de20_smo', 'smo7_smo']:
    filelist = glob.glob('rproftxt/*.'+msk+'b0.rprof.txt')
    tablelist=[]
    for file in filelist:
        gal = os.path.basename(file).split('.')[0]
        # --- Try to read the table
        try:
            rproftxt = Table.read(file, format='ascii.ecsv')
        except:
            print('reading '+file+' failed')
            continue
        rproftxt.rename_column('radius', 'r_arcs')
        rproftxt.rename_column('wtmean', 'ico_avg')
        rproftxt.rename_column('rms', 'ico_rms')
        ordlist=['r_arcs', 'r_kpc', 'r_r25', 'npix', 'ngood', 'ico_avg', 
                 'ico_rms', 'detlim', 'sigmol', 'cumlum', 'cummass']
        outtable = rproftxt[ordlist]
        # --- Add the Galaxy Name
        if len((outtable['ngood']>0).nonzero()[0]) > 0:
            nrows = len(outtable)
            gname = Column([gal]*nrows, name='Name', description='Galaxy Name')
            outtable.add_column(gname, index=0)
            rmax = outtable['r_arcs'][outtable['ngood']>0].max()
            tablelist.append(outtable[outtable['r_arcs']<=rmax])
        else:
            continue
    if len(tablelist) > 0:
        t_merge = vstack(tablelist)
        t_merge['r_arcs'].description = 'outer radius of ring in arcsec'
        t_merge['r_kpc'].description = 'outer radius of ring in kpc'
        t_merge['r_r25'].description = 'outer radius of ring as frac of opt radius'
        t_merge['npix'].description = 'number of available pixels in the ring'
        t_merge['ngood'].description = 'number of unmasked pixels in the ring'
        t_merge['ico_avg'].description = 'average face-on intensity in ring if ngood/npix > 0.1 and setting blanks to zero'
        t_merge['ico_rms'].description = 'rms face-on intensity in ring for unmasked pixels'
        t_merge['detlim'].description = '3 sigma detection limit based on emom0max'
        t_merge['sigmol'].description = 'average face-on surface density including He with alphaco=4.3'
        t_merge['cumlum'].description = 'total CO luminosity within the given radius'
        t_merge['cummass'].description = 'total molecular gas mass within the given radius'
        t_merge.write('../rprof_'+msk+'.csv',overwrite=True,delimiter=',',format='ascii.ecsv')

