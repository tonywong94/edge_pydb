#!/usr/bin/env python

# Combine the CALIFA data into binary tables.

from datetime import datetime
import glob
import os
from astropy.table import Table, Column, join, vstack
from astropy import units as u
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from edge_pydb import EdgeTable
from edge_pydb.conversion import stmass_pc2, sfr_ha, ZOH_M13, bpt_type
from edge_pydb.fitsextract import fitsextract, getlabels

# Get the orientation parameters from LEDA
ort = EdgeTable('edge_leda.csv', cols=['Name', 'ledaRA', 'ledaDE', 'ledaPA', 'ledaAxIncl'])
ort.add_index('Name')

# Get the distance from the CALIFA table
dist = EdgeTable('edge_califa.csv', cols=['Name', 'caDistP3d'])
dist.add_index('Name')

# Read the FITS data
codir = '../img_comom/fitsdata/'
cadir = 'fitsdata/'
prodtype = ['ELINES', 'SFH', 'SSP', 'indices', 'flux_elines']

for prod in prodtype:
    zsel, labels, units, nsel = getlabels(prod)
    filelist = [fn for fn in sorted(glob.glob(cadir+'*'+prod+'*.fits.gz'))
            if not os.path.basename(fn).startswith('x')]
    print('\n',filelist)
    rglist = []
    smlist = []

    for file in filelist:
        base = os.path.basename(file)
        if prod not in ['indices', 'flux_elines']:
            gal = base.split('.')[0]
        elif prod == 'indices':
            gal = base.split('.')[2]
        else:
            gal = base.split('.')[1]

        print('\nWorking on galaxy {} product {} nsel={}'.format(
            gal, prod, nsel))

        # Generate output header using CO astrometry
        cohd = fits.getheader(codir+gal+'.co.smo7_dil.snrpk.fits.gz')
        cahd = fits.getheader(cadir+'x'+base, ignore_missing_end=True)
        outhd = cahd.copy()
        for key in ['NAXIS1', 'NAXIS2', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 
                    'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2']:
            outhd[key] = cohd[key]
        cdkeys = ['CD1_1','CD1_2','CD2_1','CD2_2','CD1_3','CD2_3','CD3_1','CD3_2','CD3_3']
        for key in cdkeys:
            if key in outhd.keys():
                del outhd[key]
            if key in cahd.keys():
                del cahd[key]

        # First process the native resolution file since it has the astrometry
        hdu = fits.open(cadir+'x'+base, ignore_missing_end=True)[0]
        hdu.data[hdu.data==0] = np.nan
        hdu.header = cahd
        newim = reproject_interp(hdu, outhd, order=0, return_footprint=False)
        #fits.writeto(base.replace('fits','rg.fits'), newim, outhd, overwrite=True)
        rglabels = [s+'_rg' for s in labels]
        tab0 = fitsextract(newim, header=outhd, keepnan=True, stride=[3,3,1], 
            bunit=units, col_lbl=rglabels, zselect=zsel, ra_gc=15*ort.loc[gal]['ledaRA'],
			dec_gc=ort.loc[gal]['ledaDE'], pa=ort.loc[gal]['ledaPA'],
            inc=ort.loc[gal]['ledaAxIncl'], ortlabel='LEDA', first=True)
        gname = Column([np.string_(gal)]*len(tab0), name='Name', 
                       description='Galaxy Name')
        tab0.add_column(gname, index=0)
        rglist.append(tab0)

        # Then process the smoothed file
        hdu = fits.open(cadir+base, ignore_missing_end=True)[0]
        hdu.data[hdu.data==0] = np.nan
        hdu.header = cahd
        newim = reproject_interp(hdu, outhd, order=0, return_footprint=False)
        #fits.writeto(base.replace('fits','sm.fits'), newim, outhd, overwrite=True)
        smlabels = [s+'_sm' for s in labels]
        tab1 = fitsextract(newim, header=outhd, keepnan=True, stride=[3,3,1], 
            bunit=units, col_lbl=smlabels, zselect=zsel, ra_gc=15*ort.loc[gal]['ledaRA'],
			dec_gc=ort.loc[gal]['ledaDE'], pa=ort.loc[gal]['ledaPA'],
            inc=ort.loc[gal]['ledaAxIncl'], ortlabel='LEDA', first=True)
        gname = Column([np.string_(gal)]*len(tab1), name='Name', 
                       description='Galaxy Name')
        tab1.add_column(gname, index=0)
        smlist.append(tab1)
        
        # Add additional columns
        if prod == 'ELINES':
            sfr0, sfrext0 = sfr_ha(tab0['Halpha_rg'], tab0['Hbeta_rg'], name='sigsfr_rg')
            sfr1, sfrext1 = sfr_ha(tab1['Halpha_sm'], tab1['Hbeta_sm'], name='sigsfr_sm')
            tab0.add_column(sfr0)
            tab0.add_column(sfrext0)
            tab1.add_column(sfr1)
            tab1.add_column(sfrext1)
        elif prod == 'flux_elines':
            sfr0, sfrext0 = sfr_ha(tab0['flux_Halpha_rg'], tab0['flux_Hbeta_rg'], 
                                   name='flux_sigsfr_rg')
            sfr1, sfrext1 = sfr_ha(tab1['flux_Halpha_sm'], tab1['flux_Hbeta_sm'], 
                                   name='flux_sigsfr_sm')
            tab0.add_column(sfr0)
            tab0.add_column(sfrext0)
            tab1.add_column(sfr1)
            tab1.add_column(sfrext1)
            zoh0, zoherr0 = ZOH_M13(tab0, ext='_rg', name='ZOH_rg', err=True)
            zoh1, zoherr1 = ZOH_M13(tab1, ext='_sm', name='ZOH_sm', err=True)
            tab0.add_column(zoh0)
            tab0.add_column(zoherr0)
            tab1.add_column(zoh1)
            tab1.add_column(zoherr1)
            BPT0 = bpt_type(tab0, ext='_rg', name='BPT_rg')
            BPT1 = bpt_type(tab1, ext='_sm', name='BPT_sm')
            tab0.add_column(BPT0)
            tab1.add_column(BPT1)
        elif prod == 'SSP':
            # For stellar surface density we need distance
            star0 = stmass_pc2(tab0['mass_ssp_rg'], 
                            dist=dist.loc[gal]['caDistP3d'], name='sigstar_rg')
            star1 = stmass_pc2(tab1['mass_ssp_sm'], 
                            dist=dist.loc[gal]['caDistP3d'], name='sigstar_sm')
            avstar0 = stmass_pc2(tab0['mass_Avcor_ssp_rg'], 
                            dist=dist.loc[gal]['caDistP3d'], name='sigstar_Avcor_rg')
            avstar1 = stmass_pc2(tab1['mass_Avcor_ssp_sm'], 
                            dist=dist.loc[gal]['caDistP3d'], name='sigstar_Avcor_sm')
            tab0.add_column(star0)
            tab1.add_column(star1)
            tab0.add_column(avstar0)
            tab1.add_column(avstar1)

    if len(rglist) > 0:
        rg_merge = vstack(rglist)
    rg_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    print(rg_merge.colnames)
    print('There are',len(rg_merge),'rows in native table')

    if len(smlist) > 0:
        sm_merge = vstack(smlist)
    sm_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    print(sm_merge.colnames)
    print('There are',len(sm_merge),'rows in smoothed table')

    if (len(filelist) > 1):
        outname = 'edge'
    else:
        outname = gal
    if prod == prodtype[0]:
        rg_merge.write(outname+'.pipe3d.hdf5', path=prod+'_rg', overwrite=True, 
                serialize_meta=True, compression=True)
    else:
        rg_merge.write(outname+'.pipe3d.hdf5', path=prod+'_rg', append=True, 
                serialize_meta=True, compression=True)
    sm_merge.write(outname+'.pipe3d.hdf5', path=prod+'_sm', append=True, 
                serialize_meta=True, compression=True)

