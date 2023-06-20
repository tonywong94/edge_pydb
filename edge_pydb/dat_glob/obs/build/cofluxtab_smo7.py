#!/usr/bin/env python

# Makes tables of CO & 13CO flux information and spectra
# Run this in the directory where the moments script was run.
# 16mar2023 - get vwidth for dilated mask

from datetime import datetime
import os.path
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join, vstack
from astropy.io import ascii
from astropy import units as u

fdir = '../'   # where FITS files reside
linelist = ['co', '13co']
lbl      = ['co', 'cott']
base     = 'smo7'
altbase  = 'smo7'
matchres = True

def getvwidth(spectbl=None):
    vel    = spectbl['Velocity']
    flux   = spectbl['Flux']
    totwidth = abs(vel[-1]-vel[0])
    idx = np.where(flux != 0.00)[0]
    if len(idx) > 0:
        vmin = vel[idx[0]]
        vmax = vel[idx[-1]]
        mskwidth = abs(vmax-vmin)
    else:
        mskwidth = np.nan
    #print('Masked velocity width is {0}'.format(mskwidth))
    return totwidth, mskwidth

# Initialize the global table
gnames = ['coRactr_'+base,'coDectr_'+base,'coCtrint_'+base,
          'coNomask_'+base,'coeNomask_'+base,'coNomaskDv_'+base,
          'coDilated_'+base,'coeDilated_'+base,'coDilatedDv_'+base,
          'coSmooth_'+base,'coeSmooth_'+base,'coSmoothDv_'+base,
          'coMask2d_'+base,'coeMask2d_'+base,
          'coSNRmax_'+base,'coSNR4pix_'+base,'coSNR5pix_'+base]
if len(linelist)>1:
    gnames += ['cottNomask_'+base,'cotteNomask_'+base,'cottNomaskDv_'+base,
               'cottDilated_'+base,'cotteDilated_'+base,'cottDilatedDv_'+base,
               'cottSmooth_'+base,'cotteSmooth_'+base,'cottSmoothDv_'+base,
               'cottMask2d_'+base,'cotteMask2d_'+base]
    if matchres:
        gnames += ['cottMk12_di_'+base,'cotteMk12_di_'+base,
                   'cottMk12_sm_'+base,'cotteMk12_sm_'+base]
    gnames += ['cottSNRmax_'+base,'cottSNR4pix_'+base,'cottSNR5pix_'+base]
gtab=Table(names=gnames)
namecol=Column(['UGC05498NED01'],name='Name',dtype='str',description='Galaxy Name')
gtab.add_row()
gtab.add_column(namecol, index=0)


# Metadata for global table
gtab['coRactr_'+base].unit = 'deg'
gtab['coRactr_'+base].description = 'Reference R.A. of '+base+' CARMA cube'
gtab['coDectr_'+base].unit = 'deg'
gtab['coDectr_'+base].description = 'Reference Dec. of '+base+' CARMA cube'
gtab['coCtrint_'+base].description = 'Unmasked CO intensity at reference pixel'
for i, line in enumerate(linelist):
    gtab[lbl[i]+'Nomask_'+base].description = line+' flux from unmasked '+base+' cube'
    gtab[lbl[i]+'eNomask_'+base].description = line+' flux uncertainty from unmasked '+base+' cube'
    gtab[lbl[i]+'NomaskDv_'+base].unit = 'km / s'
    gtab[lbl[i]+'NomaskDv_'+base].description = line+' velocity width for unmasked '+base+' cube'
    gtab[lbl[i]+'Dilated_'+base].description = line+' flux from dilated-masked '+base+' cube'
    gtab[lbl[i]+'eDilated_'+base].description = line+' flux uncertainty from dilated-masked '+base+' cube'
    gtab[lbl[i]+'DilatedDv_'+base].unit = 'km / s'
    gtab[lbl[i]+'DilatedDv_'+base].description = line+' velocity width for dilated-masked '+base+' cube'
    gtab[lbl[i]+'Smooth_'+base].description = line+' flux from smoothed-masked '+base+' cube'
    gtab[lbl[i]+'eSmooth_'+base].description = line+' flux uncertainty from smoothed-masked '+base+' cube'
    gtab[lbl[i]+'SmoothDv_'+base].unit = 'km / s'
    gtab[lbl[i]+'SmoothDv_'+base].description = line+' velocity width for smooth-masked '+base+' cube'
    gtab[lbl[i]+'Mask2d_'+base].description = line+' flux from 2D masked '+base+' cube'
    gtab[lbl[i]+'eMask2d_'+base].description = line+' flux uncertainty from 2D masked '+base+' cube'
    if line == '13co' and matchres:
        gtab[lbl[i]+'Mk12_di_'+base].description = line+' flux in 12co dilated mask '+base+' cube'
        gtab[lbl[i]+'eMk12_di_'+base].description = line+' flux uncertainty in 12co dilated mask '+base+' cube'
        gtab[lbl[i]+'Mk12_sm_'+base].description = line+' flux in 12co smoothed mask '+base+' cube'
        gtab[lbl[i]+'eMk12_sm_'+base].description = line+' flux uncertainty in 12co smoothed mask '+base+' cube'
    gtab[lbl[i]+'SNRmax_'+base].description = 'Maximum SNR of '+line+' in '+base+' cube'
    gtab[lbl[i]+'SNR4pix_'+base].description = 'Number of XY pixels with '+line+' peak Tb > 4 sigma'
    gtab[lbl[i]+'SNR5pix_'+base].description = 'Number of XY pixels with '+line+' peak Tb > 5 sigma'



# Get the parameters
flist = sorted(glob.glob('*.'+linelist[0]+'.'+altbase+'_dil.flux.csv'))
tablelist=[]
for j, file in enumerate(flist):
    gal = file.split('.')[0]
    if j > 0:
        gtab.add_row()
    gtab['Name'][j] = gal
    for col in gtab.colnames:
        if col != 'Name':
            gtab[col][j] = float('NaN')
    print('Working on galaxy',gal,'lines',linelist)
    for i, line in enumerate(linelist):
        if not os.path.isfile(fdir+gal+'.'+line+'.'+altbase+'_str.mom0.fits.gz'):
            print('Missing',fdir+gal+'.'+line+'.'+altbase+'_str.mom0.fits.gz')
            continue
        hdu = fits.open(fdir+gal+'.'+line+'.'+altbase+'_str.mom0.fits.gz')
        img = hdu[0].data
        hd  = hdu[0].header
        if line == 'co':
            gtab['coRactr_'+base][j]  = hd['CRVAL1']
            gtab['coDectr_'+base][j]  = hd['CRVAL2']
            gtab['coCtrint_'+base][j] = img[int(hd['CRPIX2']-1),int(hd['CRPIX1']-1)]
        # --- No mask
        dat = Table.read(gal+'.'+line+'.'+altbase+'_str.flux.csv',format='ascii.ecsv')
        totflux = dat.meta['totflux'].split()
        gtab[lbl[i]+'Nomask_'+base][j]  = totflux[0]
        gtab[lbl[i]+'eNomask_'+base][j] = totflux[2]
        if j == 0:
            gtab[lbl[i]+'Nomask_'+base].unit  = " ".join(totflux[3:])
            gtab[lbl[i]+'eNomask_'+base].unit = " ".join(totflux[3:])
        newt = Table(dat, names=['Vel',
                                 lbl[i]+'NomaskSpec_'+base,
                                 lbl[i]+'NomaskUnc_'+base])
        del newt.meta['totflux']
        newt.columns[1].description = 'Unmasked spectrum'
        newt.columns[2].description = 'Error in Unmasked spectrum'
        if line == 'co':
            spec = newt
        else:
            spec = join(spec, newt, keys='Vel', join_type='outer')
        # --- Dilated mask
        dat = Table.read(gal+'.'+line+'.'+altbase+'_dil.flux.csv',format='ascii.ecsv')
        totflux = dat.meta['totflux'].split()
        gtab[lbl[i]+'Dilated_'+base][j]  = totflux[0]
        gtab[lbl[i]+'eDilated_'+base][j] = totflux[2]
        gtab[lbl[i]+'NomaskDv_'+base][j], gtab[lbl[i]+'DilatedDv_'+base][j] = getvwidth(dat)
        if j == 0:
            gtab[lbl[i]+'Dilated_'+base].unit  = " ".join(totflux[3:])
            gtab[lbl[i]+'eDilated_'+base].unit = " ".join(totflux[3:])
        newt = Table(dat, names=['Vel',
                                 lbl[i]+'DilatedSpec_'+base,
                                 lbl[i]+'DilatedUnc_'+base])
        del newt.meta['totflux']
        newt.columns[1].description = 'Dilated mask spectrum'
        newt.columns[2].description = 'Error in Dilated mask spectrum'
        spec = join(spec, newt, keys='Vel', join_type='outer')
        # --- Smooth mask
        dat = Table.read(gal+'.'+line+'.'+altbase+'_smo.flux.csv',format='ascii.ecsv')
        totflux = dat.meta['totflux'].split()
        gtab[lbl[i]+'Smooth_'+base][j]  = totflux[0]
        gtab[lbl[i]+'eSmooth_'+base][j] = totflux[2]
        junk, gtab[lbl[i]+'SmoothDv_'+base][j] = getvwidth(dat)
        if j == 0:
            gtab[lbl[i]+'Smooth_'+base].unit  = " ".join(totflux[3:])
            gtab[lbl[i]+'eSmooth_'+base].unit = " ".join(totflux[3:])
        tot2dflux = dat.meta['tot2dflux'].split()
        gtab[lbl[i]+'Mask2d_'+base][j]  = tot2dflux[0]
        gtab[lbl[i]+'eMask2d_'+base][j] = tot2dflux[2]
        if j == 0:
            gtab[lbl[i]+'Mask2d_'+base].unit  = " ".join(tot2dflux[3:])
            gtab[lbl[i]+'eMask2d_'+base].unit = " ".join(tot2dflux[3:])
        newt = Table(dat, names=['Vel',
                                 lbl[i]+'SmoothSpec_'+base,
                                 lbl[i]+'SmoothUnc_'+base,
                                 lbl[i]+'Mask2dSpec_'+base,
                                 lbl[i]+'Mask2dUnc_'+base])
        del newt.meta['totflux']
        del newt.meta['tot2dflux']
        newt.columns[1].description = 'Smoothed mask spectrum'
        newt.columns[2].description = 'Error in Smoothed mask spectrum'
        newt.columns[3].description = '2D mask spectrum'
        newt.columns[4].description = 'Error in 2D mask spectrum'
        spec = join(spec, newt, keys='Vel', join_type='outer')
        # --- 13CO masks from 12CO
        if line == '13co' and matchres:
            if os.path.isfile(gal+'.'+line+'.'+altbase+'_mk12_dil.flux.csv'):
                dat = Table.read(gal+'.'+line+'.'+altbase+'_mk12_dil.flux.csv',
                                 format='ascii.ecsv')
                totflux = dat.meta['totflux'].split()
                gtab[lbl[i]+'Mk12_di_'+base][j]  = totflux[0]
                gtab[lbl[i]+'eMk12_di_'+base][j] = totflux[2]
                if j == 0:
                    gtab[lbl[i]+'Mk12_di_'+base].unit  = " ".join(totflux[3:])
                    gtab[lbl[i]+'eMk12_di_'+base].unit = " ".join(totflux[3:])
                newt = Table(dat, names=['Vel',
                                         lbl[i]+'Mk12_diSpec_'+base,
                                         lbl[i]+'Mk12_diUnc_'+base])
                del newt.meta['totflux']
                newt.columns[1].description = '12CO-based Dilated mask spectrum'
                newt.columns[2].description = 'Error in 12CO-based Dilated mask spectrum'
                spec = join(spec, newt, keys='Vel', join_type='outer')
            if os.path.isfile(gal+'.'+line+'.'+altbase+'_mk12_smo.flux.csv'):
                dat = Table.read(gal+'.'+line+'.'+altbase+'_mk12_smo.flux.csv',
                                 format='ascii.ecsv')
                totflux = dat.meta['totflux'].split()
                gtab[lbl[i]+'Mk12_sm_'+base][j]  = totflux[0]
                gtab[lbl[i]+'eMk12_sm_'+base][j] = totflux[2]
                if j == 0:
                    gtab[lbl[i]+'Mk12_sm_'+base].unit  = " ".join(totflux[3:])
                    gtab[lbl[i]+'eMk12_sm_'+base].unit = " ".join(totflux[3:])
                newt = Table(dat, names=['Vel',
                                         lbl[i]+'Mk12_smSpec_'+base,
                                         lbl[i]+'Mk12_smUnc_'+base])
                del newt.meta['totflux']
                newt.columns[1].description = '12CO-based Smoothed mask spectrum'
                newt.columns[2].description = 'Error in 12CO-based Smoothed mask spectrum'
                spec = join(spec, newt, keys='Vel', join_type='outer')
        # -- Peak SNR statistics
        if os.path.isfile(fdir+gal+'.'+line+'.'+altbase+'_dil.snrpk.fits.gz'):
            hdulist = fits.open(fdir+gal+'.'+line+'.'+altbase+'_dil.snrpk.fits.gz')
            snrdata = hdulist[0].data
            gtab[lbl[i]+'SNRmax_'+base][j] = np.nanmax(snrdata)
            good=np.isfinite(snrdata)
            gtab[lbl[i]+'SNR4pix_'+base][j] = np.size(np.where(snrdata[good]>4))
            gtab[lbl[i]+'SNR5pix_'+base][j] = np.size(np.where(snrdata[good]>5))
    nrows = len(spec)
    gname = Column([gal]*nrows, name='Name', description='Galaxy Name')
    spec.add_column(gname, index=0)
    tablelist.append(spec)

# Only need 3 sigfigs for most parameters
for cname in gtab.colnames:
    if cname == 'coRactr_'+base or cname == 'coDectr_'+base:
        gtab[cname].info.format='.4f'
    elif cname != 'Name':
        gtab[cname].info.format='.3f'
    if base in cname and 'SNR' not in cname:
        gtab[cname][np.where(gtab[cname] == 0.)] = np.nan

# Take this from the last galaxy processed
bunit = hd['BUNIT']
gtab['coCtrint_'+base].unit = u.Unit(bunit,format="fits")

# Write the global table
gtab.meta['date'] = datetime.today().strftime('%Y-%m-%d')
gtab.meta['comments'] = ('Integrated CO fluxes from '+altbase+' cubes')
gtab.write('edge_coflux_'+base+'.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# Write the spectral table
spectab = vstack(tablelist)
spectab['Vel'].name = 'coVlsr_'+base
spectab['coVlsr_'+base].unit = 'km / s'
for cname in spectab.colnames[2:]:
    spectab[cname][np.where(spectab[cname] == 0.)] = np.nan
    spectab[cname].fill_value = np.nan
    spectab[cname].format='.5f'
spectab.meta['date'] = datetime.today().strftime('%Y-%m-%d')
spectab.meta['comments'] = ('Integrated CO spectra from '+altbase+' cubes')
spectab.write('edge_cospec_'+base+'.csv', format='ascii.ecsv', delimiter=',', overwrite=True)
