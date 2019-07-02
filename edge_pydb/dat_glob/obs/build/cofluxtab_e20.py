#!/usr/bin/env python

# Makes a table of CO & 13CO flux information
# Run this in the directory where the IDL moments script was run.

import os.path
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.io import ascii
from astropy import units as u

def getvwidth(fname=None):
    data = ascii.read(fname)
    vel    = data['col1']
    flux   = data['col2']
    idx = np.where(flux != 0.00)[0]
    if len(idx) > 0:
        vmin = vel[idx[0]]
        vmax = vel[idx[-1]]
        vwidth = abs(vmax-vmin)
    else:
        vwidth = float('NaN')
    print('The velocity width is {0}'.format(vwidth))
    return vwidth


linelist = ['co', '13co']
lbl      = ['co', 'cott']
base     = 'e20'
altbase  = 'e20'

# Initialize the table
tab=Table(names=['coRactr_'+base,'coDectr_'+base,'coCtrint_'+base
    ,'coBmaj_'+base,'coBmin_'+base,'coBpa_'+base
    ,'coNomask_'+base,'coeNomask_'+base,'coNomaskDv_'+base
    ,'coDilated_'+base,'coeDilated_'+base
    ,'coSmooth_'+base,'coeSmooth_'+base,'coSmoothDv_'+base
    ,'coMask2d_'+base,'coeMask2d_'+base
    ,'coSNRmax_'+base,'coSNR4pix_'+base,'coSNR5pix_'+base
    ,'cottBmaj_'+base,'cottBmin_'+base,'cottBpa_'+base
    ,'cottNomask_'+base,'cotteNomask_'+base,'cottNomaskDv_'+base
    ,'cottDilated_'+base,'cotteDilated_'+base
    ,'cottSmooth_'+base,'cotteSmooth_'+base,'cottSmoothDv_'+base
    ,'cottMask2d_'+base,'cotteMask2d_'+base
    ,'cottSNRmax_'+base,'cottSNR4pix_'+base,'cottSNR5pix_'+base])
col=Column(['UGC05498NED01'],name='Name',dtype='str',description='Galaxy Name')
tab.add_row()
tab.add_column(col, index=0)
tab['coRactr_'+base].unit = 'deg'
tab['coRactr_'+base].description = 'Reference R.A. of '+base+' CARMA cube'
tab['coDectr_'+base].unit = 'deg'
tab['coDectr_'+base].description = 'Reference Dec. of '+base+' CARMA cube'
tab['coCtrint_'+base].description = 'Unmasked CO intensity at reference pixel'
for i, line in enumerate(linelist):
    tab[lbl[i]+'Bmaj_'+base].unit = 'arcsec'
    tab[lbl[i]+'Bmaj_'+base].description = 'Beam major axis of '+base+' '+line+' cube'
    tab[lbl[i]+'Bmin_'+base].unit = 'arcsec'
    tab[lbl[i]+'Bmin_'+base].description = 'Beam minor axis of '+base+' '+line+' cube'
    tab[lbl[i]+'Bpa_'+base].unit = 'deg'
    tab[lbl[i]+'Bpa_'+base].description = 'Beam pos ang (deg E of N) of '+base+' '+line+' cube'
    tab[lbl[i]+'Nomask_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'Nomask_'+base].description = line+' flux from unmasked '+base+' cube'
    tab[lbl[i]+'eNomask_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'eNomask_'+base].description = line+' flux uncertainty from unmasked '+base+' cube'
    tab[lbl[i]+'NomaskDv_'+base].unit = 'km / s'
    tab[lbl[i]+'NomaskDv_'+base].description = line+' velocity width for unmasked '+base+' cube'
    tab[lbl[i]+'Dilated_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'Dilated_'+base].description = line+' flux from dilated-masked '+base+' cube'
    tab[lbl[i]+'eDilated_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'eDilated_'+base].description = line+' flux uncertainty from dilated-masked '+base+' cube'
    tab[lbl[i]+'Smooth_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'Smooth_'+base].description = line+' flux from smoothed-masked '+base+' cube'
    tab[lbl[i]+'eSmooth_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'eSmooth_'+base].description = line+' flux uncertainty from smoothed-masked '+base+' cube'
    tab[lbl[i]+'SmoothDv_'+base].unit = 'km / s'
    tab[lbl[i]+'SmoothDv_'+base].description = line+' velocity width for smooth-masked '+base+' cube'
    tab[lbl[i]+'Mask2d_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'Mask2d_'+base].description = line+' flux from 2D masked '+base+' cube'
    tab[lbl[i]+'eMask2d_'+base].unit = 'Jy km / s'
    tab[lbl[i]+'eMask2d_'+base].description = line+' flux uncertainty from 2D masked '+base+' cube'
    tab[lbl[i]+'SNRmax_'+base].description = 'Maximum SNR of '+line+' in '+base+' cube'
    tab[lbl[i]+'SNR4pix_'+base].description = 'Number of XY pixels with '+line+' peak Tb > 4 sigma'
    tab[lbl[i]+'SNR5pix_'+base].description = 'Number of XY pixels with '+line+' peak Tb > 5 sigma'

# Get the parameters
flist = glob.glob('*.'+line+'.'+altbase+'_str.flux.out')
for j, file in enumerate(flist):
    gal = file.split('.')[0]
    if j > 0:
        tab.add_row()
    tab['Name'][j] = gal
    for col in tab.colnames:
        if col != 'Name':
            tab[col][j] = float('NaN')
    print(file)
    for i, line in enumerate(linelist):
        if os.path.isfile(gal+'.'+line+'.'+altbase+'_str.mom0.fits') == False:
            continue
        hdu = fits.open(gal+'.'+line+'.'+altbase+'_str.mom0.fits',ignore_missing_end=True)
        img = hdu[0].data
        hd  = hdu[0].header
        if line == 'co':
            tab['coRactr_'+base][j]  = hd['CRVAL1']
            tab['coDectr_'+base][j]  = hd['CRVAL2']
            tab['coCtrint_'+base][j] = img[int(hd['CRPIX2']-1),int(hd['CRPIX1']-1)]
        if 'K.KM/S' in hd['BUNIT']:
            aspp = abs(hd['CDELT2'])*3600. * u.arcsec
            freq = hd['RESTFREQ'] * u.Hz
            convfac = (u.K).to(u.Jy, equivalencies=u.brightness_temperature(aspp**2,freq))
            print("The conversion factor to Jy/pix is %s" % convfac)
        else:
            convfac = 1.
        tab[lbl[i]+'Bmaj_'+base][j] = hd['BMAJ']*3600.
        tab[lbl[i]+'Bmin_'+base][j] = hd['BMIN']*3600.
        if 'BPA' in hd.keys():
            tab[lbl[i]+'Bpa_'+base][j] = hd['BPA']
        else:
            tab[lbl[i]+'Bpa_'+base][j] = 0.
        # --- Get the fluxes
        dat = np.genfromtxt(gal+'.'+line+'.'+altbase+'_str.flux.out',max_rows=1,
            comments=';',usecols=[1,4])
        tab[lbl[i]+'Nomask_'+base][j]   = dat[0] * convfac
        tab[lbl[i]+'eNomask_'+base][j]  = dat[1] * convfac
        tab[lbl[i]+'NomaskDv_'+base][j] = getvwidth(gal+'.'+line+'.'+altbase+'_str.flux.out')
        dat = np.genfromtxt(gal+'.'+line+'.'+altbase+'_dil.flux.out',max_rows=1,
            comments=';',usecols=[1,4])
        if dat[0] != 0.0 and dat[1] != 0.0:
            tab[lbl[i]+'Dilated_'+base][j] = dat[0] * convfac
            tab[lbl[i]+'eDilated_'+base][j] = dat[1] * convfac
        dat = np.genfromtxt(gal+'.'+line+'.'+altbase+'_smo.flux.out',max_rows=1,
            comments=';',usecols=[1,4])
        if dat[0] != 0.0 and dat[1] != 0.0:
            tab[lbl[i]+'Smooth_'+base][j] = dat[0] * convfac
            tab[lbl[i]+'eSmooth_'+base][j] = dat[1] * convfac
        tab[lbl[i]+'SmoothDv_'+base][j] = getvwidth(gal+'.'+line+'.'+altbase+'_smo.flux.out')
        if os.path.isfile(gal+'.'+line+'.'+altbase+'_mk2.flux.out'):
            dat = np.genfromtxt(gal+'.'+line+'.'+altbase+'_mk2.flux.out',max_rows=1,
                comments=';',usecols=[1,4])
            if dat[0] != 0.0 and dat[1] != 0.0:
                tab[lbl[i]+'Mask2d_'+base][j] = dat[0] * convfac
                tab[lbl[i]+'eMask2d_'+base][j] = dat[1] * convfac
        if os.path.isfile(gal+'.'+line+'.'+altbase+'_dil.snrpk.fits'):
            hdulist = fits.open(gal+'.'+line+'.'+altbase+'_dil.snrpk.fits')
            snrdata = hdulist[0].data
            tab[lbl[i]+'SNRmax_'+base][j] = np.nanmax(snrdata)
            good=np.isfinite(snrdata)
            tab[lbl[i]+'SNR4pix_'+base][j] = np.size(np.where(snrdata[good]>4))
            tab[lbl[i]+'SNR5pix_'+base][j] = np.size(np.where(snrdata[good]>5))

# Only need 3 sigfigs for most parameters
for cname in tab.colnames:
    if cname == 'coRactr_'+base or cname == 'coDectr_'+base:
        tab[cname].format='.4f'
    elif cname != 'Name':
        tab[cname].format='.3f'

# Take this from the last galaxy processed
bunit = hd['BUNIT'].replace('KM/S','km / s')
tab['coCtrint_'+base].unit = u.Unit(bunit,format="fits")

# Rename these two galaxies
tab['Name'][np.where(tab['Name']=='NGC4211N')] = 'NGC4211NED02'
tab['Name'][np.where(tab['Name']=='UGC05498')] = 'UGC05498NED01'

# Write the table
tab.write('edge_'+base+'_coflux.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

