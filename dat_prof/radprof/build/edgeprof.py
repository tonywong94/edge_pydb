#!/usr/bin/env python

import os
import numpy as np
from astropy import units as u
from astropy.table import Table, Column
from astropy.io import fits 
from astropy.wcs import WCS
from ellprof import ellprof

def edgeprof(datdir='.', gname=None, line='co', type='de20_smo', delr=3, 
             dbdir=None, alpha=4.3, fmask=0.2, blanksig=0, deproj=True, 
             replace=False, outdir='rproftxt', radimg=False):

    """
    Derive the radial profile of an EDGE galaxy from a mom-0 image. Geometric
    parameters are derived from LEDA and ringfit; distance is from CALIFA.
    Detection limit profile is taken as 3 times the larger of emom0 and
    emom0max, the latter being the noise map integrated over 200 km/s.  Even if
    the mom0 image is empty, it is assumed that emom0max is available. Need
    edge-pydb installed from https://github.com/tonywong94/edge_pydb.

    INPUTS:
    datdir   - directory where moment map fits files reside
    gname    - Galaxy name
    line     - line code in filename, e.g. 'co'
    type     - type code in filename, e.g. 'de20_smo'
    dbdir    - top level directory of EDGE database
    delr     - ring spacing in arcsec (default 3 arcsec)
    alpha    - conversion from K km/s to Msol/pc^2 (default 4.3 for Galactic CO)
    fmask    - minimum fraction of non-NaN pixels in a ring (default 0.2)
    deproj   - True to deproject surface densities & noises (default=True)
    blanksig - constant to multiply RMS image for blank replacement (default 0)
    replace  - Replace values < blanksig*RMS with blanksig*RMS? (default=False)
    radimg   - True to write out radius image (default=False)

    OUTPUTS: radial profile table (ECSV format)
    """

    # Handle the inputs
    alphaco      = alpha * u.solMass * u.s / (u.K * u.km * u.pc**2)
    dr_arcsec    = delr * u.arcsec
    mom0file     = datdir+'/'+gname+'.'+line+'.'+type+'.mom0.fits'
    emom0file    = datdir+'/'+gname+'.'+line+'.'+type+'.emom0.fits'
    emom0maxbase = gname+'.'+line+'.'+type+'.emom0max.fits'
    emom0maxfile = datdir+'/'+emom0maxbase

    # Make the output directory if required
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Center position and optical radius from LEDA
    leda   = Table.read(dbdir+'/dat_glob/external/edge_leda.csv', 
                    format='ascii.ecsv')
    idx    = np.where(leda['Name']==gname)[0][0]
    ractr  = leda['ledaRA'].quantity[idx].to(u.deg)
    dcctr  = leda['ledaDE'].quantity[idx].to(u.deg)
    r25    = leda['ledaD25'].quantity[idx].to(u.arcsec) / 2.
    
    # PA and INC from Becca's table
    rfpars = Table.read(dbdir+'/dat_glob/derived/edge_rfpars.csv',
                    format='ascii.ecsv')
    idx    = np.where(rfpars['Name']==gname)[0][0]
    pa     = rfpars['rfPA'][idx]
    inc    = rfpars['rfInc'][idx]
    
    # Distance from caZgas assuming Ho=70, Om=0.27, Ol=0.73
    califa = Table.read(dbdir+'/dat_glob/external/edge_califa.csv', 
                    format='ascii.ecsv')
    idx    = np.where(califa['Name']==gname)[0][0]
    dmpc   = califa['caDistMpc'].quantity[idx]

    # Read in the error maps
    edata2, ehdr = fits.getdata(emom0maxfile, header=True)
    if os.path.exists(emom0file):
        edata1 = fits.getdata(emom0file, header=False)
        # Choose larger of two estimates for each pixel
        edata = np.fmax(edata1, edata2)
    else:
        edata = edata2

    # Get the center position and bin size in pixels
    w = WCS(ehdr)
    x,y = w.wcs_world2pix(ractr,dcctr,0)
    xoff = x - (ehdr['crpix1']-1)
    yoff = y - (ehdr['crpix2']-1)
    center = np.array([x, y])
    aspp = round(ehdr['cdelt2'] * 3600., 3) * u.arcsec
    dr_pix = dr_arcsec/aspp
    dx = (aspp * dmpc).to(u.pc, equivalencies=u.dimensionless_angles())
    print('Radial profile for galaxy {}'.format(gname))
    print('Centered on XOFF={:.2f} YOFF={:.2f}'.format(xoff,yoff))
    print('Using a bin size of {0} pixels'.format(dr_pix))
    print('1 pixel is {0} or {1:.2f}'.format(aspp,dx))

    # Estimate the detection limit from emom0max map
    etab, rad = ellprof(edata**2, pa=pa, incline=inc, center=center, 
        binsize=dr_pix, nanmask=True)

    # Read in the moment-0 map, if present
    if os.path.exists(mom0file):
        data, hdr = fits.getdata(mom0file, header=True)
        if replace==True and os.path.exists(emom0file):
            small = np.where(data / edata1 < blanksig)
            data[small] = blanksig * edata1[small]
            #nsmall = np.size(data[data < blanksig*edata])
            #data[data < blanksig*edata] = blanksig*edata
            print('Replaced {0} values less than {1} sigma'.format(
                    len(small[0]),blanksig))
    else:
        # Otherwise feed an array of pure NaNs
        data = edata * np.nan

    # Calculate the radial profile, replacing NaNs with scaled uncertainty
    if blanksig > 0.:
        print('Replacing blanks with {0} sigma'.format(blanksig))
    tab, rad = ellprof(data, pa=pa, incline=inc, center=center, binsize=dr_pix, 
                       nanmask=False, blank=blanksig*edata, fracmask=fmask)

    # Output the radius image if requested
    if radimg:
        hdr['datamin'] = np.nanmin(rad)
        hdr['datamax'] = np.nanmax(rad)
        fits.writeto(outdir+'/'+gname+'.radius.fits.gz', rad, hdr, 
                     overwrite=True)
    
    # Do unit conversions to make master table
    freq = ehdr['restfreq'] * u.Hz
    bmaj = ehdr['bmaj']*3600. * u.arcsec
    bmin = ehdr['bmin']*3600. * u.arcsec
    omega_B = (bmaj*bmin)*2*np.pi/(8*np.log(2))
    ppbeam = omega_B/aspp**2
    #print("\nPixels per beam: {0}".format(ppbeam))
    
    # Propagation of errors to get detection limit
    sigma = np.sqrt(etab['annsum']*ppbeam)/etab['npix']
    tab['detlim'] = 3*sigma
    
    # This assumes units are either Jy/bm*km/s or K*km/s
    if 'JY/BEAM' in ehdr['bunit'].upper():
        convfac = u.Jy.to(u.K, equivalencies=
            u.brightness_temperature(omega_B, freq))
        print('Multiplying by {0} to convert to K'.format(convfac))
        tab['wtmean'] *= convfac
        tab['rms']    *= convfac
        tab['detlim'] *= convfac
        tab['annsum'] *= convfac
        tab['cumsum'] *= convfac
    elif ehdr['bunit'].upper()[0] != 'K':
        print("\nWarning: Unrecognized brightness unit")

    # Alternative radii
    tab['radius'] *= aspp.value
    tab['radius'].unit = 'arcsec'
    rkpc = (tab['radius'] * dmpc).to(u.kpc, equivalencies=u.dimensionless_angles())
    newcol = Column(rkpc, name='r_kpc')
    tab.add_column(newcol, index=1)
    rr25 = tab['radius'] / r25
    newcol = Column(rr25, name='r_r25')
    tab.add_column(newcol, index=2)

    # Set units and deproject
    tab['wtmean'].unit = u.K * u.km / u.s
    tab['rms'].unit    = u.K * u.km / u.s
    tab['detlim'].unit = u.K * u.km / u.s
    tab['annsum'].unit = tab['wtmean'].unit * u.pix
    tab['cumsum'].unit = tab['wtmean'].unit * u.pix
    if deproj==True:
        projfac = np.cos(np.radians(inc))
        print('Multiplying sigmol by {0:.2f} to project to face-on'.format(projfac))
        tab['wtmean'] *= projfac
        tab['rms']    *= projfac
        tab['detlim'] *= projfac
    else:
        projfac = 1.
    tab['cumlum']  = tab['cumsum'] * (dx**2/u.pix)
    tab['sigmol']  = alphaco * tab['wtmean']
    tab['cummass'] = alphaco * tab['cumlum']

    # Make the table informative
    for field in ['r_kpc', 'r_r25', 'wtmean', 'annsum', 'detlim', 'sigmol']:
        tab[field].format = '10.3f'
    for field in ['cumlum', 'cummass']:
        tab[field].format = '10.4e'
    tab['radius'].description = (
        'Adopting RA={0:.4f} DEC={1:.4f} PA={2:.1f} INC={3:.1f}'.format(
        ractr,dcctr,pa,inc))
    tab['r_kpc'].description = 'Adopting distance of {0:.2f}'.format(dmpc)
    tab['r_r25'].description = 'Adopting R25={0:.2f}'.format(r25)
    tab['ngood'].description = 'Requiring fill fraction of {0} npix'.format(fmask)
    tab['wtmean'].description = 'Replaced blanks with {0} times {1}'.format(
        blanksig,emom0maxbase)
    tab.write(outdir+'/'+gname+'.'+type+'b'+str(blanksig)+'.rprof.txt', 
        delimiter=',', format='ascii.ecsv', overwrite=True)
    return

