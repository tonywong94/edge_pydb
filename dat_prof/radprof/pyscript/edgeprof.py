#!/usr/bin/env python

from ellprof import ellprof
import numpy as np
from astropy import units as u
from astropy.table import Table, Column
from astropy.io import fits 
from astropy.wcs import WCS
from matplotlib import pyplot as plt
import os

def edgeprof(indir='.', gname=None, line='co', type='de20_smo', delr=3, db=None, 
    alpha=4.3, fmask=0.2, blanksig=0, deproj=True, ndblank=False):

    """
    Derive the radial profile of an EDGE galaxy using the smooth-mask mom-0.
    Detection limit profile is based on noise map integrated over 200 km/s.
    Need edge-sql-base installed from github.com/tonywong94/edge-sql-base

    INPUTS:
    indir    - directory where moment map fits files reside
    gname    - Galaxy name
    line     - line code in filename, e.g. 'co'
    type     - type code in filename, e.g. 'de20_smo'
    db       - LEDA database to get galaxy parameters
    delr     - ring spacing in arcsec
    alpha    - conversion from K km/s to Msol/pc^2 (default 4.3 for Galactic CO)
    fmask    - minimum fraction of non-NaN pixels required in a ring
    blanksig - constant to multiply RMS image for blank replacement (default=0)
    deproj   - True to deproject the mass surface density (default=True)
    ndblank  - True to replace small values with blanks (default=False)

    OUTPUTS: radial profile table (ECSV format)
    """

    alphaco = alpha * u.solMass * u.s / (u.K * u.km * u.pc**2)
    dr_arcsec = delr * u.arcsec
    mom0file = indir+'/'+gname+'.'+line+'.'+type+'.mom0.fits'
    emom0file = indir+'/'+gname+'.'+line+'.'+type+'.emom0.fits'
    emom0maxfile = indir+'/'+gname+'.'+line+'.'+type+'.emom0max.fits'

    # Center position and optical radius from LEDA
    idx   = np.where(db['Name']==gname)[0][0]
    ractr = db['ledaRA'].quantity[idx].to(u.deg)
    dcctr = db['ledaDE'].quantity[idx].to(u.deg)
    r25   = db['ledaD25'].quantity[idx].to(u.arcsec) / 2.
    # PA and INC from Becca's table
    pa    = db['coPA'][idx]
    inc   = db['coInc'][idx]
    # Distance from caZgas assuming Ho=70, Om=0.27, Ol=0.73
    dmpc  = db['caDistMpc'].quantity[idx]
    print('{6} Adopting RA={0} DEC={1} PA={2} INC={3} DIST={4} R25={5}'.format(
        ractr,dcctr,pa,inc,dmpc,r25,gname))

    # Read in the moment-0 map
    data, hdr = fits.getdata(mom0file, header=True)
    edata1 = fits.getdata(emom0file, header=False)
    if ndblank==True:
        nsmall = np.size(data[data < blanksig*edata1])
        data[data < blanksig*edata1] = float('NaN')
        print('\nBlanked {0} values less than {1} sigma'.format(nsmall,blanksig))
    edata2, ehdr = fits.getdata(emom0maxfile, header=True)
    edata = np.fmax(edata1, edata2)
    w = WCS(mom0file)
    x,y = w.wcs_world2pix(ractr,dcctr,0)
    xoff = x - (hdr['crpix1']-1)
    yoff = y - (hdr['crpix2']-1)
    center = np.array([x, y])
    aspp = round(hdr['cdelt2'] * 3600., 3) * u.arcsec
    dr_pix = dr_arcsec/aspp
    dx = (aspp * dmpc).to(u.pc, equivalencies=u.dimensionless_angles())
    print('Radial profile centered on XOFF={0} YOFF={1}'.format(xoff,yoff))
    print('Using a bin size of {0} pixels'.format(dr_pix))
    print('1 pixel is {0} or {1}'.format(aspp,dx))

    # Estimate the detection limit from emom0max map
    etab, rad = ellprof(edata**2, pa=pa, incline=inc, center=center, binsize=dr_pix, 
        nanmask=True)

    # Calculate the radial profile
    if blanksig > 0.:
        tab, rad = ellprof(data, pa=pa, incline=inc, center=center, binsize=dr_pix, 
            nanmask=False, blank=blanksig*edata, fracmask=fmask)
        print('\nReplacing blanks with {0} sigma'.format(blanksig))
    else:
        tab, rad = ellprof(data, pa=pa, incline=inc, center=center, binsize=dr_pix, 
            nanmask=False, blank=0, fracmask=fmask)

    # Output the radius image
    hdr['datamin'] = np.nanmin(rad)
    hdr['datamax'] = np.nanmax(rad)
    fits.writeto(gname+'.radius.fits.gz', rad, hdr, overwrite=True)
    
    # Do unit conversions to make master table
    freq = hdr['restfreq'] * u.Hz
    bmaj = hdr['bmaj']*3600. * u.arcsec
    bmin = hdr['bmin']*3600. * u.arcsec
    omega_B = (bmaj*bmin)*2*np.pi/(8*np.log(2))
    ppbeam = omega_B/aspp**2
    print("\nPixels per beam: {0}".format(ppbeam))
    sigma = np.sqrt(etab['annsum']*ppbeam)/etab['npix']
    tab['detlim'] = 3*sigma
    # This assumes units are either Jy/bm*km/s or K*km/s
    if 'JY/BEAM' in hdr['bunit'].upper():
        convfac = u.Jy.to(u.K, equivalencies=
            u.brightness_temperature(omega_B, freq))
        print('\nMultiplying by {0} to convert to K'.format(convfac))
        tab['wtmean'] *= convfac
        tab['rms']    *= convfac
        tab['detlim'] *= convfac
        tab['annsum'] *= convfac
        tab['cumsum'] *= convfac
    elif hdr['bunit'].upper()[0] != 'K':
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
        print('\nMultiplying sigmol by {0} to project to face-on'.format(projfac))
    else:
        projfac = 1.
    tab['cumlum']  = tab['cumsum'] * (dx**2/u.pix)
    tab['sigmol']  = alphaco * tab['wtmean'] * projfac
    tab['cummass'] = alphaco * tab['cumlum']

    # Make the table informative
    for field in ['r_kpc', 'r_r25', 'wtmean', 'annsum', 'detlim', 'sigmol']:
        tab[field].format = '10.3f'
    tab['radius'].description = 'Adopting RA={0} DEC={1} PA={2} INC={3}'.format(
        ractr,dcctr,pa,inc)
    tab['r_kpc'].description = 'Adopting distance of {0}'.format(dmpc)
    tab['r_r25'].description = 'Adopting R25={0}'.format(r25)
    tab['ngood'].description = 'Requiring fill fraction of {0} npix'.format(fmask)
    tab['wtmean'].description = 'Replaced blanks with {0} times {1}'.format(
        blanksig,emom0maxfile)
    print(tab)
    tab.write(gname+'.'+type+'b'+str(blanksig)+'.rprof.txt', 
        delimiter=',', format='ascii.ecsv', overwrite=True)
    return

# Plot radial profile for one galaxy
def pltcoprof(gal, ax, rtype='arcsec', itype='wtmean', imtype='smo', 
              ebar=True, mod0=None, modh=None):
    print('Plotting {0}...'.format(gal))
    b0file = gal+'.'+imtype+'b0.rprof.txt'
    b1file = gal+'.'+imtype+'b1.rprof.txt'
    if rtype == 'kpc':
        rmax = 15.
    elif rtype == 'r25':
        rmax = 1.05
    else:
        rmax = 60.
    if os.path.isfile(b0file):
        b0dat = Table.read(b0file, format='ascii.ecsv')
        b1dat = Table.read(b1file, format='ascii.ecsv')
        if rtype == 'kpc':
            rad = b0dat['r_kpc']
        elif rtype == 'r25':
            rad = b0dat['r_r25']
        else:
            rad = b0dat['radius']
        plt.plot(rad, b0dat[itype], color='b', ls='-', marker='^', ms=7)
        plt.plot(rad, b1dat[itype], color='g', ls='--', marker=None)
        plt.fill_between(rad,b0dat[itype], b1dat[itype], 
            facecolor='green', alpha=0.5)
        if ebar == True:
            plt.errorbar(rad, b0dat[itype], yerr=b0dat['rms'], 
                ecolor='dimgray', capsize=0, 
                zorder=1, marker=None, ls='None', lw=1, label=None)
        if (itype == 'sigmol' and rtype == 'kpc' and np.isfinite(mod0) and modh > 0):
            print('Scale length {0} with normalization {1}'.format(modh,mod0))
            imod = mod0/4.36 * np.exp(-rad/modh)
            plt.plot(rad, imod, color='m', ls='-', marker=None, lw=3)
        rcoeff = np.interp(0.5*b0dat['cummass'][-1], b0dat['cummass'], rad)
        plt.plot(rad, b0dat['detlim'], color='k', lw=2, ls='--', marker=None)
        ax.axvline(x=rcoeff, lw=4, color='r', alpha=0.5)
    ax.set_yscale('log')
    if itype == 'wtmean':
        ax.set_ylim(10**-1.5, 10**3)
    else:
        ax.set_ylim(10**-0.5, 10**4)
    ax.set_xlim(0,rmax)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(14)
    ax.text(0.95, 0.9, gal, horizontalalignment='right', fontsize=16,
        verticalalignment='center', transform=ax.transAxes)
    return
