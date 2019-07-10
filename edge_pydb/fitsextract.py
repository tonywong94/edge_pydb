import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy.wcs import WCS
from edge_pydb.conversion import gc_polr

def fitsextract(input, header=None, stride=[1,1,1], keepref=True, keepnan=True, 
                zselect=None, col_lbl='imdat', ra_gc=None, dec_gc=None,
                pa=0, inc=0, ortlabel='default', bunit=None, first=False):

    """
    Sample data from an image into an AstroPy table indexed by coordinates.
    Simple approach taking every nth pixel along each image axis.
    Pseudocubes are handled as separate images and are detected by a blank 
    value for CTYPE3 in the header.

    Parameters
    ----------
    input : str or `~numpy.ndarray`
        The input data to reproject. This can be:
            * The name of a FITS file
            * A numpy array (in which case header must be separately given)
    header : `~astropy.io.fits.Header` object
        Header corresponding to the input array.  Must be provided if the
        input is a numpy array.
    stride : tuple of ints, optional
        step size to select pixels along each axis.  Axes are ordered using
        the FITS convention, not numpy convention (i.e. velaxis last).
        Default is [1,1,1] to keep all pixels.
        Note: stride in z is ignored for pseudocubes.
    keepref : bool, optional
        If dropping pixels, try to ensure that the reference pixel is kept.
        Default is True (keep the reference pixel).
    keepnan : bool, optional
        If False, the output table drops rows which are all-NaN.
        Default is True (keep the NaNs).
    zselect : list of ints, optional
        Indices of planes in cube/pseudocube to be selected.
        Default is to keep all planes.
    col_lbl : string or list of strings, optional
        Column label for the data values, can be list corresponding to each 
        plane to be selected, for CALIFA pseudocubes.
        Default is "imdat"+possible integer.
    ra_gc : float, optional
        Right ascension of galaxy center.  Used to determine
        polar coordinates of each sample in the plane of the galaxy.
        Default is reference RA value of the image.
    dec_gc : float, optional
        Declination of galaxy center.  Used to determine
        polar coordinates of each sample in the plane of the galaxy.
        Default is reference DEC value of the image.
    pa : float, optional
        Position angle of the galaxy disk in degrees E of N.  Used to determine
        polar coordinates of each sample in the plane of the galaxy.
        Default is 0 (due north).
    inc : float, optional
        Inclination of the galaxy disk in degrees from face-on.  Used to
        determine polar coordinates of each sample in the plane of the galaxy.
        Default is 0 (face-on).
    ortlabel : string, optional
        String label describing the origin of the orientation parameters.
    bunit : string or list of strings, optional
        Astropy units for data values, can be list corresponding to each plane
        to be selected, for CALIFA pseudocubes.
        Default is obtained from BUNIT in the header.
    first : bool, optional
        True to write the coordinate columns, which only need to be done once
        per table.  When combining multiple FITS images into a table, fitsextract
        should be called initially with first=True and then subsequently with
        first=False.
        Default is False.

    Returns
    -------
    tab : `~astropy.Table`
        The selected pixels as a 1-D table.
    """

    # Read the inputs
    if isinstance(input, np.ndarray):
        if header is None:
            raise TypeError("Header must be given if input is not FITS file")
        else:
            hdr = header
            data = input
    else:
        hdu   = fits.open(input)[0]
        data  = hdu.data
        hdr   = hdu.header if header is None else header
    bunit = hdr['BUNIT'] if bunit is None else bunit
    w     = WCS(hdr)
    print('RA ref is',w.wcs.crval[0])
    print('DEC ref is',w.wcs.crval[1])
    ndim  = len(data.shape)
    iscube = (ndim > 2 and data.shape[ndim-3] > 1)
    if 'CTYPE3' in hdr.keys():
    	pseudo = (hdr['CTYPE3'] == '')
    else:
    	pseudo = False

    # Create the coordinate columns
    if iscube and not pseudo:
        print('This is a data cube of shape', data.shape)
        data = np.squeeze(data)
        if len(data.shape) > 3:
            raise ('Data cannot be squeezed to three dimensions')
        naxis = 3
        nz,ny,nx = data.shape
        ix,iy,iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                               indexing='ij')
        tab = Table([np.ravel(ix), np.ravel(iy), np.ravel(iz)],
                  names=('ix','iy','iz'),dtype=('i4','i4','i4'))
        # Get the pixel coordinates as tuples
        wcsin = (np.array([tab['ix'],tab['iy'],tab['iz']])).T
    else:
        print('This is an image of shape', data.shape)
        data = np.squeeze(data)
        naxis = 2
        if pseudo:
            nz,ny,nx = data.shape
        else:
            ny,nx = data.shape
        ix,iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        tab = Table([np.ravel(ix), np.ravel(iy)],
                    names=('ix','iy'),dtype=('i4','i4'))
        # Get the pixel coordinates as tuples
        wcsin = (np.array([tab['ix'],tab['iy']])).T
    wfix = w.sub(naxis)
    if first:
        wcsout = wfix.wcs_pix2world(wcsin,0)
        col_ra = Column(wcsout.T[0]-w.wcs.crval[0], name='ra_off',  dtype='f4', 
                        unit='deg', format='.6f')
        col_dc = Column(wcsout.T[1]-w.wcs.crval[1], name='dec_off', dtype='f4', 
                        unit='deg', format='.6f')
        if ra_gc is None:
            ra_gc = w.wcs.crval[0]
        if dec_gc is None:
            dec_gc = w.wcs.crval[1]
        r, theta = gc_polr(wcsout.T[0], wcsout.T[1], ra_gc, dec_gc, pa, inc)
        col_r = Column(r, name='rad_arc', dtype='f4', unit='arcsec', format='.3f',
            description='radius based on {}'.format(ortlabel))
        col_th = Column(theta, name='azi_ang', dtype='f4', unit='deg', format='.3f',
            description='azang based on {}'.format(ortlabel))
        tab.add_columns([col_ra,col_dc,col_r,col_th])
        if iscube and not pseudo:
            col_vel = Column(wcsout.T[2]/1000., name='vel', dtype='f4', unit='km/s')
            tab.add_column(col_vel)

    # Flatten the cube into a 1-D table
    # Use order = 'F' because the velocity axis is first
    # For pseudocubes, each plane becomes a separate column
    if pseudo:
        zsel = range(nz) if zselect is None else zselect
        if not isinstance(bunit, list):
            bunit = [bunit]*len(zsel)
        if not isinstance(col_lbl, list):
            col_lbl = [col_lbl+str(i) for i in range(len(zsel))]
        for iz, sel in enumerate(zsel):
            col_data = Column(np.ravel(data[sel],order='F'), 
                              name=col_lbl[iz], dtype='f4', unit=bunit[iz])
            tab.add_column(col_data)
    else:
        if isinstance(bunit, list):
            bunit = bunit[0]
        if isinstance(col_lbl, list):
            col_lbl = col_lbl[0]
        col_data = Column(np.ravel(data,order='F'), name=col_lbl, dtype='f4', unit=bunit)
        tab.add_column(col_data)

    # Use stride to select the desired rows from the full table
    idx = ['ix', 'iy', 'iz']
    rem = [0, 0, 0]
    select = [[],[],[]]
    if keepref:
        for i in range(naxis):
            crpix = wfix.wcs.crpix[i]
            if crpix < 1 or crpix > hdr['naxis'+str(i+1)] or not crpix.is_integer():
                print('Cannot use keepref on axis {}: crpix={}'.format(i+1,crpix))
                continue
            else:
                print('Axis {}: crpix={}'.format(i+1,crpix))
                rem[i] = int(crpix-1) % stride[i]
        print('Remainder: ',rem)
    for i in range(naxis):
        select[i]=np.where(tab[idx[i]] % stride[i] == rem[i])[0]
    xy = np.intersect1d(select[0], select[1])
    if iscube and not pseudo:
        xyz = np.intersect1d(xy, select[2])
        if len(xyz) < len(tab):
            newtab = tab[xyz]
            tab = newtab
    else:
        if len(xy) < len(tab):
            newtab = tab[xy]
            tab = newtab

    # Remove NaN rows if desired
    if not keepnan:
        if not pseudo:
            newtab = tab[~np.isnan(tab[cname])]
            tab = newtab
        else:
            df = tab.to_pandas()
            df.dropna(how='all', subset=col_lbl)
            tab = Table.from_pandas(df)
    return tab

# -----------------------------------------------------
# getlabels: Interpret the FITS output from Pipe3D
# -----------------------------------------------------

def getlabels(product):
    if product == 'ELINES':
        nz = 20
        zsel = range(nz)
        bright = ['[OII]3727','[OIII]5007','[OIII]4959','Hbeta',
                  'Halpha',   '[NII]6583', '[NII]6548', '[SII]6731', '[SII]6717']
        elbl = ['e_'+txt for txt in bright]
        lbl  = ['Havel', 'Vdisp'] + bright + elbl
        units = lbl.copy()
        units[0:2] = ['km/s', 'Angstrom']
        units[2:]  = ['10^-16 erg cm^-2 s^-1']*(nz-2)
    elif product == 'flux_elines':
        # We only select the bright lines that were in ELINES
        nz = 408
        flux  = [0, 26, 27, 28, 45, 46, 47, 49, 50]
        nline = len(flux)
        vel   = list(np.array(flux)+51)
        disp  = list(np.array(flux)+102)
        ew    = list(np.array(flux)+153)
        eflux = list(np.array(flux)+204)
        evel  = list(np.array(flux)+255)
        edisp = list(np.array(flux)+306)
        eew   = list(np.array(flux)+357)
        zsel  = flux + vel + disp + ew + eflux + evel + edisp + eew
        flbl  = ['flux_[OII]3727', 'flux_[OIII]5007', 'flux_[OIII]4959', 
                 'flux_Hbeta',     'flux_Halpha',     'flux_[NII]6583', 
                 'flux_[NII]6548', 'flux_[SII]6717',  'flux_[SII]6731']
        vlbl  = [ w.replace('flux', 'vel')    for w in flbl ]
        dlbl  = [ w.replace('flux', 'disp')   for w in flbl ]
        wlbl  = [ w.replace('flux', 'EW')     for w in flbl ]
        eflbl = [ w.replace('flux', 'e_flux') for w in flbl ]
        evlbl = [ w.replace('flux', 'e_vel')  for w in flbl ]
        edlbl = [ w.replace('flux', 'e_disp') for w in flbl ]
        ewlbl = [ w.replace('flux', 'e_EW')   for w in flbl ]
        lbl   = flbl + vlbl + dlbl + wlbl + eflbl + evlbl + edlbl + ewlbl
        units = lbl.copy()
        units[0*nline:1*nline] = ['10^-16 erg cm^-2 s^-1']*nline
        units[1*nline:2*nline] = ['km/s']*nline
        units[2*nline:3*nline] = ['Angstrom']*nline
        units[3*nline:4*nline] = ['Angstrom']*nline
        units[4*nline:5*nline] = ['10^-16 erg cm^-2 s^-1']*nline
        units[5*nline:6*nline] = ['km/s']*nline
        units[6*nline:7*nline] = ['Angstrom']*nline
        units[7*nline:8*nline] = ['Angstrom']*nline
    elif product == 'indices':
        nz = 18
        zsel = range(nz)
        albl = ['Hdel_idx',   'Hbet_idx',   'Mgb_idx', 
                'Fe5270_idx', 'Fe5335_idx', 'D4000_idx', 
                'Hdmod_idx',  'Hgam_idx',   'SN_idx']
        elbl = ['e_'+txt for txt in albl]
        lbl = albl + elbl
        units = ['Angstrom']*len(lbl)
    elif product == 'SFH':
        # We only select the age bins
        nz = 398
        z_age = list(range(156,195))
        n_age = len(z_age)
        z_err = list(range(355,394))
        n_err = len(z_err)
        zsel  = z_age + z_err
        ages  = ['0.0010', '0.0030', '0.0040', '0.0056', '0.0089', '0.0100',
                 '0.0126', '0.0141', '0.0178', '0.0199', '0.0251', '0.0316',
                 '0.0398', '0.0562', '0.0630', '0.0631', '0.0708', '0.1000',
                 '0.1122', '0.1259', '0.1585', '0.1995', '0.2818', '0.3548',
                 '0.5012', '0.7079', '0.8913', '1.1220', '1.2589', '1.4125',
                 '1.9953', '2.5119', '3.5481', '4.4668', '6.3096', '7.9433', 
                 '10.000', '12.5893', '14.1254']
        albl = ['lumfrac_age_'+age for age in ages]
        elbl = ['e_lumfrac_age_'+age for age in ages]
        lbl  = albl + elbl
        units = ['fraction']*len(lbl)
    elif product == 'SSP':
        nz = 20
        zsel = range(nz)
        lbl = ['Vcont_ssp',   'cont_segm','cont_dezon','medflx_ssp', 
               'e_medflx_ssp','age_lwt',  'age_mwt',   'e_age_lwt', 
               'ZH_lwt',      'ZH_mwt',   'e_ZH_lwt',  'Av_ssp', 
               'e_Av_ssp',    'vel_ssp',  'e_vel_ssp', 'vdisp_ssp', 
               'e_vdisp_ssp', 'ML_ssp',   'mass_ssp',  'mass_Avcor_ssp']
        units = ['10^-16 erg cm^-2 s^-1 AA^-1', 'none', 'none', '10^-16 erg cm^-2 s^-1 AA^-1', 
                 '10^-16 erg cm^-2 s^-1', 'dex(yr)', 'dex(yr)', 'fraction', 
                 'dex', 'dex', 'fraction', 'mag', 
                 'mag', 'km/s', 'km/s', 'km/s', 
                 'km/s', 'solMass/solLum', 'dex(solMass/arcsec^2)', 'dex(solMass/arcsec^2)']
    return zsel, lbl, units, len(zsel)
