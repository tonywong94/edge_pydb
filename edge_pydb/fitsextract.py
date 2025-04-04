import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy.wcs import WCS
from edge_pydb.conversion import gc_polr
from .hexgrid import hex_sampler


def fitsextract(input, header=None, stride=[1,1,1], keepref=True, keepnan=True, 
                zselect=None, col_lbl='imdat', ra_gc=None, dec_gc=None,
                pa=0, inc=0, ortlabel='default', bunit=None, first=False,
                use_hexgrid=False, sidelen=2, starting_angle=0, precision=2, 
                header_hex=None, hexgrid_output=None):

    """
    Sample data from an image into an AstroPy table indexed by coordinates.
    Simple approach taking every nth pixel along each image axis.
    Pseudocubes are handled as separate images and are detected by a blank 
    value for CTYPE3 in the header.

    Parameters
    ----------
    input : str or :class:`~numpy.ndarray`
        The input image or cube to turn into a table. This can be:
            * The name of a FITS file
            * A numpy array (in which case header must be separately given)
    header : :class:`~astropy.io.fits.Header` object
        Header corresponding to the input array.  Must be provided if the
        input is a numpy array.
    stride : tuple of ints, optional
        Step size to select pixels along each axis.  Axes are ordered using
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
    tab : :class:`~astropy.Table`
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

    # Create the coordinate columns.  First case is for PPV cubes.
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
        tab['ix'].description = '0-based pixel index in x direction'
        tab['iy'].description = '0-based pixel index in y direction'
        tab['iz'].description = '0-based pixel index in z direction'
        # Get the pixel coordinates as tuples
        wcsin = (np.array([tab['ix'],tab['iy'],tab['iz']])).T
    # Next case is for 2D images or pseudo-cubes
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
        tab['ix'].description = '0-based pixel index in x direction'
        tab['iy'].description = '0-based pixel index in y direction'
        # Get the pixel coordinates as tuples
        wcsin = (np.array([tab['ix'],tab['iy']])).T
    # Reduce WCS to 2-D for pseudocubes
    wfix = w.sub(naxis)
    
    # Calculate the coordinate columns
    if first:
        wcsout = wfix.wcs_pix2world(wcsin,0)
        col_ra = Column(wcsout.T[0], name='ra_abs',  dtype='f4', 
                        unit='deg', format='.6f', 
                        description='sample ra coord')
        col_dc = Column(wcsout.T[1], name='dec_abs', dtype='f4', 
                        unit='deg', format='.6f', 
                        description='sample dec coord')
        col_raoff = Column(wcsout.T[0]-w.wcs.crval[0], name='ra_off',  dtype='f4', 
                        unit='deg', format='.6f', 
                        description='ra offset from ref pixel')
        col_dcoff = Column(wcsout.T[1]-w.wcs.crval[1], name='dec_off', dtype='f4', 
                        unit='deg', format='.6f', 
                        description='dec offset from ref pixel')
        if ra_gc is None or ma.is_masked(ra_gc):
            ra_gc = w.wcs.crval[0]
        if dec_gc is None or ma.is_masked(dec_gc):
            dec_gc = w.wcs.crval[1]
        if not (inc>0):
            inc = 0.
        if ma.is_masked(pa) or ~np.isfinite(pa):
            pa = 0.
        r, theta = gc_polr(wcsout.T[0], wcsout.T[1], ra_gc, dec_gc, pa, inc)
        col_r = Column(r, name='rad_arc', dtype='f4', unit='arcsec', format='.3f',
            description='radius based on {}'.format(ortlabel))
        col_th = Column(theta, name='azi_ang', dtype='f4', unit='deg', format='.3f',
            description='azang based on {}'.format(ortlabel))
        tab.add_columns([col_ra,col_dc,col_raoff,col_dcoff,col_r,col_th])
        if iscube and not pseudo:
            col_vel = Column(wcsout.T[2]/1000., name='vel', dtype='f4', 
                unit='km/s', description='velocity in LSRK frame using radio def')
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
            try:
                desc = hdr['DESC_'+str(sel)].strip()
            except:
                desc = ''
            col_data = Column(np.ravel(data[sel],order='F'), 
                              name=col_lbl[iz], dtype='f4', 
                              description=desc, unit=bunit[iz])
            tab.add_column(col_data)
    else:
        if isinstance(bunit, list):
            bunit = bunit[0]
        if isinstance(col_lbl, list):
            col_lbl = col_lbl[0]
        col_data = Column(np.ravel(data,order='F'), name=col_lbl, dtype='f4', unit=bunit)
        tab.add_column(col_data)

    idx = ['ix', 'iy', 'iz']
    rem = [0, 0, 0]
    select = [[],[],[]]
    if not use_hexgrid:
        # Use stride to select the desired rows from the full table
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
    else:   # for hexgrid
        if iscube and not pseudo:
            iz_data = []
            tab_length = 0
            zlist = np.where(np.unique(tab['iz']) % stride[2] == rem[2])[0]
            for iz in zlist:
                if len(tab[tab['iz'] == iz]) == 0:
                    continue
                sample = hex_sampler(tab[tab['iz'] == iz], sidelen, keepref, 
                                     wfix.wcs.crpix[:2]-1, w.wcs.crval[0], 
                                     w.wcs.crval[1], ra_gc, dec_gc, pa, inc,
                                     starting_angle, precision, hexgrid_output)
                iz_data.append(sample)
                sample['iz'] = iz
                tab_length += len(sample)
            tab = tab[:tab_length]
            init = 0
            for tabs in iz_data:
                tab[init:(init+len(tabs))] = tabs
                init += len(tabs)
        else:
            sample = hex_sampler(tab, sidelen, keepref, wfix.wcs.crpix[:2]-1, 
                                 w.wcs.crval[0], w.wcs.crval[1], ra_gc, dec_gc, 
                                 pa, inc, starting_angle, precision, hexgrid_output)
            # Copy over the column descriptions
            if len(tab.columns) == len(sample.columns):
                for i in range(len(sample.columns)):
                    sample.columns[i].description = tab.columns[i].description
            tab = sample

    # Remove NaN rows if desired
    if not keepnan:
        if not pseudo:
            newtab = tab[~np.isnan(tab[col_lbl])]
            tab = newtab
        else:
            df = tab.to_pandas()
            df.dropna(how='all', subset=col_lbl)
            tab = Table.from_pandas(df)
    return tab

# -----------------------------------------------------
# getlabels: Interpret the FITS output from Pipe3D
# -----------------------------------------------------

def getlabels(product, p3dstruct='califa'):
    if product == 'ELINES':
        if p3dstruct == 'califa':
            nz = 20
            has_errors = True
            fluxlike = list(range(11))[2:]
        elif p3dstruct in ['manga', 'ecalifa']:
            nz = 11
            has_errors = False
            fluxlike = list(range(nz))[2:]
        zsel = range(nz)
        bright = ['[OII]3727', '[OIII]5007', '[OIII]4959',
                  'Hbeta'    , 'Halpha'    , '[NII]6583', 
                  '[NII]6548', '[SII]6731' , '[SII]6717']
        if has_errors:
            elbl = ['e_'+txt for txt in bright]
        else:
            elbl = []
        lbl  = ['Havel', 'Vdisp'] + bright + elbl
        units = lbl.copy()
        units[0:2] = ['km/s', 'Angstrom']
        units[2:]  = ['10^-16 erg cm^-2 s^-1']*(nz-2)
    elif product == 'flux_elines':
        # We select the bright lines that are also in ELINES, plus [OI]6300
        has_errors = True
        flux  = [0, 26, 27, 28, 41, 45, 46, 47, 49, 50]
        if p3dstruct == 'califa':
            nz = 408
        elif p3dstruct == 'ecalifa':
            nz = 432
        elif p3dstruct == 'manga':
            nz = 456
        elif p3dstruct == 'amusing':
            nz = 240
            flux  = [1, 2, 3, 19, 20, 21, 22, 24, 25]
        nline = len(flux)
        nfelines = nz // 8
        vel   = list(np.array(flux)+nfelines)
        disp  = list(np.array(flux)+nfelines*2)
        ew    = list(np.array(flux)+nfelines*3)
        eflux = list(np.array(flux)+nfelines*4)
        evel  = list(np.array(flux)+nfelines*5)
        edisp = list(np.array(flux)+nfelines*6)
        eew   = list(np.array(flux)+nfelines*7)
        zsel  = flux + vel + disp + ew + eflux + evel + edisp + eew
        flbl  = ['flux_[OII]3727', 'flux_[OIII]5007', 'flux_[OIII]4959', 
                 'flux_Hbeta',     'flux_[OI]6300',   'flux_Halpha',   
                 'flux_[NII]6583', 'flux_[NII]6548', 'flux_[SII]6717', 
                 'flux_[SII]6731']
        if p3dstruct == 'amusing':
            flbl = flbl[1:]
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
        fluxlike = flux.copy()
    elif product == 'indices':
        nz = 18
        has_errors = True
        if p3dstruct == 'ecalifa':
            idx = [22, 9, 13, 14, 15, 0, 33, 23, 34]
            e_idx  = list(np.array(idx)+35)
            zsel = idx + e_idx
        else:
            zsel = range(nz)
        albl = ['Hdel_idx',   'Hbet_idx',   'Mgb_idx', 
                'Fe5270_idx', 'Fe5335_idx', 'D4000_idx', 
                'Hdmod_idx',  'Hgam_idx',   'SN_idx']
        aunits = ['Angstrom']*8 + ['10^-16 erg cm^-2 s^-1']
        elbl = ['e_'+txt for txt in albl]
        lbl = albl + elbl
        units = aunits + aunits
        fluxlike = []
    elif product == 'SFH':
        if p3dstruct == 'califa':
            nz = 398    # 2*(39*4 + 39 + 4)
            has_errors = True
            # Note ages are in string rather than float order!
            ages = ['0.0010', '0.0040', '0.0030', '0.0056', '0.0089', '0.0126', '0.0141', 
                    '0.0178', '0.0199', '0.0100', '0.0251', '0.0316', '0.0398', '0.0562', 
                    '0.0631', '0.0630', '0.0708', '0.1122', '0.1259', '0.1585', '0.1995', 
                    '0.1000', '0.2818', '0.3548', '0.5012', '0.7079', '0.8913','10.0000', 
                    '1.1220','12.5893', '1.2589','14.1254', '1.4125', '1.9953', '2.5119', 
                    '3.5481', '4.4668', '6.3096', '7.9433']
            mets = ['0.0037', '0.0076', '0.0190', '0.0315']
        elif p3dstruct == 'amusing':
            nz = 199    # (39*4 + 39 + 4)
            has_errors = False
            # Note ages are in string rather than float order!
            ages = ['0.0010', '0.0040', '0.0030', '0.0056', '0.0089', '0.0126', '0.0141', 
                    '0.0178', '0.0199', '0.0100', '0.0251', '0.0316', '0.0398', '0.0562', 
                    '0.0631', '0.0630', '0.0708', '0.1122', '0.1259', '0.1585', '0.1995', 
                    '0.1000', '0.2818', '0.3548', '0.5012', '0.7079', '0.8913','10.0000', 
                    '1.1220','12.5893', '1.2589','14.1254', '1.4125', '1.9953', '2.5119', 
                    '3.5481', '4.4668', '6.3096', '7.9433']
            mets = ['0.0037', '0.0076', '0.0190', '0.0315']
        elif p3dstruct in ['manga', 'ecalifa']:
            nz = 319    # 39*7 + 39 + 7
            has_errors = False
            ages = ['0.0010', '0.0023', '0.0038', '0.0057', '0.0080', '0.0115', '0.0150', 
                    '0.0200', '0.0260', '0.0330', '0.0425', '0.0535', '0.0700', '0.0900', 
                    '0.1100', '0.1400', '0.1800', '0.2250', '0.2750', '0.3500', '0.4500', 
                    '0.5500', '0.6500', '0.8500', '1.1000', '1.3000', '1.6000', '2.0000', 
                    '2.5000', '3.0000', '3.7500', '4.5000', '5.2500', '6.2500', '7.5000', 
                    '8.5000','10.2500','12.0000','13.5000']
            mets = ['0.0001', '0.0005', '0.0020', '0.0080', '0.0170', '0.0300', '0.0400']
        zsel = range(nz)
        albl = []
        elbl = []
        for age in ages:
            for met in mets:
                albl.append('lumfrac_age_'+age+'_met_'+met)
                if has_errors:
                    elbl.append('e_lumfrac_age_'+age+'_met_'+met)
        ages_sort = sorted(ages, key=lambda x: float(x))
        for age in ages_sort:
            albl.append('lumfrac_age_'+age)
            if has_errors:
                elbl.append('e_lumfrac_age_'+age)
        for met in mets:
            albl.append('lumfrac_met_'+met)
            if has_errors:
                elbl.append('e_lumfrac_met_'+met)
        lbl  = albl + elbl
        units = ['fraction']*len(lbl)
        fluxlike = []
    elif product == 'SSP':
        if p3dstruct in ['califa', 'amusing']:
            nz = 20
            has_errors = False
        elif p3dstruct in ['manga', 'ecalifa']:
            nz = 21
            has_errors = True
        zsel = range(nz)
        lbl = ['Vcont_ssp',   'cont_segm','cont_dezon','medflx_ssp', 
               'e_medflx_ssp','age_lwt',  'age_mwt',   'e_age_lwt', 
               'ZH_lwt',      'ZH_mwt',   'e_ZH_lwt',  'Av_ssp', 
               'e_Av_ssp',    'vel_ssp',  'e_vel_ssp', 'vdisp_ssp', 
               'e_vdisp_ssp', 'ML_ssp',   'mass_ssp',  'mass_Avcor_ssp']
        units = ['10^-16 erg cm^-2 s^-1 AA^-1', None, None, 
                 '10^-16 erg cm^-2 s^-1 AA^-1', '10^-16 erg cm^-2 s^-1 AA^-1', 
                 'dex(yr)', 'dex(yr)', 'fraction', 'dex', 'dex', 'fraction', 
                 'mag', 'mag', 'km/s', 'km/s', 'km/s', 'km/s', 
                 'solMass/solLum', 'dex(solMass/pixel^2)', 'dex(solMass/pixel^2)']
        fluxlike = [0]
        if has_errors:
            lbl += ['e_mass_ssp']
            units += ['dex(solMass/pixel^2)']
    return zsel, lbl, units, len(zsel), has_errors, fluxlike

