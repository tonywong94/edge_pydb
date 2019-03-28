import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy.wcs import WCS

def fitsextract(input, stride=[1,1,1], keepref=True, keepnan=True, header=None, 
                lbl='index',bunit=None, col_name=None, zselect=None, suffix=''):
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
    pseudo = (hdr['ctype3'] == '')
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
    wcsout = wfix.wcs_pix2world(wcsin,0)
    col_ra = Column(wcsout.T[0]-w.wcs.crval[0], name='ra_off',  dtype='f4', unit='deg')
    col_dc = Column(wcsout.T[1]-w.wcs.crval[1], name='dec_off', dtype='f4', unit='deg')
    tab.add_columns([col_ra,col_dc])
    if iscube and not pseudo:
        col_vel = Column(wcsout.T[2]/1000., name='vel', dtype='f4', unit='km/s')
        tab.add_column(col_vel)
    # Use order = 'F' because the velocity axis is first
    if pseudo:
        zsel = range(nz) if zselect is None else zselect
        for i in zsel:
            col_data = Column(np.ravel(data[i],order='F'), 
                              name=hdr[lbl+str(i)]+suffix, dtype='f4', unit=bunit)
            tab.add_column(col_data)
    else:
        cname = 'imgdata' if col_name is None else col_name
        col_data = Column(np.ravel(data,order='F'), name=cname, dtype='f4', unit=bunit)
        tab.add_column(col_data)
    # Select the desired rows from the full table
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
    return tab

# getlabels: Interpret the FITS output from Pipe3D

def getlabels(product):
    if product == 'ELINES':
        nz = 20
        zsel = range(nz)
        lbl = ['Havel', 'Vdisp', 
               'OII_3727', 'OIII_5007', 'OIII_4959', 'Hbeta', 
               'Halpha', 'NII_6583', 'NII_6548', 'SII_6731', 'SII_6717',
               'e_OII_3727', 'e_OIII_5007', 'e_OIII_4959', 'e_Hbeta', 
               'e_Halpha', 'e_NII_6583', 'e_NII_6548', 'e_SII_6731', 'e_SII_6717']
        units = lbl.copy()
        units[0:1] = ['km/s', 'Angstrom']
        units[2:] = ['10^-16 erg cm^-2 s^-1']*(nz-2)
    elif product == 'flux_elines':
        # We only select the bright lines that were in ELINES
        nz = 408
        zsel = []
        lbl = []
        units = lbl.copy()
    elif product == 'indices':
        nz = 18
        zsel = range(nz)
        lbl = []
        units = ['Angstrom']*nz
    elif product == 'sfh':
        # We only select the age bins
        nz = 398
        z_age = list(range(156,195))
        n_age = len(z_age)
        z_err = list(range(355,382))
        n_err = len(z_err)
        zsel = z_age + z_err
        ages = ['0.0010', '0.0030', '0.0040', '0.0056', '0.0089', '0.0100',
                '0.0126', '0.0141', '0.0178', '0.0199', '0.0251', '0.0316',
                '0.0398', '0.0562', '0.0630', '0.0631', '0.0708', '0.1000',
                '0.1122', '0.1259', '0.1585', '0.1995', '0.2818', '0.3548',
                '0.5012', '0.7079', '0.8913', '1.1220', '1.2589', '1.4125',
                '1.9953', '2.5119', '3.5481', '4.4668', '6.3096', '7.9433', 
                '10.000', '12.5893', '14.1254']
        albl = ['lumfrac_age_'+age for age in ages]
        elbl = ['e_lumfrac_age_'+age for age in ages]
        lbl = albl + elbl
        units = ['fraction']*len(lbl)
    elif product == 'ssp':
        nz = 20
        zsel = range(nz)
        lbl = ['Vcont_ssp', 'cont_segm', 'cont_dezon', 'medflx_ssp', 
               'e_medflx_ssp', 'age_lwt', 'age_mwt', 'e_age_lwt', 
               'ZH_lwt', 'ZH_mwt', 'e_ZH_lwt', 'Av_ssp', 
               'e_Av_ssp', 'vel_ssp', 'e_vel_ssp', 'vdisp_ssp', 
               'e_vdisp_ssp', 'ML_ssp', 'mass_ssp', 'mass_Avcor_ssp']
        units = ['10^-16 erg cm^-2 s^-1 AA^-1', 'none', 'none', '10^-16 erg cm^-2 s^-1 AA^-1', 
                 '10^-16 erg cm^-2 s^-1', 'dex(yr)', 'dex(yr)', 'fraction', 
                 'dex', 'dex', 'fraction', 'mag', 
                 'mag', 'km/s', 'km/s', 'km/s', 
                 'km/s', 'solMass/solLum', 'dex(solMass/arcsec^2)', 'dex(solMass/arcsec^2)']
    return zsel, lbl, units, len(zsel)
