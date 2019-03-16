import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy.wcs import WCS

def fitsextract(filename, stride=[1,1,1], keepref=True, keepnan=True):
    hdu   = fits.open(filename)[0]
    data  = hdu.data
    bunit = hdu.header['BUNIT']
    w     = WCS(hdu)
    print('RA ref is',w.wcs.crval[0])
    print('DEC ref is',w.wcs.crval[1])
    ndim  = len(data.shape)
    iscube = (ndim > 2 and data.shape[ndim-3] > 1)
    if iscube:
        print('This is a data cube')
        data = np.squeeze(data)
        if len(data.shape) > 3:
            raise ('Data cannot be squeezed to three dimensions')
        naxis = 3
        ix,iy,iz = np.meshgrid(np.arange(data.shape[2]),np.arange(data.shape[1]),
                             np.arange(data.shape[0]),indexing='ij')
        tab = Table([np.ravel(ix),np.ravel(iy),np.ravel(iz)],
                  names=('ix','iy','iz'),dtype=('i4','i4','i4'))
        # Get the pixel coordinates as tuples
        wcsin = (np.array([tab['ix'],tab['iy'],tab['iz']])).T
    else:
        print('This is an image')
        data = np.squeeze(data)
        naxis = 2
        ix,iy = np.meshgrid(np.arange(data.shape[1]),np.arange(data.shape[0]),
                            indexing='ij')
        tab = Table([np.ravel(ix),np.ravel(iy)],
                    names=('ix','iy'),dtype=('i4','i4'))
        # Get the pixel coordinates as tuples
        wcsin = (np.array([tab['ix'],tab['iy']])).T
    wfix = w.sub(naxis)
    wcsout = wfix.wcs_pix2world(wcsin,0)
    col_ra = Column(wcsout.T[0]-w.wcs.crval[0], name='ra_off',  dtype='f4', unit='deg')
    col_dc = Column(wcsout.T[1]-w.wcs.crval[1], name='dec_off', dtype='f4', unit='deg')
    tab.add_columns([col_ra,col_dc])
    if iscube:
        col_vel = Column(wcsout.T[2]/1000., name='vel', dtype='f4', unit='km/s')
        tab.add_column(col_vel)
    # Use order = 'F' because the velocity axis is first
    col_data = Column(np.ravel(data,order='F'), name='imgdata', dtype='f4', unit=bunit)
    tab.add_column(col_data)
    # Select the desired rows from the full table
    idx = ['ix', 'iy', 'iz']
    rem = [0, 0, 0]
    select = [[],[],[]]
    if keepref:
        for i in range(naxis):
            crpix = wfix.wcs.crpix[i]
            if crpix < 1 or crpix > data.shape[naxis-i-1] or not crpix.is_integer():
                print('Cannot use keepref on axis {}: crpix={}'.format(i+1,crpix))
                continue
            else:
                print('Axis {}: crpix={}'.format(i+1,crpix))
                rem[i] = int(crpix-1) % stride[i]
        print('Remainder: ',rem)
    for i in range(naxis):
        select[i]=np.where(tab[idx[i]] % stride[i] == rem[i])[0]
    xy = np.intersect1d(select[0], select[1])
    if iscube:
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
        newtab = tab[~np.isnan(tab['imgdata'])]
        tab = newtab
    return tab

