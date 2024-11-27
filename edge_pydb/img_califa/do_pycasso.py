#!/usr/bin/env python

# Combine the CALIFA PyCASSO data (GMe) into binary tables.

from datetime import datetime
import glob
import os
import numpy as np
from astropy import units as u
from astropy.table import Table, Column, join, vstack
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from reproject import reproject_interp
from edge_pydb import EdgeTable
from edge_pydb.conversion import stmass_pc2, sfr_ha, ZOH_M13, bpt_type, get_AHa
from edge_pydb.fitsextract import fitsextract, getlabels
np.seterr(divide='ignore', invalid='ignore')

def do_pycasso(outfile='NGC4047.pipe3d.hdf5', filelist=None, 
             comomdir=None, colabel='co.smo7',
             ext='', ortpar='edge_leda.csv', distpar='edge_califa.csv',
             distcol='caDistP3d', hexgrid=False, allpix=False, debug=False, 
             discard_cdmatrix=False, append=True, overwrite=True):
    """
    Extract Pipe3D products into an HDF5 database.  This script assumes
    there are 5 Pipe3D output files per galaxy.

    Parameters
    ----------
    outfile : str
        Name of the output filename.  Appended to if it exists.
    filelist : list of str
        List of input files
    comomdir : str
        Path to the directory where CO moments FITS files reside.  Default is
        process PyCASSO data only without regridding to match CO.
    colabel : str
        Identifier for reference line in the CO FITS filenames
    ext : str
        Suffix to add to column names, e.g. '_sm'
    ortpar : filename
        Name of the EdgeTable which has LEDA orientation parameters for the sample
    hexgrid : boolean
        True to sample on a hexagonal grid (experimental)
    allpix : boolean
        True to dump every pixel, otherwise every 3rd pixel in x and y is used.
    debug : boolean
        True to generate some additional output
    append : boolean
        True to append to an existing file.  This is the default (write to the
        file created by do_comom.py).
    overwrite : boolean
        True to overwrite existing tables.  This is the default (replace same table
        but do not delete other tables in the file).
    """
    if allpix:
        stride = [1,1,1]
    else:
        stride = [3,3,1]

    # cuts for when to apply BD correction
    hacut = 0.06    # 1e-16 erg / (cm2 s)
    hbcut = 0.04    # 1e-16 erg / (cm2 s)
    ahalo = 0       # mag
    ahahi = 6       # mag

    # FITS keywords important for astrometry
    wcskeys = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 
               'CDELT1', 'CDELT2']
    cdkeys = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CD1_3', 'CD2_3',
                'CD3_1', 'CD3_2', 'CD3_3']
    dimkeys = ['NAXIS1', 'NAXIS2']

    # Get the orientation parameters from LEDA
    orttbl = EdgeTable(ortpar)
    orttbl.add_index('ID') 

    # Read the FITS data
    # The columns to save are defined in fitsextract.py
    prodtype = ['sigma_star', 'sigma_star_ini', 'L_5635', 'log_age_flux',
                'log_age_mass', 'log_Z_flux', 'log_Z_mass', 'sigma_sfr',
                'x_young', 'tau_V', 'v_star', 'v_star_disp',
                'adev', 'nlambda_clip', 'chi2', 'badpix',
                'zones', 'sn_pix', 'sn_zone']
    units = ['solMass/pc^2', 'solMass/pc^2', 'solLum pc^-2 AA^-1', 'dex(yr)',
                'dex(yr)', 'dex', 'dex', 'solMass pc^-2 Gyr^-1',
                '', '', 'km/s', 'km/s',
                'pct', 'pct', '', '',
                '', '', '']

    if len(filelist) == 0:
        raise RuntimeError('Error: filelist is empty!')

    tlist = []

    for fname in filelist:
        if not os.path.exists(fname):
            print('####### Cannot find',fname)
            continue
        else:
            califa_id = int(os.path.basename(fname)[1:5])
            gal = orttbl.loc[califa_id]['Name']

        for i_prod, prod in enumerate(prodtype):

            hdu = fits.open(fname)[i_prod]
            print('\nWorking on galaxy {} product {}'.format(gal,prod))
            cahd = hdu.header.copy()

            # Read in CO template
            if comomdir is not None:
                cofile = os.path.join(comomdir,gal+'.'+colabel+'_dil.snrpk.fits.gz')
                if not os.path.exists(cofile):
                    print('####### Cannot find',cofile)
                    continue
                cohd = fits.getheader(cofile)
                # Copy the CALIFA header and replace wcskeys with CO values
                for key in dimkeys+wcskeys:
                    if key in cohd.keys():
                        cahd[key] = cohd[key]
                # Need to discard CD matrix which would override the new wcskeys
                if 'CDELT1' in cohd.keys() and 'CDELT2' in cohd.keys():
                    for key in cdkeys:
                        if key in cahd.keys():
                            del cahd[key]
                # Optionally discard CD matrix in CALIFA files and fall back on CDELTs
                if discard_cdmatrix:
                    for key in cdkeys:
                        if key in hdu.header.keys():
                            del hdu.header[key]
                if debug:
                    print('\nINPUT',WCS(hdu.header))
                    print('\nCO data',WCS(cohd))
                    print('\nOUTPUT',WCS(cahd))
                newim = reproject_interp(hdu, cahd, order=0, return_footprint=False)
                if debug:
                    fits.writeto(fname.replace('.fits','.rg.fits'), newim, cahd, 
                                 overwrite=True)
            else:
                newim = hdu.data
                if hdu.header['bitpix'] < 0:
                    newim[newim == 0.] = np.nan

            if i_prod == 0:
                print("RA, DEC, PA, INC:",orttbl.loc[califa_id]['ledaRA'],
                      orttbl.loc[califa_id]['ledaDE'], orttbl.loc[califa_id]['ledaPA'],
                      orttbl.loc[califa_id]['ledaAxIncl'])
                galtab = fitsextract(newim, header=cahd, keepnan=True, stride=stride, 
                               bunit=units[i_prod], col_lbl=prod, 
                               ra_gc=15*orttbl.loc[califa_id]['ledaRA'],
                               dec_gc=orttbl.loc[califa_id]['ledaDE'], 
                               pa=orttbl.loc[califa_id]['ledaPA'],
                               inc=orttbl.loc[califa_id]['ledaAxIncl'], 
                               ortlabel='LEDA', first=True, use_hexgrid=hexgrid)
                gname = Column([np.string_(gal)]*len(galtab), name='Name', 
                               description='Galaxy Name')
                galtab.add_column(gname, index=0)
            else:
                addtb = fitsextract(newim, header=cahd, keepnan=True, stride=stride, 
                                bunit=units[i_prod], col_lbl=prod, use_hexgrid=hexgrid)
                jointb = join(galtab, addtb, keys=['ix','iy'])
                galtab = jointb

        tlist.append(galtab)

    if len(tlist) > 0:
        t_merge = vstack(tlist)
    t_merge['sigma_star'].description = 'Stellar mass surface density'
    t_merge['sigma_star_ini'].description = 'Initial stellar mass surface density'
    t_merge['L_5635'].description = 'Luminosity surface density in normalization window'
    t_merge['log_age_flux'].description = 'Mean log of stellar age lum weighted'
    t_merge['log_age_mass'].description = 'Mean log of stellar age mass weighted'
    t_merge['log_Z_flux'].description = 'Mean log of stellar met lum weighted'
    t_merge['log_Z_mass'].description = 'Mean log of stellar met mass weighted'
    t_merge['sigma_sfr'].description = 'Star formation surface density last 32 Myr'
    t_merge['x_young'].description = 'Lum fraction of stellar pops under 32 Myr'
    t_merge['tau_V'].description = 'Attenuation coefficient for dust screen model'
    t_merge['v_star'].description = 'Line of sight stellar velocity'
    t_merge['v_star_disp'].description = 'Line of sight stellar velocity dispersion'
    t_merge['adev'].description = 'Mean model deviation'
    t_merge['nlambda_clip'].description = '% of wavelengths clipped by fitting algorithm'
    t_merge['chi2'].description = 'Fit statistic'
    t_merge['badpix'].description = 'Masked pixels'
    t_merge['zones'].description = 'Voronoi segmentation zones'
    t_merge['sn_pix'].description = 'SNR in individual pixels'
    t_merge['sn_zone'].description = 'SNR in zones'
    t_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    if debug:
        print(t_merge.colnames)
        print('There are',len(t_merge),'rows in merged table')

    t_merge.write(outfile, path='starlight'+ext, overwrite=overwrite, 
            append=append, serialize_meta=True, compression=True)

    return


if __name__ == "__main__":
    # All DR3 galaxies, PyCASSO only
    filelist = sorted(glob.glob('fits_pycasso/K*_gsd6e.fits.gz'))
    do_pycasso(filelist=filelist, outfile='dr3_allpix.pycasso.hdf5', 
              append=False, allpix=True)
