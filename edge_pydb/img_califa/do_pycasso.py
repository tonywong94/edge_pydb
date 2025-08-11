#!/usr/bin/env python

# Combine the CALIFA PyCASSO data (GMe) into binary tables.

from datetime import datetime
import glob
import os
import numpy as np
from astropy import units as u
from astropy.table import Table, Column, join, vstack
from astropy.io import fits
from edge_pydb import EdgeTable
from edge_pydb.fitsextract import fitsextract, getlabels
np.seterr(divide='ignore', invalid='ignore')

def do_pycasso(outfile='dr3_allpix.pycasso.hdf5', filelist=None, 
               ext='', ortpar='edge_leda.csv', distpar='edge_califa.csv',
               distcol='caDistP3d', hexgrid=False, allpix=False, debug=False, 
               append=True, overwrite=True, regrid_to_p3d=False):
    """
    Extract Pipe3D products into an HDF5 database.  This script assumes
    there are 5 Pipe3D output files per galaxy.

    Parameters
    ----------
    outfile : str
        Name of the output filename.
    filelist : list of str
        List of input files
    ext : str
        Suffix to add to column names, e.g. '_sm'
    ortpar : filename
        Name of the EdgeTable which has LEDA orientation parameters for the sample
    distpar : filename
        Name of the EdgeTable which has distances for converting Sigma_*.
    distcol : str
        Name of the distance column in 'distpar' to use.  Default is 'caDistP3d'
        taken from 'DL' column in get_proc_elines_CALIFA.csv.
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

    # Get the IDs from CALIFA and orientation parameters from LEDA
    orttbl = EdgeTable(ortpar)
    orttbl.add_index('Name') 
    disttbl = EdgeTable(distpar)
    disttbl.add_index('ID') 

    # Read the FITS data
    # The columns to save are defined in fitsextract.py
    prodtype = ['sigma_star', 'sigma_star_ini', 'L_5635', 'log_age_flux',
                'log_age_mass', 'log_Z_flux', 'log_Z_mass', 'sigma_sfr',
                'x_young', 'tau_V', 'v_0', 'v_d',
                'adev', 'nl_clip', 'chi2', 'zones', 'sn', 'sn_zone']
    units = ['solMass/pc^2', 'solMass/pc^2', 'solLum pc^-2 AA^-1', 'dex(yr)',
                'dex(yr)', 'dex', 'dex', 'solMass pc^-2 Gyr^-1',
                '', '', 'km/s', 'km/s',
                'pct', 'pct', '', '', '', '']

    if len(filelist) == 0:
        raise RuntimeError('Error: filelist is empty!')

    tlist = []

    for fname in filelist:
        if not os.path.exists(fname):
            print('####### Cannot find',fname)
            continue
        else:
            califa_id = int(os.path.basename(fname)[1:5])
            gal = disttbl.loc[califa_id]['Name']

        hdul = fits.open(fname)
        msk = hdul['BADPIX'].data

        for i_prod, prod in enumerate(prodtype):

            #hdu = fits.open(fname)[i_prod]
            print('\nWorking on galaxy {} product {}'.format(gal,prod))

            cahd  = hdul[prod].header
            newim = hdul[prod].data
            if cahd['bitpix'] < 0:
                newim[msk>0] = np.nan
#                 newim[newim == 0.] = np.nan

#             if regrid_to_p3d:  # Not yet working
#                 cadat = reproject_interp((cadat,p3dhd), WCS(w_cahd), order=interp_order,
#                     shape_out=(cadat.shape[0],w_cahd['NAXIS2'],w_cahd['NAXIS1']),
#                     return_footprint=False)

            if i_prod == 0:
                print("RA, DEC, PA, INC:",orttbl.loc[gal]['ledaRA'],
                      orttbl.loc[gal]['ledaDE'], orttbl.loc[gal]['ledaPA'],
                      orttbl.loc[gal]['ledaAxIncl'])
                galtab = fitsextract(newim, header=cahd, keepnan=True, stride=stride, 
                               bunit=units[i_prod], col_lbl=prod, 
                               ra_gc=orttbl.loc[gal]['ledaRA'],
                               dec_gc=orttbl.loc[gal]['ledaDE'], 
                               pa=orttbl.loc[gal]['ledaPA'],
                               inc=orttbl.loc[gal]['ledaAxIncl'], 
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
    t_merge['x_young'].description = 'Luminosity fraction of stellar pops under 32 Myr'
    t_merge['tau_V'].description = 'Attenuation coefficient for dust screen model'
    t_merge['v_0'].description = 'Line of sight stellar velocity'
    t_merge['v_d'].description = 'Line of sight stellar velocity dispersion'
    t_merge['adev'].description = 'Mean model deviation'
    t_merge['nl_clip'].description = '% of wavelengths clipped by fitting algorithm'
    t_merge['chi2'].description = 'Fit statistic'
    t_merge['zones'].description = 'Voronoi segmentation zones'
    t_merge['sn'].description = 'SNR in individual pixels'
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
