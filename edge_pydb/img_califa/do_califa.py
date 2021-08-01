#!/usr/bin/env python

# Combine the CALIFA data into binary tables.
# We now use Salpeter IMF for SFR calculation (as of 05-Dec-2020)

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

def do_califa(outfile='NGC4047.pipe3d.hdf5', gallist=['NGC4047'], 
             fitsdir='fits_natv_edge', comomdir=None, colabel='co.smo7',
             ext='', nsm=2, ortpar='edge_leda.csv', distpar='edge_califa.csv',
             distcol='caDistP3d', hexgrid=False, allpix=False, debug=False, 
             discard_cdmatrix=False, append=True, overwrite=True):
    """
    Extract Pipe3D products into an HDF5 database.  This script assumes
    there are 5 Pipe3D output files per galaxy.

    Parameters
    ----------
    outfile : str
        Name of the output filename.  Appended to if it exists.
    gallist : list of str
        List of galaxy names
    fitsdir : str
        Path to the directory where CALIFA FITS files reside
    comomdir : str
        Path to the directory where CO moments FITS files reside.  Default is
        process CALIFA data only without regridding to match CO.
    colabel : str
        Identifier for reference line in the CO FITS filenames
    ext : str
        Suffix to add to column names, e.g. '_sm'
    nsm : int
        Stddev of Gaussian smoothing kernel in pixels to use for Balmer decrement
        noise reduction.
    ortpar : filename
        Name of the EdgeTable which has LEDA orientation parameters for the sample
    distpar : filename
        Name of the EdgeTable which has distances for converting \Sigma_*.
    distcol : str
        Name of the distance column in 'distpar' to use.  Default is 'caDistP3d'
        taken from 'DL' column in get_proc_elines_CALIFA.csv.
    hexgrid : boolean
        True to sample on a hexagonal grid (experimental)
    allpix : boolean
        True to dump every pixel, otherwise every 3rd pixel in x and y is used.
    debug : boolean
        True to generate some additional output
    discard_cdmatrix : boolean
        True to disregard CD matrix in CALIFA files.  Use with care since this
        relies on the CDELT1 and CDELT2 being correct.
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
    orttbl.add_index('Name') 

    # Get the distance from the CALIFA table
    disttbl = EdgeTable(distpar)
    disttbl.add_index('Name')

    # Read the FITS data
    # The columns to save are defined in fitsextract.py
    prodtype = ['ELINES', 'SFH', 'SSP', 'indices', 'flux_elines']
    leadstr  = ['', '', '', 'indices.CS.', 'flux_elines.']
    tailstr  = ['.ELINES', '.SFH', '.SSP', '', '']
    tailstr  = [s+'.cube.fits.gz' for s in tailstr]

    for i_prod, prod in enumerate(prodtype):
        zsel, labels, units, nsel = getlabels(prod)
        default_len = len(zsel)
        tlist = []

        if len(gallist) == 0:
            raise RuntimeError('Error: gallist is empty!')

        for gal in gallist:
            print('\nWorking on galaxy {} product {} nsel={}'.format(
                gal, prod, nsel))

            # Read in Pipe3D output
            cafile = os.path.join(fitsdir,leadstr[i_prod]+gal+tailstr[i_prod])
            if not os.path.exists(cafile):
                print('####### Cannot find',cafile)
                continue          
            hdu = fits.open(cafile, ignore_missing_end=True)[0]
            cahd = hdu.header.copy()
            # Blanking of CTYPE3 so that fitsextract treats cubes as pseudocubes
            cahd['CTYPE3'] = ''
            # Set CDELT3 to 1 since that will be its value in template
            for key in ['CDELT3', 'CD3_3']:
                if key in cahd.keys():
                    cahd[key] = 1.

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
                    fits.writeto(cafile.replace('.fits','.rg.fits'), newim, cahd, 
                                 overwrite=True)
            else:
                newim = hdu.data

            # Set up output table
            nz = newim.shape[0]
            if debug:
                print('nz=',nz)
            col_lbl = [s+ext for s in labels]

            # Add smoothed Ha and Hb columns for extinction estimates
            if prod == 'ELINES' or prod == 'flux_elines':
                kernel = Gaussian2DKernel(nsm)
                if prod == 'ELINES':
                    hb_idx = 5
                    ha_idx = 6
                    col_lbl += ['Hbeta_sm'+str(nsm)+ext, 'Halpha_sm'+str(nsm)+ext]
                    cahd['DESC_20'] = ' Hbeta after {} pix smooth'.format(str(nsm))
                    cahd['DESC_21'] = ' Halpha after {} pix smooth'.format(str(nsm))
                else:
                    hb_idx = 28
                    ha_idx = 45
                    col_lbl += ['flux_Hbeta_sm'+str(nsm)+ext, 'flux_Halpha_sm'+str(nsm)+ext]
                    cahd['NAME408'] = ' Hbeta after {} pix smooth'.format(str(nsm))
                    cahd['NAME409'] = ' Halpha after {} pix smooth'.format(str(nsm))
                hb_conv = convolve(newim[hb_idx,:,:], kernel, preserve_nan=True)
                ha_conv = convolve(newim[ha_idx,:,:], kernel, preserve_nan=True)
                newim = np.concatenate((newim, hb_conv[np.newaxis], ha_conv[np.newaxis]))
                if len(zsel) == default_len:
                    zsel = list(zsel) + [nz, nz+1]
                if len(units) == default_len:
                    units += ['10^-16 erg cm^-2 s^-1', '10^-16 erg cm^-2 s^-1']
            tab0 = fitsextract(newim, header=cahd, keepnan=True, stride=stride, 
                               bunit=units, col_lbl=col_lbl, zselect=zsel, 
                               ra_gc=15*orttbl.loc[gal]['ledaRA'],
                               dec_gc=orttbl.loc[gal]['ledaDE'], 
                               pa=orttbl.loc[gal]['ledaPA'],
                               inc=orttbl.loc[gal]['ledaAxIncl'], 
                               ortlabel='LEDA', first=True, use_hexgrid=hexgrid)
            gname = Column([np.string_(gal)]*len(tab0), name='Name', 
                           description='Galaxy Name')
            tab0.add_column(gname, index=0)
        
            # Add additional columns
            if prod == 'ELINES' or prod == 'flux_elines':
                if prod == 'ELINES':
                    prfx = ''
                else:
                    prfx = 'flux_'

                # sfr0 is SFR from Halpha without extinction correction
                sfr0 = sfr_ha(tab0[prfx+'Halpha'+ext], imf='salpeter', 
                                 name=prfx+'sigsfr0'+ext)
                e_sfr0 = Column(sfr0 *
                    abs(tab0['e_'+prfx+'Halpha'+ext]/tab0[prfx+'Halpha'+ext]), 
                    name='e_'+prfx+'sigsfr0'+ext, dtype='f4', unit=sfr0.unit,
                    description='error of uncorrected SFR surface density')
                tab0.add_columns([sfr0, e_sfr0])

                # Balmer decrement corrected SFR
                sfr_cor, A_Ha, e_sfr_cor, e_A_Ha = sfr_ha(
                            tab0[prfx+'Halpha'+ext], 
                            flux_hb=tab0[prfx+'Hbeta'+ext], 
                            e_flux_ha=tab0['e_'+prfx+'Halpha'+ext],
                            e_flux_hb=tab0['e_'+prfx+'Hbeta'+ext], 
                            imf='salpeter', 
                            name=prfx+'sigsfr_corr'+ext)
                tab0.add_columns([sfr_cor, e_sfr_cor, A_Ha, e_A_Ha])

                # Halpha extinction and SFR after smoothing and clipping
                A_Ha_smo = Column(get_AHa(tab0[prfx+'Halpha_sm'+str(nsm)+ext], 
                            tab0[prfx+'Hbeta_sm'+str(nsm)+ext], np.log10), 
                            name=prfx+'AHa_smooth'+str(nsm)+ext, dtype='f4', unit='mag',
                            description='Ha extinction after {} pix smooth'.format(str(nsm)))
                clip = ((tab0[prfx+'Halpha_sm'+str(nsm)+ext] < hacut) | 
                        (tab0[prfx+'Hbeta_sm'+str(nsm)+ext] < hbcut) | 
                        (A_Ha_smo > ahahi) | (A_Ha_smo < ahalo))
                sfr_smo = Column(sfr0 * 10**(0.4*A_Ha_smo),
                            name=prfx+'sigsfr_adopt'+ext, dtype='f4', unit=sfr0.unit,
                            description='smooth+clip BD corrected SFR surface density')
                sfr_smo[clip] = sfr0[clip]
                # This would clean up the smoothed A_Ha map but obscure the selection
                # A_Ha_smo[clip] = np.nan
                tab0.add_columns([A_Ha_smo, sfr_smo])

                # BPT requires flux_elines since EW(Ha) is part of classification
                if prod == 'flux_elines':
                    BPT0, BPT0sf, p_BPT0 = bpt_type(tab0, ext=ext, name='BPT'+ext, prob=True)
                    tab0.add_columns([BPT0, p_BPT0, BPT0sf])
                    #
                    zoh0, zoherr0 = ZOH_M13(tab0, ext=ext, name='ZOH'+ext, err=True)
                    tab0.add_columns([zoh0, zoherr0])

            elif prod == 'SSP':
                # For stellar surface density we need distance
                star0 = stmass_pc2(tab0['mass_ssp'+ext], dz=tab0['cont_dezon'+ext],
                                dist=disttbl.loc[gal][distcol], name='sigstar'+ext)
                avstar0 = stmass_pc2(tab0['mass_Avcor_ssp'+ext], dz=tab0['cont_dezon'+ext],
                                dist=disttbl.loc[gal][distcol], name='sigstar_Avcor'+ext)
                avstar0.description += ' dust corrected'
                ferr0 = Column(abs(tab0['e_medflx_ssp'+ext]/tab0['medflx_ssp'+ext]), 
                    name='fe_medflx'+ext, dtype='f4', unit='fraction',
                    description='fractional error in continuum flux')
                tab0.add_columns([star0, avstar0, ferr0])

            tlist.append(tab0)

        if len(tlist) > 0:
            t_merge = vstack(tlist)
        t_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
        if debug:
            print(t_merge.colnames)
            print('There are',len(t_merge),'rows in merged table')

        if prod == prodtype[0]:
            t_merge.write(outfile, path=prod+ext, overwrite=overwrite, 
                    append=append, serialize_meta=True, compression=True)
        else:
            t_merge.write(outfile, path=prod+ext, overwrite=overwrite, 
                    append=True, serialize_meta=True, compression=True)
    return


if __name__ == "__main__":
    # NGC4047, CALIFA only
    do_califa(outfile='NGC4047_allpix.pipe3d.hdf5', append=False, allpix=True)
    # NGC4047, append to 2d_smo7
    do_califa(outfile='../img_comom/NGC4047.2d_smo7.hdf5', fitsdir='fits_smo7_edge',
              comomdir='../img_comom/fitsdata', ext='_sm', nsm=3, append=True)
    do_califa(outfile='../img_comom/NGC4047_hex.2d_smo7.hdf5', fitsdir='fits_smo7_edge',
              comomdir='../img_comom/fitsdata', ext='_sm', nsm=3, append=True,
              hexgrid=True)
    # All EDGE125 galaxies, CALIFA only
    gallist = [os.path.basename(file).split('.')[0] for file in 
               sorted(glob.glob('fits_smo7_edge/[A-Z]*.SSP.cube.fits.gz'))]
    do_califa(gallist=gallist, outfile='edge_carma_allpix.pipe3d.hdf5', 
              append=False, allpix=True)
    # All EDGE125 galaxies, append to 2d_smo7
    do_califa(gallist=gallist, outfile='../img_comom/edge_carma.2d_smo7.hdf5', 
              fitsdir='fits_smo7_edge', comomdir='../img_comom/fitsdata', 
              ext='_sm', nsm=3, append=True)
    do_califa(gallist=gallist, outfile='../img_comom/edge_carma_allpix.2d_smo7.hdf5', 
              fitsdir='fits_smo7_edge', comomdir='../img_comom/fitsdata', 
              ext='_sm', nsm=3, append=True, allpix=True)
    # EDGE125 hexgrid
#     do_califa(gallist=gallist, outfile='../img_comom/edge_carma_hex.2d_smo7.hdf5', 
#               fitsdir='fits_smo7_edge', comomdir='../img_comom/fitsdata', 
#               ext='_sm', nsm=3, hexgrid=True, append=True)
    # ACA galaxies, 9" resolution
    gallist = [os.path.basename(file).split('.')[0] for file in 
               sorted(glob.glob('fits_smo9_aca/[A-Z]*.SSP.cube.fits.gz'))]
    do_califa(gallist=gallist, outfile='edge_aca_allpix.pipe3d.hdf5', 
              fitsdir='fits_natv_aca', append=False, allpix=True)
    do_califa(gallist=gallist, outfile='../img_comom/edge_aca.2d_smo9.hdf5', 
              colabel='CO.smo9', comomdir='../img_comom/acadata', fitsdir='fits_smo9_aca', 
              ext='_sm', nsm=4, append=True)
    do_califa(gallist=gallist, outfile='../img_comom/edge_aca_allpix.2d_smo9.hdf5',
              colabel='CO.smo9', comomdir='../img_comom/acadata', fitsdir='fits_smo9_aca', 
              ext='_sm', nsm=4, append=True, allpix=True)
