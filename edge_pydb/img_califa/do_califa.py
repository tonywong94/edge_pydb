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

def do_califa(outname='NGC4047', gallist=['NGC4047'], linelbl='co', seq='smo7',
             hexgrid=False, allpix=False, debug=False, califa_natv='fitsdata',
             califa_smo='fitsdata', comom='../img_comom/fitsdata', nfiles=5, 
             astrom='x', ortpar='edge_leda.csv', distpar='edge_califa.csv',
             distcol='caDistP3d', discard_cdmatrix=False):
    """
    Extract Pipe3D products into an HDF5 database.  This script assumes
    there are 5 native resolution and 5 smoothed resolution files per galaxy.

    Parameters
    ----------
    outname : str
        Prefix of the output filename
    gallist : list of str
        List of galaxy names
    linelbl : str
        Identifier for reference line in the CO FITS filenames
    seq : str
        Identifier, generally to indicate smoothing resolution for CO
    hexgrid : boolean
        True to sample on a hexagonal grid (experimental)
    allpix : boolean
        True to dump every pixel, otherwise every 3rd pixel in x and y is used.
    debug : boolean
        True to generate some additional output
    califa_natv : str
        Path to the directory where native res CALIFA FITS files reside
    califa_smo : str
        Path to the directory where smoothed res CALIFA FITS files reside
    comom : str
        Path to the directory where CO moments FITS files reside
    nfiles : int
        Number of Pipe3D files per galaxy.  Should be 5 (old) or 1 (packed).
    astrom : str
        String at start of filename for native resolution images with astrometry.
        This is ignored in nfiles=1.
    ortpar : filename
        Name of the EdgeTable which has LEDA orientation parameters for the sample
    distpar : filename
        Name of the EdgeTable which has distances for converting \Sigma_*.
    distcol : str
        Name of the distance column in 'distpar' to use.  Default is 'caDistP3d'
        taken from 'DL' column in get_proc_elines_CALIFA.csv.
    discard_cdmatrix : boolean
        True to disregard CD matrix in CALIFA files.  Use with care since this
        relies on the CDELT1 and CDELT2 being correct.
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
        rglist = []
        smlist = []

        if len(gallist) == 0:
            raise RuntimeError('Error: gallist is empty!')

        for gal in gallist:
            print('\nWorking on galaxy {} product {} nsel={}'.format(
                gal, prod, nsel))

            # Generate output header using CO astrometry
            cofile = os.path.join(comom,
                    gal+'.'+linelbl+'.'+seq+'_dil.snrpk.fits.gz')
            if not os.path.exists(cofile):
                print('####### Cannot find',cofile)
                continue
            cohd = fits.getheader(cofile)
            # CALIFA files with x in name have optical astrometry
            if nfiles == 5:
                cafile = os.path.join(califa_natv,
                    astrom+leadstr[i_prod]+gal+tailstr[i_prod])
            else:
                cafile = os.path.join(califa_natv, gal+'.Pipe3D.cube.fits.gz')
            if not os.path.exists(cafile):
                print('####### Cannot find',cafile)
                continue          
            if nfiles == 5:
                hdus = fits.open(cafile, ignore_missing_end=True)
                cahd = hdus[0].header
            else:
                hdus = fits.open(cafile)
                # The header for the selected extension
                cahd = hdus[hdus.index_of(prod)].header
                # Blanking of CTYPE3 so that fitsextract treats as pseudocube
                cahd['CTYPE3'] = ''
                # Use HDU 0 'ORG_HDR' astrometry when possible
                cahd0 = hdus[0].header
                for key in cdkeys+wcskeys:
                    if key in cahd0.keys():
                        cahd[key] = cahd0[key]
                # Set CDELT3 to 1 since that will be its value in template
                for key in ['CDELT3', 'CD3_3']:
                	if key in cahd.keys():
                		cahd[key] = 1.
            # Copy the CALIFA header and replace wcskeys with CO values
            outhd = cahd.copy()
            for key in dimkeys+wcskeys:
                if key in cohd.keys():
                    outhd[key] = cohd[key]
            # Need to discard CD matrix which would override the new wcskeys
            if 'CDELT1' in cohd.keys() and 'CDELT2' in cohd.keys():
                for key in cdkeys:
                    if key in outhd.keys():
                        del outhd[key]
            # Optionally discard CD matrix in CALIFA files and fall back on CDELTs
            if discard_cdmatrix:
                for key in cdkeys:
                    if key in cahd.keys():
                        del cahd[key]

            # First process the native resolution file (tab0) with astrometry
            if nfiles == 5:
                hdu = hdus[0]
            else:
                hdu = hdus[hdus.index_of(prod)]
            if debug:
                print('\nINPUT',WCS(hdu.header))
                print('\nOUTPUT',WCS(outhd))
            newim = reproject_interp(hdu, outhd, order=0, return_footprint=False)
            nz = newim.shape[0]
            if debug:
                print('nz=',nz)
                if nfiles == 5:
                	fits.writeto(cafile.replace('fits','rg.fits'), newim, outhd, overwrite=True)
            rglabels = [s+'_rg' for s in labels]
            # Add smoothed Ha and Hb columns for extinction estimates
            if prod == 'ELINES' or prod == 'flux_elines':
                kernel = Gaussian2DKernel(3)
                if prod == 'ELINES':
                    hb_idx = 5
                    ha_idx = 6
                    rglabels += ['Hbeta_sm3_rg', 'Halpha_sm3_rg']
                    outhd['DESC_20'] = ' Hbeta after 3as smooth'
                    outhd['DESC_21'] = ' Halpha after 3as smooth'
                else:
                    hb_idx = 28
                    ha_idx = 45
                    rglabels += ['flux_Hbeta_sm3_rg', 'flux_Halpha_sm3_rg']
                    outhd['NAME408'] = ' Hbeta after 3as smooth'
                    outhd['NAME409'] = ' Halpha after 3as smooth'
                hb_conv = convolve(newim[hb_idx,:,:], kernel, preserve_nan=True)
                ha_conv = convolve(newim[ha_idx,:,:], kernel, preserve_nan=True)
                newim = np.concatenate((newim, hb_conv[np.newaxis], ha_conv[np.newaxis]))
                if len(zsel) == default_len:
                    zsel = list(zsel) + [nz, nz+1]
                if len(units) == default_len:
                    units += ['10^-16 erg cm^-2 s^-1', '10^-16 erg cm^-2 s^-1']
            tab0 = fitsextract(newim, header=outhd, keepnan=True, stride=stride, 
                               bunit=units, col_lbl=rglabels, zselect=zsel, 
                               ra_gc=15*orttbl.loc[gal]['ledaRA'],
                               dec_gc=orttbl.loc[gal]['ledaDE'], 
                               pa=orttbl.loc[gal]['ledaPA'],
                               inc=orttbl.loc[gal]['ledaAxIncl'], 
                               ortlabel='LEDA', first=True, use_hexgrid=hexgrid)
            gname = Column([np.string_(gal)]*len(tab0), name='Name', 
                           description='Galaxy Name')
            tab0.add_column(gname, index=0)
            rglist.append(tab0)

            # Then process the smoothed file (tab1)
            smofile = os.path.join(califa_smo,
                                   leadstr[i_prod]+gal+tailstr[i_prod])
            hdu = fits.open(smofile, ignore_missing_end=True)[0]
            hdu.header = cahd
            newim = reproject_interp(hdu, outhd, order=0, return_footprint=False)
            if debug:
                fits.writeto(smofile.replace('fits','sm.fits'), newim, outhd, overwrite=True)
            smlabels = [s+'_sm' for s in labels]
            # Add smoothed Ha and Hb for extinction estimates
            if prod == 'ELINES' or prod == 'flux_elines':
                kernel = Gaussian2DKernel(5)
                if prod == 'ELINES':
                    hb_idx = 5
                    ha_idx = 6
                    smlabels += ['Hbeta_sm5_sm', 'Halpha_sm5_sm']
                    outhd['DESC_20'] = ' Hbeta after 5as smooth'
                    outhd['DESC_21'] = ' Halpha after 5as smooth'
                else:
                    hb_idx = 28
                    ha_idx = 45
                    smlabels += ['flux_Hbeta_sm5_sm', 'flux_Halpha_sm5_sm']
                    outhd['NAME408'] = ' Hbeta after 5as smooth'
                    outhd['NAME409'] = ' Halpha after 5as smooth'
                hb_conv = convolve(newim[hb_idx,:,:], kernel, preserve_nan=True)
                ha_conv = convolve(newim[ha_idx,:,:], kernel, preserve_nan=True)
                newim = np.concatenate((newim, hb_conv[np.newaxis], ha_conv[np.newaxis]))
            tab1 = fitsextract(newim, header=outhd, keepnan=True, stride=stride, 
                               bunit=units, col_lbl=smlabels, zselect=zsel, 
                               ra_gc=15*orttbl.loc[gal]['ledaRA'],
                               dec_gc=orttbl.loc[gal]['ledaDE'], 
                               pa=orttbl.loc[gal]['ledaPA'],
                               inc=orttbl.loc[gal]['ledaAxIncl'], 
                               ortlabel='LEDA', first=True, use_hexgrid=hexgrid)
            gname = Column([np.string_(gal)]*len(tab1), name='Name', 
                           description='Galaxy Name')
            tab1.add_column(gname, index=0)
            smlist.append(tab1)
        
            # Add additional columns
            if prod == 'ELINES' or prod == 'flux_elines':
                if prod == 'ELINES':
                    prfx = ''
                else:
                    prfx = 'flux_'
                #
                # Native resolution
                # sfr0 is SFR from Halpha without extinction correction
                sfr0_rg = sfr_ha(tab0[prfx+'Halpha_rg'], imf='salpeter', 
                                 name=prfx+'sigsfr0_rg')
                e_sfr0_rg = Column(sfr0_rg *
                    abs(tab0['e_'+prfx+'Halpha_rg']/tab0[prfx+'Halpha_rg']), 
                    name='e_'+prfx+'sigsfr0_rg', dtype='f4', unit=sfr0_rg.unit,
                    description='error of uncorrected SFR surface density')
                tab0.add_columns([sfr0_rg, e_sfr0_rg])
                # Balmer decrement corrected SFR
                sfr_rg, sfrext_rg, e_sfr_rg, e_sfrext_rg = sfr_ha(tab0[prfx+'Halpha_rg'], 
                    flux_hb=tab0[prfx+'Hbeta_rg'], e_flux_ha=tab0['e_'+prfx+'Halpha_rg'],
                    e_flux_hb=tab0['e_'+prfx+'Hbeta_rg'], imf='salpeter', 
                    name=prfx+'sigsfr_corr_rg')
                tab0.add_columns([sfr_rg, e_sfr_rg, sfrext_rg, e_sfrext_rg])
                # Halpha extinction and SFR after 3" smoothing and clipping
                A_Ha3_rg = Column(get_AHa(tab0[prfx+'Halpha_sm3_rg'], 
                            tab0[prfx+'Hbeta_sm3_rg'], np.log10), 
                            name=prfx+'AHa_smooth3_rg', dtype='f4', unit='mag',
                            description='Ha extinction after 3as smooth')
                clip = ((tab0[prfx+'Halpha_sm3_rg'] < hacut) | 
                        (tab0[prfx+'Hbeta_sm3_rg'] < hbcut) | 
                        (A_Ha3_rg > ahahi) | (A_Ha3_rg < ahalo))
                sfr3_rg = Column(sfr0_rg * 10**(0.4*A_Ha3_rg),
                            name=prfx+'sigsfr_adopt_rg', dtype='f4', unit=sfr0_rg.unit,
                            description='smooth+clip BD corrected SFR surface density')
                sfr3_rg[clip] = sfr0_rg[clip]
                # A_Ha3_rg[clip] = np.nan
                tab0.add_columns([A_Ha3_rg, sfr3_rg])
                #
                # Smoothed resolution
                # sfr0 is SFR from Halpha without extinction correction
                sfr0_sm = sfr_ha(tab1[prfx+'Halpha_sm'], imf='salpeter', 
                                 name=prfx+'sigsfr0_sm')
                e_sfr0_sm = Column(sfr0_sm * 
                    abs(tab1['e_'+prfx+'Halpha_sm']/tab1[prfx+'Halpha_sm']), 
                    name='e_'+prfx+'sigsfr0_sm', dtype='f4', unit=sfr0_sm.unit,
                    description='error of uncorrected SFR surface density')
                tab1.add_columns([sfr0_sm, e_sfr0_sm])
                # Balmer decrement corrected SFR
                sfr_sm, sfrext_sm, e_sfr_sm, e_sfrext_sm = sfr_ha(tab1[prfx+'Halpha_sm'],
                    flux_hb=tab1[prfx+'Hbeta_sm'], e_flux_ha=tab1['e_'+prfx+'Halpha_sm'],
                    e_flux_hb=tab1['e_'+prfx+'Hbeta_sm'], imf='salpeter', 
                    name=prfx+'sigsfr_corr_sm')
                tab1.add_columns([sfr_sm, e_sfr_sm, sfrext_sm, e_sfrext_sm])
                # Halpha extinction and SFR after 5" smoothing and clipping
                A_Ha5_sm = Column(get_AHa(tab1[prfx+'Halpha_sm5_sm'], 
                            tab1[prfx+'Hbeta_sm5_sm'], np.log10), 
                            name=prfx+'AHa_smooth5_sm', dtype='f4', unit='mag',
                            description='Ha extinction after 5as smooth')
                clip = ((tab1[prfx+'Halpha_sm5_sm'] < hacut) | 
                        (tab1[prfx+'Hbeta_sm5_sm'] < hbcut) | 
                        (A_Ha5_sm > ahahi) | (A_Ha5_sm < ahalo))
                sfr5_sm = Column(sfr0_sm * 10**(0.4*A_Ha5_sm),
                            name=prfx+'sigsfr_adopt_sm', dtype='f4', unit=sfr0_rg.unit,
                            description='smooth+clip BD corrected SFR surface density')
                sfr5_sm[clip] = sfr0_sm[clip]
                # A_Ha5_sm[clip] = np.nan
                tab1.add_columns([A_Ha5_sm, sfr5_sm])
                #
                # BPT requires flux_elines since EW(Ha) is part of classification
                if prod == 'flux_elines':
                    BPT0, BPT0sf, p_BPT0 = bpt_type(tab0, ext='_rg', name='BPT_rg', prob=True)
                    tab0.add_columns([BPT0, p_BPT0, BPT0sf])
                    BPT1, BPT1sf, p_BPT1 = bpt_type(tab1, ext='_sm', name='BPT_sm', prob=True)
                    tab1.add_columns([BPT1, p_BPT1, BPT1sf])
                    #
                    zoh0, zoherr0 = ZOH_M13(tab0, ext='_rg', name='ZOH_rg', err=True)
                    tab0.add_columns([zoh0, zoherr0])
                    zoh1, zoherr1 = ZOH_M13(tab1, ext='_sm', name='ZOH_sm', err=True)
                    tab1.add_columns([zoh1, zoherr1])
            elif prod == 'SSP':
                # For stellar surface density we need distance
                star0 = stmass_pc2(tab0['mass_ssp_rg'], dz=tab0['cont_dezon_rg'],
                                dist=disttbl.loc[gal][distcol], name='sigstar_rg')
                star1 = stmass_pc2(tab1['mass_ssp_sm'], dz=tab1['cont_dezon_sm'],
                                dist=disttbl.loc[gal][distcol], name='sigstar_sm')
                avstar0 = stmass_pc2(tab0['mass_Avcor_ssp_rg'], dz=tab0['cont_dezon_rg'],
                                dist=disttbl.loc[gal][distcol], name='sigstar_Avcor_rg')
                avstar0.description += ' dust corrected'
                avstar1 = stmass_pc2(tab1['mass_Avcor_ssp_sm'], dz=tab1['cont_dezon_sm'],
                                dist=disttbl.loc[gal][distcol], name='sigstar_Avcor_sm')
                avstar1.description += ' dust corrected'
                ferr0 = Column(abs(tab0['e_medflx_ssp_rg']/tab0['medflx_ssp_rg']), 
                    name='fe_medflx_rg', dtype='f4', unit='fraction',
                    description='fractional error in continuum flux')
                ferr1 = Column(abs(tab1['e_medflx_ssp_sm']/tab1['medflx_ssp_sm']), 
                    name='fe_medflx_sm', dtype='f4', unit='fraction',
                    description='fractional error in continuum flux')
                tab0.add_columns([star0, avstar0, ferr0])
                tab1.add_columns([star1, avstar1, ferr1])

        if len(rglist) > 0:
            rg_merge = vstack(rglist)
        rg_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
        if debug:
            print(rg_merge.colnames)
            print('There are',len(rg_merge),'rows in native table')

        if len(smlist) > 0:
            sm_merge = vstack(smlist)
        sm_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
        if debug:
            print(sm_merge.colnames)
            print('There are',len(sm_merge),'rows in smoothed table')

        if prod == prodtype[0]:
            rg_merge.write(outname+'.pipe3d.hdf5', path=prod+'_rg', overwrite=True, 
                    serialize_meta=True, compression=True)
        else:
            rg_merge.write(outname+'.pipe3d.hdf5', path=prod+'_rg', append=True, 
                    serialize_meta=True, compression=True)
        sm_merge.write(outname+'.pipe3d.hdf5', path=prod+'_sm', append=True, 
                    serialize_meta=True, compression=True)
    return


if __name__ == "__main__":
    # NGC4047 only
    do_califa()
    do_califa(hexgrid=True, outname='NGC4047_hex')
    # All EDGE125 galaxies
    gallist = [os.path.basename(file).split('.')[0] for file in 
               sorted(glob.glob('fitsdata/[A-Z]*.SSP.cube.fits.gz'))]
    do_califa(gallist=gallist, outname='edge')
    do_califa(gallist=gallist, hexgrid=True, outname='edge_hex')
    # EDGE125 allpix data
    do_califa(gallist=gallist, outname='edge_allpix', allpix=True)
    # ACA galaxies, native resolution
    gallist = [os.path.basename(file).split('.')[0] for file in 
               sorted(glob.glob('aca_packed/*.Pipe3D.cube.fits.gz'))]
    do_califa(gallist=gallist, outname='edge_aca', seq='natv', linelbl='co21', 
    		califa_natv='aca_packed', califa_smo='aca_conv9', comom='../img_comom/acadata', 
    		nfiles=1)

