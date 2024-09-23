#!/usr/bin/env python

# Combine the CALIFA data into binary tables.
# We now use Salpeter IMF for SFR calculation (as of 05-Dec-2020)

from datetime import datetime
import glob
import os
import numpy as np
from astropy import units as u
from astropy.table import Table, Column, join, vstack, hstack
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from reproject import reproject_interp
from edge_pydb import EdgeTable
from edge_pydb.conversion import stmass_pc2, sfr_ha, ZOH_M13, bpt_type, get_AHa
from CO_conversion_factor.alphaCO import predict_alphaCO10_B13, predict_alphaCO_SL24
from edge_pydb.fitsextract import fitsextract, getlabels
from pyFIT3D.modelling.stellar import SSPModels
np.seterr(divide='ignore', invalid='ignore')

def do_pipe3d(outfile='NGC4047.pipe3d.hdf5', gallist=None, fitsdir=None, 
              p3dstruct='califa', packed=True, comomdir=None, 
              cotempl='GNAME.co_dil.snrpk.fits.gz', ssptable='gsd01_156.fits', 
              interp_order=1, ext='', nsm=2, ortpar='edge_leda.csv', 
              distpar='edge_califa.csv', ortlabel='LEDA', coln_ra='ledaRA', 
              coln_dc='ledaDE', coln_pa='ledaPA', coln_inc='ledaAxIncl',
              coln_dmpc='caDistP3d', hexgrid=False, allpix=False, debug=False, 
              keepnan=True, blankzero=True, prob=True, discard_cdmatrix=False, 
              append=True, overwrite=True, 
              prodtype=['ELINES', 'SFH', 'SSP', 'indices', 'flux_elines'],
              leadstr=['', '', '', 'indices.CS.', 'flux_elines.'],
              tailstr=['.ELINES','.SFH','.SSP','',''], tailx='.cube.fits.gz'):
    """
    Extract Pipe3D products into an HDF5 database.  This script now handles packed
    Pipe3D files and MaNGA data.
    
    Explanation of WCS handling: the Pipe3D header ('p3dhd') taken from the
    packed file or the flux_elines file is copied into 'cawcshd' which
    has 3rd axis keywords tweaked to allow fitsextract to treat the cubes as 
    pseudo-cubes.  If a CO template is read in, its header ('cohd') is used to
    replace the WCS keywords (including NAXIS1 and NAXIS2) in 'cawcshd', which
    becomes the target header for reprojection.  'w_cahd' is a copy of 'cawcshd' 
    which has additional descriptor keywords from the particular Pipe3D products 
    (ELINES, etc.) included, and is provided to fitsextract.

    Parameters
    ----------
    outfile : filename
        Name of the output HDF5 filename.  Appended to if it exists.
    gallist : list of str
        List of galaxy names
    fitsdir : dirname
        Path to the directory where CALIFA FITS files reside
    p3dstruct : str
        Structure convention for the Pipe3D FITS file.
        'califa': ELINES has errors, flux_elines has 51 emission lines, SFH has 39
                  ages 4 mets and uncertainties, SSP lacks 'e_mass_ssp'
        'manga' : ELINES lacks errors, flux_elines has 57 emission lines, SFH has
                  39 ages 7 mets no uncertainties, SSP has 'e_mass_ssp'
    packed : boolean
        True if Pipe3D outputs for a galaxy are in a single multi-extension FITS file.
    comomdir : str
        Path to the directory where CO moments FITS files reside.  Default is None,
        i.e. process CALIFA data only without regridding to match CO.
    cotempl : str
        Name of CO file with the astrometry which will be adopted for the IFU sampling.
        Note that the CO template files should use the CDELT1 and CDELT2 keywords
        and not CD1_1 and CD2_2.
    ssptable : str
        Name of SSP models table for deriving mass-to-luminosity ratio
    interp_order : int
        Interpolation order for reproject.  Default is 1 (bilinear).  Use 0 for 
        nearest neighbor interpolation, where new values are drawn from existing ones.
    ext : str
        Suffix to add to column names, e.g. '_sm'
    nsm : int
        Stddev of Gaussian smoothing kernel in pixels to use for Balmer decrement
        noise reduction.
    ortpar : filename
        Name of the EdgeTable which has orientation parameters for the sample, 
        including center position, inclination, and position angle.  This can also 
        be given as the path to a regular astropy-compatible Table.
    distpar : filename
        Name of the EdgeTable which has distances for converting \Sigma_*.
    ortlabel : str
        String for labeling source of the orientation parameters.  Default is 'LEDA'.
    coln_ra : str
        Name of the RA column in 'ortpar' to use.  Default is 'ledaRA'.
    coln_dc : str
        Name of the DEC column in 'ortpar' to use.  Default is 'ledaDE'.
    coln_pa : str
        Name of the PA column in 'ortpar' to use.  Default is 'ledaPA'.
    coln_inc : str
        Name of the INC column in 'ortpar' to use.  Default is 'ledaAxIncl'.
    coln_dmpc : str
        Name of the distance column in 'distpar' to use.  Default is 'caDistP3d'
        taken from 'DL' column in get_proc_elines_CALIFA.csv.
    hexgrid : boolean
        True to sample on a hexagonal grid (experimental)
    allpix : boolean
        True to dump every pixel, otherwise every 3rd pixel in x and y is used.
    debug : boolean
        True to generate some additional output
    keepnan : boolean
        True to replace zeroes in the FITS file with NaN.  This is the default.
    blankzero : boolean
        True to replace zeroes in the FITS file with NaN.  This is the default.
    prob : boolean
        True to add BPT probability column to flux_elines table.
    discard_cdmatrix : boolean
        True to disregard CD matrix in CALIFA files.  Use with care since this
        relies on the CDELT1 and CDELT2 being correct.
    append : boolean
        True to append to an existing file.  This is the default (write to the
        file created by do_comom.py).
    overwrite : boolean
        True to overwrite existing tables.  This is the default (replace same table
        but do not delete other tables in the file).
    prodtype : list of str
        List of Pipe3D products in the order to be analyzed.
    leadstr : list of str
        File name convention for Pipe3D products.  These strings precede the galaxy name.
        Not needed if packed=True.
    tailstr : list of str
        File name convention for Pipe3D products.  These strings follow the galaxy name.
        Not needed if packed=True.
    tailx : str
        Common ending string for the unpacked file names.
        Not needed if packed=True.
    """
    if allpix:
        stride = [1,1,1]
    else:
        stride = [3,3,1]

    if gallist is None or len(gallist) == 0:
        raise RuntimeError('Error: gallist is empty!')

    # cuts for when to apply BD correction
    ahalo = 0       # mag
    ahahi = 6       # mag

    # FITS keywords important for astrometry
    wcskeys = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 
               'CDELT1', 'CDELT2']
    cdkeys = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CD1_3', 'CD2_3',
                'CD3_1', 'CD3_2', 'CD3_3']
    dimkeys = ['NAXIS1', 'NAXIS2']

    # Get the orientation parameters and distances from global tables
    try:
        orttbl = EdgeTable(ortpar)
    except:
        orttbl = Table.read(ortpar, format='ascii.ecsv')
    orttbl.add_index('Name') 
    try:
        disttbl = EdgeTable(distpar)
    except:
        disttbl = Table.read(distpar, format='ascii.ecsv')
    disttbl.add_index('Name') 

    # Required file for SFH lum to mass conversion
    models = SSPModels(ssptable)
    print('Number of model steps:',models.n_models)
    nlumcols = models.n_models

    tablist = []

    # BEGIN loop over galaxies
    for i_gal, gname in enumerate(gallist):
        print('\nCanonical galaxy name is {}'.format(gname))

        if gname not in orttbl['Name']:
            print('\nERROR: Did not find galaxy',gname,'in',ortpar)
            continue
        elif gname not in disttbl['Name']:
            print('\nERROR: Did not find galaxy',gname,'in',distpar)
            continue
        else:
            print("RA={:.3f} deg Dec={:.3f} deg PA={:.1f} deg Inc={:.1f} deg Dist={:.1f} Mpc".format(
                  orttbl.loc[gname][coln_ra],
                  orttbl.loc[gname][coln_dc],
                  orttbl.loc[gname][coln_pa],
                  orttbl.loc[gname][coln_inc],
                  disttbl.loc[gname][coln_dmpc]))

        # Read in Pipe3D output to get astrometry
        if packed:
            if p3dstruct=='manga':
                p3d_file = os.path.join(fitsdir, 'manga-'+gname+'.Pipe3D.cube.fits.gz')
            else:
                p3d_file = os.path.join(fitsdir, gname+'.Pipe3D.cube.fits.gz')
        else:
            p3d_file = os.path.join(fitsdir, 'flux_elines.'+gname+tailx)
        if not os.path.exists(p3d_file):
            print('####### Cannot find',p3d_file)
            continue          
        hdul = fits.open(p3d_file, ignore_missing_end=True)
        p3dhd = hdul[0].header
        cawcshd = p3dhd.copy()
        # Arcseconds per pixel in Pipe3D output
        pixsca = round(3600*WCS(cawcshd).pixel_scale_matrix[1][1],2) * u.arcsec
        print('The pixel scale is', pixsca)
        # Blanking of CTYPE3 so that fitsextract treats cubes as pseudocubes
        cawcshd['CTYPE3'] = ''
        # Set CDELT3 to 1 since that will be its value in template
        for key in ['CDELT3', 'CD3_3']:
            if key in cawcshd.keys():
                cawcshd[key] = 1.
            if key in p3dhd.keys():
                p3dhd[key] = 1.

        # Read in CO template and prepare target header
        if comomdir is not None:
            cofile = os.path.join(comomdir,cotempl.replace('GNAME',gname))
            if not os.path.exists(cofile):
                # Fudge for almaquest file naming convention
                gname2 = gname.replace('-','_')
                cofile = os.path.join(comomdir,cotempl.replace('GNAME',gname2))
                if not os.path.exists(cofile):
                    print('####### Cannot find',cofile)
                    continue
            cohd = fits.getheader(cofile)
            # Copy the CALIFA header and replace wcskeys with CO values
            for key in dimkeys+wcskeys:
                if key in cohd.keys():
                    cawcshd[key] = cohd[key]
            # Need to discard CD matrix in tgt hdr which would override new wcskeys
            if 'CDELT1' in cohd.keys() and 'CDELT2' in cohd.keys():
                for key in cdkeys:
                    if key in cawcshd.keys():
                        del cawcshd[key]
            # Optionally discard CD matrix in P3D files and fall back on CDELTs
            if discard_cdmatrix:
                for key in cdkeys:
                    if key in p3dhd.keys():
                        del p3dhd[key]
            if debug:
                print('\nINPUT',WCS(p3dhd))
                print('\nCO data',WCS(cohd))
                print('\nOUTPUT',WCS(cawcshd))

        # BEGIN loop over products
        for i_prod, prod in enumerate(prodtype):
            print('\nWorking on extension {}'.format(prod))
            zsel, labels, units, nsel, has_errors = getlabels(prod, p3dstruct=p3dstruct)
            default_len = len(zsel)

            # Read the extension header and data
            if packed:
                cahd = hdul[prod].header
                cadat = hdul[prod].data
            else:
                cafile = os.path.join(fitsdir,leadstr[i_prod]+gname+tailstr[i_prod]+tailx)
                hdu = fits.open(cafile, ignore_missing_end=True)[0]
                cahd = hdu.header
                cadat = hdu.data
            w_cahd = cawcshd.copy()
            # Extension header may have useful documentation
            desckeys = [key for key in list(cahd.keys()) if key.startswith('DESC')]
            for key in desckeys:
                w_cahd[key] = cahd[key]
            # Fixes issue where header 0 has no NAXIS1 or NAXIS2
            if comomdir is None:
                for key in ['NAXIS2', 'NAXIS1']:
                    if key in w_cahd.keys():
                        del w_cahd[key]
                    w_cahd.insert('NAXIS', (key, cahd[key]), after=True)
            if blankzero:
                cadat[cadat==0] = np.nan

            # Regrid to the CO template
            if comomdir is not None:
                newim = reproject_interp((cadat,p3dhd), WCS(cawcshd), order=interp_order,
                    shape_out=(cahd['NAXIS3'],cawcshd['NAXIS2'],cawcshd['NAXIS1']),
                    return_footprint=False)
                if debug:
                    print(WCS(cawcshd))
                    fits.writeto(gname+'.'+prod+'.rg.fits', newim, w_cahd, overwrite=True)
            else:
                newim = cadat
        
            # Set up output table
            nz = newim.shape[0]
            if debug:
                print('nz=',nz,'zsel=',zsel)
            col_lbl = [s+ext for s in labels]

            # Add smoothed Ha and Hb columns for extinction estimates
            if prod == 'ELINES' or prod == 'flux_elines':
                kernel = Gaussian2DKernel(nsm)
                if prod == 'ELINES':
                    hb_idx = 5
                    ha_idx = 6
                    if has_errors:
                        e_hb_idx = hb_idx + (nz-2)//2
                        e_ha_idx = ha_idx + (nz-2)//2
                    col_lbl += ['Hbeta_sm'+str(nsm)+ext, 'Halpha_sm'+str(nsm)+ext]
                    cahd['DESC_20'] = ' Hbeta after {} pix smooth'.format(str(nsm))
                    cahd['DESC_21'] = ' Halpha after {} pix smooth'.format(str(nsm))
                else:
                    hb_idx = 28
                    ha_idx = 45
                    if has_errors:
                        e_hb_idx = hb_idx + nz//2
                        e_ha_idx = ha_idx + nz//2
                    col_lbl += ['flux_Hbeta_sm'+str(nsm)+ext, 'flux_Halpha_sm'+str(nsm)+ext]
                # Estimate of noise following convolution, see 2021RNAAS...5...39K
#                 e_kernel = Gaussian2DKernel(nsm/np.sqrt(2))
#                 base_res = fwhm / np.sqrt(8*np.log(2))
#                 scalefac = 4 * np.pi * base_res**2 * nsm**2 / (base_res**2 + nsm**2)
#                 e_hb_conv = scalefac * convolve(newim[e_hb_idx,:,:], e_kernel, preserve_nan=True)
#                 e_ha_conv = scalefac * convolve(newim[e_ha_idx,:,:], e_kernel, preserve_nan=True)
                hb_conv = convolve(newim[hb_idx,:,:], kernel, preserve_nan=True)
                ha_conv = convolve(newim[ha_idx,:,:], kernel, preserve_nan=True)
                newim = np.concatenate((newim, hb_conv[np.newaxis], ha_conv[np.newaxis]))
                print('Reshaped newim:',newim.shape)
                if len(zsel) == default_len:
                    zsel = list(zsel) + [nz, nz+1]
                if len(units) == default_len:
                    units += ['10^-16 erg cm^-2 s^-1', '10^-16 erg cm^-2 s^-1']

            tab0 = fitsextract(newim, header=w_cahd, keepnan=keepnan, stride=stride, 
                               bunit=units, col_lbl=col_lbl, zselect=zsel, 
                               ra_gc=orttbl.loc[gname][coln_ra],
                               dec_gc=orttbl.loc[gname][coln_dc], 
                               pa=orttbl.loc[gname][coln_pa],
                               inc=orttbl.loc[gname][coln_inc], 
                               ortlabel=ortlabel, first=True, use_hexgrid=hexgrid)
            if debug:
                print(tab0.colnames)
            gname_coln = Column([np.string_(gname)]*len(tab0), name='Name', 
                           description='Galaxy Name')
            tab0.add_column(gname_coln, index=0)
        
            # Add additional columns
            if prod == 'ELINES' or prod == 'flux_elines':
                if prod == 'ELINES':
                    prfx = ''
                else:
                    prfx = 'flux_'
                    # Provide labels for flux_elines columns
                    for linecol in labels:
                        if linecol.startswith('e_'):
                            linetype = linecol.split('_')[1]
                            linename = linecol.split('_')[2]
                            prelbl = 'error in '
                        else:
                            linetype = linecol.split('_')[0]
                            linename = linecol.split('_')[1]
                            prelbl = ''
                        if linetype == 'flux':
                            suffix = 'intensity'
                        elif linetype == 'vel':
                            suffix = 'velocity'
                        elif linetype == 'disp':
                            suffix = 'velocity dispersion'
                        elif linetype == 'EW':
                            suffix = 'equivalent width'
                        tab0[linecol+ext].description=prelbl+linename+' '+suffix
                    tab0['flux_Hbeta_sm'+str(nsm)+ext].description=\
                         'Hbeta intensity after {} pix smooth'.format(str(nsm))
                    tab0['flux_Halpha_sm'+str(nsm)+ext].description=\
                         'Halpha intensity after {} pix smooth'.format(str(nsm))

                # sfr0 is SFR from Halpha without extinction correction
                sfr0 = sfr_ha(tab0[prfx+'Halpha'+ext], imf='salpeter', 
                                 pixsca=pixsca, name=prfx+'sigsfr0'+ext)
                tab0.add_column(sfr0)

                # Balmer decrement corrected SFR
                if has_errors:
                    e_sfr0 = Column(sfr0 *
                        abs(tab0['e_'+prfx+'Halpha'+ext]/tab0[prfx+'Halpha'+ext]), 
                        name='e_'+prfx+'sigsfr0'+ext, dtype='f4', unit=sfr0.unit,
                        description='error of uncorrected SFR surface density')
                    sfr_cor, A_Ha, e_sfr_cor, e_A_Ha = sfr_ha(
                                tab0[prfx+'Halpha'+ext], 
                                flux_hb=tab0[prfx+'Hbeta'+ext], 
                                e_flux_ha=tab0['e_'+prfx+'Halpha'+ext],
                                e_flux_hb=tab0['e_'+prfx+'Hbeta'+ext], 
                                imf='salpeter', pixsca=pixsca,
                                name=prfx+'sigsfr_corr'+ext)
                    # For negative extinction we assume A=0
                    sfr_cor[A_Ha < ahalo]   = sfr0[A_Ha < ahalo]
                    e_sfr_cor[A_Ha < ahalo] = e_sfr0[A_Ha < ahalo]
                    # For high extinction we blank the value
                    sfr_cor[A_Ha > ahahi]   = np.nan
                    e_sfr_cor[A_Ha > ahahi] = np.nan
                    tab0.add_columns([e_sfr0, sfr_cor, e_sfr_cor, A_Ha, e_A_Ha])
                else:
                    sfr_cor, A_Ha, = sfr_ha(
                                tab0[prfx+'Halpha'+ext], 
                                flux_hb=tab0[prfx+'Hbeta'+ext], 
                                imf='salpeter', pixsca=pixsca,
                                name=prfx+'sigsfr_corr'+ext)
                    # For negative extinction we assume A=0
                    sfr_cor[A_Ha < ahalo]   = sfr0[A_Ha < ahalo]
                    # For high extinction we blank the value
                    sfr_cor[A_Ha > ahahi]   = np.nan
                    tab0.add_columns([sfr_cor, A_Ha])

                # Halpha extinction and SFR after smoothing and clipping
                A_Ha_smo = Column(get_AHa(tab0[prfx+'Halpha_sm'+str(nsm)+ext], 
                            tab0[prfx+'Hbeta_sm'+str(nsm)+ext], np.log10), 
                            name=prfx+'AHa_smooth'+str(nsm)+ext, dtype='f4', unit='mag',
                            description='Ha extinction after {} pix smooth'.format(str(nsm)))
                sfr_smo = Column(sfr0 * 10**(0.4*A_Ha_smo),
                            name=prfx+'sigsfr_adopt'+ext, dtype='f4', unit=sfr0.unit,
                            description='smooth+clip BD corrected SFR surface density')
                # For negative extinction we assume A=0
                sfr_smo[A_Ha_smo < ahalo] = sfr0[A_Ha_smo < ahalo]
                # For high extinction we blank the value
                sfr_smo[A_Ha_smo > ahahi] = np.nan
                tab0.add_columns([A_Ha_smo, sfr_smo])

                # BPT requires flux_elines since EW(Ha) is part of classification
                if prod == 'flux_elines':
                    if prob:
                        BPT0, BPT0sf, p_BPT0 = bpt_type(tab0, ext=ext, name='BPT'+ext, 
                                                    prob=prob)
                        tab0.add_columns([BPT0, p_BPT0, BPT0sf])
                    else:
                        BPT0, BPT0sf = bpt_type(tab0, ext=ext, name='BPT'+ext, 
                                                    prob=prob)
                        tab0.add_columns([BPT0, BPT0sf])
                    #
                    zoh0, zoherr0 = ZOH_M13(tab0, ext=ext, name='ZOH'+ext, err=True)
                    tab0.add_columns([zoh0, zoherr0])
                    zoh2, zoherr2 = ZOH_M13(tab0, ext=ext, name='ZOH_N2'+ext, 
                                            method='n2', err=True)
                    tab0.add_columns([zoh2, zoherr2])
                    zoh_pp04 = ZOH_M13(tab0, ext=ext, name='ZOH_PP04'+ext, 
                                            method='o3n2_pp04', err=False)
                    tab0.add_column(zoh_pp04)
                    Zprime = 10**(zoh0 - 8.69)
                    # Scaling for alphaCO from Schinnerer & Leroy 2024
                    try:
                        alpha, f_term, g_term, rco = predict_alphaCO_SL24(
                               Zprime=Zprime, Sigma_star=star0, return_all_terms=True)
                        alphsca_SL24 = Column(f_term * g_term,
                                              name='alphsca_SL24'+ext, dtype='f4',
                                              description='alphaCO scaling factor from SL24')
                        tab0.add_column(alphsca_SL24)
                    except NameError:
                        pass
                    # Scaling for alphaCO from Bolatto, Wolfire, Leroy 2013
                    try:
                        if append==True:
                            cotbl  = Table.read(outfile, path='comom_dil')
                            comom0 = cotbl[cotbl['Name']==gname]['mom0_12']
                            alpha2 = predict_alphaCO10_B13(Zprime=Zprime,
                                            WCO10kpc=comom0, Sigmaelsekpc=star0+9)
                            alphsca_B13 = Column(alpha2.value/4.3,
                                                 name='alphsca_B13'+ext, dtype='f4',
                                                 description='alphaCO scaling factor from B13')
                            tab0.add_column(alphsca_B13)
                    except NameError:
                        pass

            elif prod == 'SFH':
                # For star formation history also calculate mass fractions
                # Multiply the luminosity fraction by M/L ratio and re-normalize
                lumcols = Table(tab0.columns[9:nlumcols+9])
                df_lum = lumcols.to_pandas()
                df_mass = df_lum.multiply(models.mass_to_light, axis='columns')
                df_norm = df_mass.divide(df_mass.sum(axis=1), axis='index')
                df_norm.columns = [x.replace('lum','mass') for x in list(df_norm.columns)]
                # Add aggregated mass fraction columns to table
                agecols = [s.split('_')[2] for s in df_norm.columns.values]
                metcols = [s.split('_')[4] for s in df_norm.columns.values]
                df_age = df_norm.groupby(agecols, sort=False, axis=1).sum(min_count=1)
                df_age = df_age.reindex(sorted(df_age.columns, key=float), axis=1)
                # Total the mass fractions < 33 Myr for later SFR calculation
                sublist = (df_age.columns.values.astype(float) < 0.035)
                df_age.columns = ['massfrac_age_'+x+ext for x in list(df_age.columns)]
                f_young = np.array(df_age[df_age.columns.values[sublist]].sum(axis=1, 
                                 min_count=1).astype(np.float32))
                df_met = df_norm.groupby(metcols, axis=1).sum(min_count=1)
                df_met.columns = ['massfrac_met_'+x+ext for x in list(df_met.columns)]
                naggcols = len(df_age.columns) + len(df_met.columns)
                print('Number of aggregated columns:', naggcols)
                t_mass_age = Table.from_pandas(df_age.astype(np.float32))
                t_mass_met = Table.from_pandas(df_met.astype(np.float32))
                indexcols  = Table(tab0.columns[:9])
                lumaggcols = Table(tab0.columns[nlumcols+9:nlumcols+naggcols+9])
                if has_errors:
                    erraggcols = Table(tab0.columns[2*nlumcols+naggcols+9:])
                    tab0 = hstack([indexcols, lumaggcols, erraggcols,
                               t_mass_age.filled(np.nan), 
                               t_mass_met.filled(np.nan)], join_type='exact')
                else:
                    tab0 = hstack([indexcols, lumaggcols,
                               t_mass_age.filled(np.nan), 
                               t_mass_met.filled(np.nan)], join_type='exact')
                tab0.add_column(f_young, name='f_young')
                tab0['f_young'].description='total mass fraction < 33 Myr'
                for i_col in range(naggcols):
                    newname=lumaggcols.columns[i_col].name.replace('lum','mass')
                    newdesc=lumaggcols.columns[i_col].description.replace('Luminosity','Mass')
                    tab0[newname].description = newdesc
                    tab0[newname].unit = 'fraction'

            elif prod == 'SSP':
                # For stellar surface density we need distance
                star0 = stmass_pc2(tab0['mass_ssp'+ext], dz=tab0['cont_dezon'+ext],
                                dist=disttbl.loc[gname][coln_dmpc], pixsca=pixsca,
                                name='sigstar'+ext)
                avstar0 = stmass_pc2(tab0['mass_Avcor_ssp'+ext], dz=tab0['cont_dezon'+ext],
                                dist=disttbl.loc[gname][coln_dmpc], pixsca=pixsca,
                                name='sigstar_Avcor'+ext)
                avstar0.description += ' dust corrected'
                ferr0 = Column(abs(tab0['e_medflx_ssp'+ext]/tab0['medflx_ssp'+ext]), 
                               name='fe_medflx'+ext, dtype='f4', unit='fraction',
                               description='fractional error in continuum flux')
                tab0.add_columns([star0, avstar0, ferr0])
                # Add the SSP-based SFR if SFH was run
                try:
                    ssp_sfr = Column(f_young * star0 / (0.033*u.Gyr),
                                name='sigsfr_ssp'+ext, dtype='f4',
                                description='Sigma_SFR from < 33 Myr SSP')
                    avssp_sfr = Column(f_young * avstar0 / (0.033*u.Gyr),
                                name='sigsfr_Avcor_ssp'+ext, dtype='f4',
                                description='Sigma_SFR Av-corrected from < 33 Myr SSP')
                    tab0.add_columns([ssp_sfr, avssp_sfr])
                except NameError:
                    pass
            if i_gal == 0:
                tablist.append(tab0)
            else:
                tablist[i_prod] = vstack([tablist[i_prod],tab0], join_type='exact')

        # END loop over products
    # END loop over galaxies

    for i_prod, prod in enumerate(prodtype):
        tablist[i_prod].meta['date'] = datetime.today().strftime('%Y-%m-%d')
        if debug:
            print(prod, tablist[i_prod].colnames)
            print('There are',len(tablist[i_prod]),'rows in merged table')
        tablist[i_prod].write(outfile, path=prod+ext, overwrite=overwrite, 
                              append=True, serialize_meta=True, compression=True)

# Below code was for when the products were the outer loop and galaxies the inner
#     if len(tlist) > 0:
#         t_merge = vstack(tlist)
#     t_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
#     if debug:
#         print(t_merge.colnames)
#         print('There are',len(t_merge),'rows in merged table')
# 
#         if prod == prodtype[0]:
#             t_merge.write(outfile, path=prod+ext, overwrite=overwrite, 
#                     append=append, serialize_meta=True, compression=True)
#         else:
#     t_merge.write(outfile, path=prod+ext, overwrite=overwrite, 
#             append=True, serialize_meta=True, compression=True)

    return


if __name__ == "__main__":
    # NGC4047, CALIFA only
    do_pipe3d(outfile='NGC4047_allpix.pipe3d.hdf5', append=False, allpix=True,
              fitsdir='fits_natv_carma', packed=False, gallist=['NGC4047'])
    # NGC4047, append to 2d_smo7
    do_pipe3d(outfile='../img_comom/NGC4047.2d_smo7.hdf5', append=True, 
              comomdir='../img_comom/fitsdata', ext='_sm', nsm=3, 
              cotempl='GNAME.co.smo7_dil.snrpk.fits.gz', packed=False,
              gallist=['NGC4047'], fitsdir='fits_smo7_carma')
    do_pipe3d(outfile='../img_comom/NGC4047_hex.2d_smo7.hdf5', append=True, 
              comomdir='../img_comom/fitsdata', ext='_sm', nsm=3, 
              cotempl='GNAME.co.smo7_dil.snrpk.fits.gz', packed=False,
              gallist=['NGC4047'], fitsdir='fits_smo7_carma', hexgrid=True)

    # All EDGE125 galaxies, CALIFA only
    carma125 = [os.path.basename(file).split('.')[0] for file in 
                sorted(glob.glob('fits_smo7_carma/[A-Z]*.SSP.cube.fits.gz'))]
    print(carma125)
    do_pipe3d(gallist=carma125, outfile='edge_carma_allpix.pipe3d.hdf5', 
              fitsdir='fits_natv_carma', packed=False, append=False, allpix=True)
    # All EDGE125 galaxies, append to 2d_smo7
    do_pipe3d(outfile='../img_comom/edge_carma.2d_smo7.hdf5', append=True,
              comomdir='../img_comom/fitsdata', ext='_sm', nsm=3,
              cotempl='GNAME.co.smo7_dil.snrpk.fits.gz', packed=False,
              gallist=carma125, fitsdir='fits_smo7_carma')

    # ACA-60 galaxies, 12" resolution
    aca60   = [os.path.basename(file).split('.')[0] for file in 
               sorted(glob.glob('fits_smo12_aca/[A-Z]*.SSP.cube.fits.gz'))]
    do_pipe3d(gallist=aca60, outfile='edge_aca_allpix.pipe3d.hdf5', 
              fitsdir='fits_natv_aca', packed=False, append=False, allpix=True)
    do_pipe3d(gallist=aca60, outfile='../img_comom/edge_aca.2d_smo12.hdf5', 
              cotempl='GNAME.co21.smo12_dil.snrpk.fits.gz', comomdir='../img_comom/aca12', 
              fitsdir='fits_smo12_aca', nsm=4, ext='_sm', packed=False, append=True)

    # ALMaQUEST galaxies, native resolution, append to 2d_preregrid
    aquest = [os.path.basename(file).split('.')[0].split('-',1)[1] for file in sorted(
              glob.glob('fits_natv_aq/*.Pipe3D.cube.fits.gz'))]
    do_pipe3d(gallist=aquest, outfile='../img_comom/almaquest.2d_preregrid.hdf5',
              fitsdir='fits_natv_aq', p3dstruct='manga', packed=True,
              comomdir='../img_comom/aquest_comom_fits/', 
              cotempl='manga_GNAME.co.preregrid_dil.snrpk.fits.gz',
              ssptable='/Users/tonywong/Work/projects/MaNGA/mastar_237.fits', nsm=3,
              ortpar='/Users/tonywong/Work/projects/MaNGA/MaNGA_props_pipe3d.csv', 
              distpar='/Users/tonywong/Work/projects/MaNGA/MaNGA_props_pipe3d.csv', 
              coln_ra='objra', coln_dc='objdec', coln_pa='nsa_sersic_phi',
              coln_inc='nsa_inclination', coln_dmpc='nsa_z_dMpc', append=True)

    # All of CALIFA DR3
    # Select the 646 galaxies in CALIFA DR3 V500
    dr3 = Table.read('../dat_glob/califa/build/QCflags_std_V500_DR3.csv',
                     format='ascii.no_header')
    dr3.keep_columns(['col2','col18'])
    dr3.rename_columns(['col2','col18'],['Name','in_dr3'])
    dr3only = dr3[ dr3['in_dr3']==1 ]
    califa = sorted(list(dr3only['Name']))
    do_pipe3d(outfile='dr3_allpix.pipe3d.hdf5', gallist=califa, fitsdir='fits_califa',
              p3dstruct='califa', packed=True, ortpar='edge_leda.csv', 
              distpar='edge_califa.csv', prob=False, append=False, allpix=True)

    # All of MaNGA
    gallist = [os.path.basename(file).split('.')[0].split('-',1)[1] for file in sorted(
                glob.glob('fits_manga/*.Pipe3D.cube.fits.gz'))]
    do_pipe3d(gallist=gallist, outfile='manga_allpix.pipe3d.hdf5', allpix=True,
              fitsdir='fits_manga', p3dstruct='manga', packed=True,
              ssptable='/Users/tonywong/Work/projects/MaNGA/mastar_237.fits', nsm=3,
              ortpar='/Users/tonywong/Work/projects/MaNGA/MaNGA_props_pipe3d.csv', 
              distpar='/Users/tonywong/Work/projects/MaNGA/MaNGA_props_pipe3d.csv', 
              coln_ra='objra', coln_dc='objdec', coln_pa='nsa_sersic_phi',
              coln_inc='nsa_inclination', coln_dmpc='nsa_z_dMpc', prob=False, append=False)
    