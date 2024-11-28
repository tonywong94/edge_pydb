#!/usr/bin/env python

# Combine the CO moment maps into binary tables.

from datetime import datetime
import glob
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column, join, vstack
from edge_pydb import EdgeTable
from edge_pydb.conversion import msd_co
from edge_pydb.fitsextract import fitsextract
from reproject import reproject_interp
from CO_conversion_factor.alphaCO import predict_alphaCO10_B13, predict_alphaCO_SL24

def do_comom(outfile='NGC4047.2d_smo7.hdf5', gallist=['NGC4047'], seq='smo7', 
             lines=['12','13'], linelbl=['co','13co'], msktyp=['str', 'dil', 'smo'], 
             alphaco=4.3, hexgrid=False, allpix=False, fitsdir='fitsdata', 
             ortpar='edge_leda.csv', ortlabel='LEDA', coln_ra='ledaRA', 
             coln_dc='ledaDE', coln_pa='ledaPA', coln_inc='ledaAxIncl', deproj='inc', 
             p3d_dir=None, interp_order=1, p3dtempl='flux_elines.GNAME.cube.fits.gz', 
             zoh_col='ZOH_PP04_cobm', append=True, overwrite=True, manganame=False):
    """
    Extract 2D molecular line data into an HDF5 database.  This script assumes
    standardized naming conventions, for example:
        UGC10710.co.smo7_smo.emom2.fits.gz = 
            ${galaxy}.${linelbl}.${seq}_${msktyp}.${ftype}.fits.gz
    The possible values for ${ftype} need to be defined within 'dotypes'.
    For the 13CO line, masked moments using the corresponding 12CO mask are used.
    
    The <p3d_dir> parameter allows the maps to be regridded to match
    the Pipe3D data, otherwise the CO astrometric grid is used.

    Parameters
    ----------
    outfile : str
        Name of the output HDF5 filename.  Appended to if it exists.
    gallist : list of str
        List of galaxy names
    seq : str
        Identifier, generally to indicate smoothing resolution
    lines : list of str
        How different lines will be identified in the database
    linelbl : list of str
        How different lines are identified in the FITS file names
    msktyp : list of str
        The types of masks to include.  Each mask is a separate path in the HDF5 file.
        Only the following choices are supported:
            'str' : no masking, for mom-0 only
            'dil' : dilated mask, for all moments and peak SNR image
            'smo' : smoothed and dilated mask, for all moments
    alphaco : float
        CO to H2 conversion factor, in Msol/pc2/(K km/s).  Default=4.3
    hexgrid : boolean
        True to sample on a hexagonal grid (experimental)
    allpix : boolean
        True to dump every pixel, otherwise every 3rd pixel in x and y is used.
    fitsdir : str
        Path to the directory where the CO FITS files reside
    ortpar : filename
        Name of the EdgeTable which has orientation parameters for the sample, 
        including center position, inclination, and position angle.  This can also 
        be given as the path to a regular astropy-compatible Table.
    ortlabel : str
        String for labeling the origin of the orientation parameters.  Default is 'LEDA'.
    coln_ra : str
        Name of the RA column in 'ortpar' to use.  Default is 'ledaRA'.
    coln_dc : str
        Name of the DEC column in 'ortpar' to use.  Default is 'ledaDE'.
    coln_pa : str
        Name of the PA column in 'ortpar' to use.  Default is 'ledaPA'.
    coln_inc : str
        Name of the INC column in 'ortpar' to use.  Default is 'ledaAxIncl'.
    deproj : str
        Whether to interpret the inc column as an inclination in degrees (deproj='inc') 
        or an axis ratio b/a (deproj='axrat') when determining cosi.  Default is 'inc'.
    p3d_dir : str
        Path to the directory where Pipe3D files reside.  Default is None (do not
        regrid CO maps to match Pipe3D).
    interp_order : int
        Interpolation order for reproject.  Default = 1 (linear).
    p3dtempl : str
        File name of regridding template, where GNAME is replaced by the galaxy name.
    zoh_col : str
        Column name to use for metallicity from flux_elines table.
    append : boolean
        True (default) to append to an existing file.  Use False to create/overwrite.
        False also suppresses calculation of the metallicity-dependent conversion factor.
    overwrite : boolean
        True to overwrite existing tables.  Default is True (replace the table
        but do not delete other tables in the file).
    manganame : boolean
        For economy MaNGA galaxies are labeled as (e.g.) '8952-6104' in the Pipe3D
        table whereas the ALMA file names use a longer string like 'manga_8952_6104'. 
        If True, this translation is made to allow the names to be matched.
    """
    if allpix:
        stride = [1,1,1]
    else:
        stride = [3,3,1]

    # Get the orientation parameters from LEDA
    try:
        orttbl = EdgeTable(ortpar)
    except:
        orttbl = Table.read(ortpar, format='ascii.ecsv')
    orttbl.add_index('Name') 

    for i_msk, msk in enumerate(msktyp):
        tablelist=[]
        if msk == 'str':
            dotypes = ['mom0',   'e_mom0']
            unit    = ['K km/s', 'K km/s']
        if msk == 'dil':
            dotypes = ['snrpk', 'mom0',   'e_mom0', 'mom1', 'e_mom1', 'mom2', 'e_mom2']
            unit    = ['',      'K km/s', 'K km/s', 'km/s', 'km/s',   'km/s', 'km/s']
        if msk == 'smo':
            dotypes = ['mom0',   'e_mom0', 'mom1', 'e_mom1', 'mom2', 'e_mom2']
            unit    = ['K km/s', 'K km/s', 'km/s', 'km/s',   'km/s', 'km/s']
        for gal in gallist:
            # snrpk.fits is produced for dilated mask only
#             if msk == 'smo':
#                 file0 = os.path.join(fitsdir,
#                         gal+'.'+linelbl[0]+'.'+seq+'_dil.snrpk.fits.gz')
#             else:
            file0 = os.path.join(fitsdir,
                    gal+'.'+linelbl[0]+'.'+seq+'_'+msk+'.'+dotypes[0]+'.fits.gz')
            if not os.path.exists(file0):
                continue
            # Regrid to the Pipe3D template
            hdul = fits.open(file0, ignore_missing_end=True)
            newhd = hdul[0].header.copy()
            if p3d_dir is not None:
                p3dfile = os.path.join(p3d_dir,p3dtempl.replace('GNAME',gal))
                if not os.path.exists(p3dfile):
                    # Fudge for almaquest file naming convention
                    gal2 = gal.replace('_','-')
                    p3dfile = os.path.join(p3d_dir,p3dtempl.replace('GNAME',gal2))
                    if not os.path.exists(p3dfile):
                        print('####### Cannot find',p3dfile)
                        continue
                p3dhd = fits.getheader(p3dfile)
                hd2d = WCS(p3dhd).celestial.to_header()
                for key in hd2d.keys():
                    newhd[key] = hd2d[key]
                # For packed files NAXIS1 and NAXIS2 may be missing
                for key in ['NAXIS1', 'NAXIS2']:
                    if key in p3dhd.keys():
                        newhd[key] = p3dhd[key]
                    else:
                        p3dhd1 = fits.getheader(p3dfile, 1)
                        newhd[key] = p3dhd1[key]
                newim = reproject_interp(hdul[0], newhd, order=interp_order,
                                         return_footprint=False)
#                 if debug:
#                     print(WCS(p3dhd))
#                     fits.writeto(file0.replace('.fits','_rg.fits'), newim, newhd, overwrite=True)
            else:
                newim = hdul[0].data
            hdul.close()
            if manganame:
                gname = gal.replace('manga_','').replace('_','-')
            else:
                gname = gal
            if deproj == 'inc':
                adopt_incl = orttbl.loc[gname][coln_inc]
                adopt_cosi = np.cos(np.radians(adopt_incl))
            elif deproj == 'axrat':
                adopt_cosi = orttbl.loc[gname][coln_inc]
                adopt_incl = np.degrees(np.arccos(adopt_cosi))
            print('Adopted inclination, axis ratio is {} deg, {}'.format(adopt_incl, adopt_cosi))
            for i_line, line in enumerate(lines):
                for i_mtype, mtype in enumerate(dotypes):
                    # --- Read the first image (should be snrpk or mom0)
                    if i_line == 0 and i_mtype == 0:
                        print('Reading',file0)
                        galtab = fitsextract(newim, header=newhd, bunit=unit[0], 
                                col_lbl=dotypes[0]+'_'+line,
                                keepnan=True, stride=stride,
                                ra_gc=orttbl.loc[gname][coln_ra],
                                dec_gc=orttbl.loc[gname][coln_dc],
                                pa=orttbl.loc[gname][coln_pa],
                                inc=adopt_incl,
                                ortlabel=ortlabel, first=True,
                                use_hexgrid=hexgrid)
                        gnamecol = Column([np.string_(gname)]*len(galtab), name='Name', description='Galaxy Name')
                        galtab.add_column(gnamecol, index=0)
                        print(galtab[20:50])
                    # --- Read the subsequent images
                    else:
                        ftype = mtype.replace('e_m','em',1)
                        if line != '13':
                            getfile = os.path.join(fitsdir,
                                gal+'.'+linelbl[i_line]+'.'+seq+'_'+msk+'.'+ftype+'.fits.gz')
                        elif msk == 'str' or mtype == 'snrpk':
                            getfile = os.path.join(fitsdir,
                                gal+'.'+linelbl[i_line]+'.'+seq+'_'+msk+'.'+ftype+'.fits.gz')
                        else:
                            getfile = os.path.join(fitsdir,
                                gal+'.'+linelbl[i_line]+'.'+seq+'_mk12_'+msk+'.'+ftype+'.fits.gz')
                        if os.path.exists(getfile):
                            print('Reading',getfile)
                            if p3d_dir is None:
                                addtb = fitsextract(getfile, bunit=unit[i_mtype], 
                                                    col_lbl=mtype+'_'+line, 
                                                    keepnan=True, stride=stride, 
                                                    use_hexgrid=hexgrid)
                            else:
                                hdul = fits.open(getfile, ignore_missing_end=True)
                                newim = reproject_interp(hdul[0], newhd, order=interp_order,
                                                         return_footprint=False)
                                addtb = fitsextract(newim, header=newhd, bunit=unit[i_mtype], 
                                                    col_lbl=mtype+'_'+line, 
                                                    keepnan=True, stride=stride, 
                                                    use_hexgrid=hexgrid)
                                hdul.close()
                            jointb = join(galtab, addtb, keys=['ix','iy'])
                            galtab = jointb
                        else:
                            newcol = Column(data=[np.nan]*len(galtab), name=mtype+'_'+line, 
                                            unit=unit[i_mtype], dtype='f4')
                            galtab.add_column(newcol)
                # Add the H2 column density, with and without deprojection
                if line == '12':
                    cosi = Column([adopt_cosi]*len(galtab), name='cosi', 
                            description='factor to deproject to face-on using {}'.format(coln_inc), dtype='f4')
                    sigmol = msd_co(galtab['mom0_12'], name='sigmol', alphaco=alphaco)
                    e_sigmol = msd_co(galtab['e_mom0_12'], name='e_sigmol', alphaco=alphaco)
                    galtab.add_columns([sigmol, e_sigmol, cosi])
                    if append and msk == 'dil':
                        # Scaling for alphaCO from Bolatto, Wolfire, Leroy 2013
                        try:
                            ssptab   = Table.read(outfile, path='SSP')
                            fluxtab  = Table.read(outfile, path='flux_elines')
                            star0    = ssptab[ssptab['Name']==gname]['sigstar']
                            zoh0     = fluxtab[fluxtab['Name']==gname][zoh_col]
                            Zprime   = 10**(zoh0 - 8.69)
                            # Bolatto+13, iterative mode, based on metallicity,
                            # kpc-scale CO brightness, and stellar surface density (+9 for HI)
                            # A minimum alpha_CO is imposed by the optically thin limit and 30 K
                            alpha2   = predict_alphaCO10_B13(Zprime=Zprime,
                                            WCO10kpc=galtab['mom0_12'], Sigmaelsekpc=star0+9)
                            alphsca_B13 = Column(alpha2.value/4.3, name='alphsca_B13', unit=None,
                                   dtype='f4', description='alphaCO scaling factor from B13')
                            galtab.add_column(alphsca_B13)
                            # Scaling for alphaCO from Schinnerer & Leroy 2024
                            # SFR term not yet included (so only valid for J=1-0)
                            alpha, f_term, g_term, rco = predict_alphaCO_SL24(
                                   Zprime=Zprime, Sigma_star=star0, return_all_terms=True)
                            alphsca_SL24 = Column(f_term * g_term, name='alphsca_SL24', unit=None,
                                   dtype='f4', description='alphaCO scaling factor from SL24')
                            galtab.add_column(alphsca_SL24)
                        except:
                            print('Paths missing from output file: SSP, flux_elines')
            tablelist.append(galtab)

        if len(tablelist) > 0:
            t_merge = vstack(tablelist)
            for i_line, line in enumerate(lines):
                if 'snrpk' in dotypes:
                    t_merge['snrpk_'+line].description = linelbl[i_line]+' peak signal to noise ratio'
                t_merge['mom0_'+line].description = linelbl[i_line]+' integrated intensity using {} mask'.format(msk)
                t_merge['e_mom0_'+line].description = linelbl[i_line]+' error in mom0 assuming {} mask'.format(msk)
                if msk != 'str':
                    t_merge['mom1_'+line].description = linelbl[i_line]+' intensity wgtd mean velocity using {} mask'.format(msk)
                    t_merge['e_mom1_'+line].description = linelbl[i_line]+' error in mom1 assuming {} mask'.format(msk)
                    t_merge['mom2_'+line].description = linelbl[i_line]+' intensity wgtd vel disp using {} mask'.format(msk)
                    t_merge['e_mom2_'+line].description = linelbl[i_line]+' error in mom2 assuming {} mask'.format(msk)
                if line == '12':
                    t_merge['sigmol'].description = 'apparent H2+He surf density not deprojected'
                    t_merge['e_sigmol'].description = 'error in sigmol not deprojected'
            t_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')   
            print(t_merge[20:50])

        if i_msk == 0:
            t_merge.write(outfile, path='comom_'+msk, append=append,
                    overwrite=overwrite, serialize_meta=True, compression=True)
        else:
            t_merge.write(outfile, path='comom_'+msk, append=True, 
                    overwrite=overwrite, serialize_meta=True, compression=True)
    return


if __name__ == "__main__":

    ## CO Astrometric Grid (run first)
    # NGC4047 only
    do_comom(outfile='NGC4047.2d_smo7.hdf5', append=False)
    do_comom(outfile='NGC4047_hex.2d_smo7.hdf5', hexgrid=True, append=False)
