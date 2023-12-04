#!/usr/bin/env python

# Combine the CO moment maps into binary tables.

from datetime import datetime
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack
from edge_pydb import EdgeTable
from edge_pydb.conversion import msd_co
from edge_pydb.fitsextract import fitsextract


def do_comom(outname='NGC4047', gallist=['NGC4047'], seq='smo7', lines=['12','13'],
             linelbl=['co','13co'], msktyp=['str', 'dil', 'smo'], alphaco=4.35, 
             hexgrid=False, allpix=False, fitsdir='fitsdata', ortpar='edge_leda.csv',
             append=False, overwrite=True):
    """
    Extract 2D molecular line data into an HDF5 database.  This script assumes
    standardized naming conventions, for example:
        UGC10710.co.smo7_smo.emom2.fits.gz = 
            ${galaxy}.${linelbl}.${seq}_${msktyp}.${ftype}.fits.gz
    The possible values for ${ftype} need to be defined within 'dotypes'.

    Parameters
    ----------
    outname : str
        Prefix of the output filename
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
    alphaco : float
        CO to H2 conversion factor, in Msol/pc2/(K km/s).  Default=4.3
    hexgrid : boolean
        True to sample on a hexagonal grid (experimental)
    allpix : boolean
        True to dump every pixel, otherwise every 3rd pixel in x and y is used.
    fitsdir : str
        Path to the directory where FITS files reside
    ortpar : filename
        Name of the EdgeTable which has LEDA orientation parameters for the sample
    append : boolean
        True to append to an existing file.  Default is False (create/overwrite).
    overwrite : boolean
        True to overwrite existing tables.  Default is True (replace the table
        but do not delete other tables in the file).
    """
    if allpix:
        stride = [1,1,1]
    else:
        stride = [3,3,1]

    # Get the orientation parameters from LEDA
    orttbl = EdgeTable(ortpar)
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
            dotypes = ['snrpk', 'mom0',   'e_mom0', 'mom1', 'e_mom1', 'mom2', 'e_mom2']
            unit    = ['',      'K km/s', 'K km/s', 'km/s', 'km/s',   'km/s', 'km/s']
        for gal in gallist:
            # snrpk.fits is produced for dilated mask only
            if msk == 'smo':
                file0 = os.path.join(fitsdir,
                        gal+'.'+linelbl[0]+'.'+seq+'_dil.snrpk.fits.gz')
            else:
                file0 = os.path.join(fitsdir,
                        gal+'.'+linelbl[0]+'.'+seq+'_'+msk+'.'+dotypes[0]+'.fits.gz')
            if not os.path.exists(file0):
                continue
            adopt_incl = orttbl.loc[gal]['ledaAxIncl']
            print('Adopted inclination is {} deg'.format(adopt_incl))
            for i_line, line in enumerate(lines):
                for i_mtype, mtype in enumerate(dotypes):
                    # --- Read the first image (should be snrpk or mom0)
                    if i_line == 0 and i_mtype == 0:
                        print('Reading',file0)
                        galtab = fitsextract(file0, bunit=unit[0], 
                                col_lbl=dotypes[0]+'_'+line,
                                keepnan=True, stride=stride,
                                ra_gc=15*orttbl.loc[gal]['ledaRA'],
                                dec_gc=orttbl.loc[gal]['ledaDE'],
                                pa=orttbl.loc[gal]['ledaPA'],
                                inc=adopt_incl,
                                ortlabel='LEDA', first=True,
                                use_hexgrid=hexgrid)
                        gname = Column([np.string_(gal)]*len(galtab), name='Name', description='Galaxy Name')
                        galtab.add_column(gname, index=0)
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
                            addtb = fitsextract(getfile, bunit=unit[i_mtype], col_lbl=mtype+'_'+line, 
                                            keepnan=True, stride=stride, use_hexgrid=hexgrid)
                            jointb = join(galtab, addtb, keys=['ix','iy'])
                            galtab = jointb
                        else:
                            newcol = Column(data=[np.nan]*len(galtab), name=mtype+'_'+line, 
                                            unit=unit[i_mtype], dtype='f4')
                            galtab.add_column(newcol)
                # Add the H2 column density, with and without deprojection
                if line == '12':
                    cosi = Column([np.cos(np.radians(adopt_incl))]*len(galtab), name='cosi', 
                            description='factor to deproject to face-on using ledaAxIncl', dtype='f4')
                    sigmol = msd_co(galtab['mom0_12'], name='sigmol', alphaco=alphaco)
                    e_sigmol = msd_co(galtab['e_mom0_12'], name='e_sigmol', alphaco=alphaco)
                    galtab.add_columns([sigmol, e_sigmol, cosi])
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
            t_merge.write(outname+'.2d_'+seq+'.hdf5', path='comom_'+msk, append=append,
                    overwrite=overwrite, serialize_meta=True, compression=True)
        else:
            t_merge.write(outname+'.2d_'+seq+'.hdf5', path='comom_'+msk, append=True, 
                    overwrite=overwrite, serialize_meta=True, compression=True)
    return


if __name__ == "__main__":
    # NGC4047 only
    do_comom(outname='NGC4047')
    do_comom(hexgrid=True, outname='NGC4047_hex')

    # All EDGE125 galaxies, 7" resolution
    gallist = [os.path.basename(file).split('.')[0] for file in 
               sorted(glob.glob('fitsdata/*.co.smo7_dil.snrpk.fits.gz'))]
    do_comom(gallist=gallist, outname='edge_carma')
    do_comom(gallist=gallist, outname='edge_carma_allpix', allpix=True)
    # EDGE125 hexgrid
    do_comom(gallist=gallist, outname='edge_carma_hex', hexgrid=True)
    
    # ACA galaxies, 12" resolution, using alphaco=6.6 instead of 4.3
    # R21 = 0.65 (Leroy+22)
#     gallist = [os.path.basename(file).split('.')[0] for file in 
#                sorted(glob.glob('aca12/*_dil.snrpk.fits.gz'))]
#     do_comom(gallist=gallist, outname='edge_aca', seq='smo12', lines=['12'],
#              linelbl=['co21'], msktyp=['str', 'dil'], alphaco=6.6, fitsdir='aca12')
#     do_comom(gallist=gallist, outname='edge_aca_allpix', seq='smo12', lines=['12'],
#              linelbl=['co21'], msktyp=['str', 'dil'], alphaco=6.6, fitsdir='aca12', 
#              allpix=True)

