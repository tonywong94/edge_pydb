#!/usr/bin/env python

# Combine the CO cubes into binary tables.

from datetime import datetime
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack
from edge_pydb import EdgeTable
from edge_pydb.fitsextract import fitsextract

def do_cocube(outname='NGC4047', gallist=['NGC4047'], seq='smo7', lines=['12','13'], 
             linelbl=['co','13co'], colm=['data3d','rms3d','dilmsk3d','smomsk3d'], 
             unit=['K','K','',''], colmlbl=['msk.K','_dil.ecube','_dil.mask','_smo.mask'], 
             hexgrid=False, allpix=False, fitsdir='fitsdata', ortpar='edge_leda.csv'):
    """
    Extract 3D molecular line data into an HDF5 database.  This script assumes
    standardized naming conventions, for example:
        UGC10710.co.smo7_dil.mask.fits.gz = 
            ${galaxy}.${linelbl}.${seq}${colmlbl}.fits.gz
    The data cube is assumed to be uncompressed, the errors and masks should be
    gzip compressed.

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
    colm : list of str
        How different file types will be identified in the database
    unit : list of str
        The brightness unit for each member of colm
    colmlbl : list of str
        How different file types are identified in the FITS file names
    hexgrid : boolean
        True to sample on a hexagonal grid (experimental)
    allpix : boolean
        True to dump every pixel, otherwise every 3rd pixel in x and y is used.
    fitsdir : str
        Path to the directory where FITS files reside
    ortpar : filename
        Name of the EdgeTable which has LEDA orientation parameters for the sample
    """
    if allpix:
        stride = [1,1,1]
    else:
        stride = [3,3,1]

    # Get the orientation parameters from LEDA
    orttbl = EdgeTable(ortpar)
    orttbl.add_index('Name') 

    tablelist=[]
    for gal in gallist:
        file0 = os.path.join(fitsdir, gal+'.'+linelbl[0]+'.'+seq+colmlbl[0]+'.fits')
        if not os.path.exists(file0):
            print('####### Cannot find',file0)
            continue
        for i_line, line in enumerate(lines):
            for i_col in range(len(colm)):
                # --- Read the first image (main cube data)
                if i_line == 0 and i_col == 0:
                    print('Reading',file0)
                    galtab = fitsextract(file0, bunit=unit[i_col], 
                            col_lbl=colm[i_col]+'_'+line,
                            keepnan=True, stride=stride,
                            ra_gc=15*orttbl.loc[gal]['ledaRA'],
                            dec_gc=orttbl.loc[gal]['ledaDE'],
                            pa=orttbl.loc[gal]['ledaPA'],
                            inc=orttbl.loc[gal]['ledaAxIncl'],
                            ortlabel='LEDA', first=True,
                            use_hexgrid=hexgrid)
                    gname = Column([np.string_(gal)]*len(galtab), name='Name', description='Galaxy Name')
                    galtab.add_column(gname, index=0)
                    print(galtab[20:50])
                # --- Read the subsequent images (assumed to be gzipped if i_col>0)
                else:
                    if i_col == 0:
                        getfile = os.path.join(fitsdir,
                            gal+'.'+linelbl[i_line]+'.'+seq+colmlbl[i_col]+'.fits')
                    else:
                        getfile = os.path.join(fitsdir,
                            gal+'.'+linelbl[i_line]+'.'+seq+colmlbl[i_col]+'.fits.gz')
                    if os.path.exists(getfile):
                        print('Reading',getfile)
                        addtb = fitsextract(getfile, bunit=unit[i_col], 
                                    col_lbl=colm[i_col]+'_'+line, 
                                    keepnan=True, stride=stride, 
                                    use_hexgrid=hexgrid)
                        jointb = join(galtab, addtb, keys=['ix','iy','iz'])
                        galtab = jointb
                    else:
                        print('####### Cannot find',getfile)
                        newcol = Column(data=[np.nan]*len(galtab), 
                                        name=colm[i_col]+'_'+line, 
                                        unit=unit[i_col], dtype='f4')
                        galtab.add_column(newcol)
        tablelist.append(galtab)

    if len(tablelist) > 0:
        t_merge = vstack(tablelist)
        print(t_merge[20:50])
        for i_line, line in enumerate(lines):
            t_merge[colm[0]+'_'+line].description = linelbl[i_line]+' brightness temperature in cube'
            t_merge[colm[1]+'_'+line].description = linelbl[i_line]+' estimated 1-sigma channel noise'
            t_merge[colm[2]+'_'+line].description = linelbl[i_line]+' mask value for dilated mask'
            t_merge[colm[3]+'_'+line].description = linelbl[i_line]+' mask value for smoothed mask'
        t_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')
        t_merge.meta['comments'] = 'Sampled CO and 13CO data cubes from EDGE'
        t_merge.write(outname+'.cocube_'+seq+'.hdf5', path='data', overwrite=True, 
                serialize_meta=True, compression=True)
    return


if __name__ == "__main__":
    # NGC4047 only
    do_cocube()
    # All EDGE125 galaxies
    gallist = [os.path.basename(file).split('.')[0] for file in 
               sorted(glob.glob('fitsdata/*.co.smo7msk.K.fits'))]
    do_cocube(gallist=gallist, outname='edge_carma')
    # EDGE125 hexgrid - not yet working
    # do_cocube(gallist=gallist, outname='edge_carma_hex', hexgrid=True)
    # EDGE125 allpix data
    do_cocube(gallist=gallist, outname='edge_carma_allpix', allpix=True)
