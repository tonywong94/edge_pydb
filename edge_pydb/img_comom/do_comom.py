#!/usr/bin/env python

# Combine the 7" CO moment maps into binary tables.

from datetime import datetime
import glob
import os
import numpy as np
from astropy.table import Table, Column, join, vstack
from edge_pydb import EdgeTable
from edge_pydb.conversion import msd_co
from edge_pydb.fitsextract import fitsextract

seq = 'smo7'
#seq = 'de20'
msktyp = ['str', 'dil', 'smo']
lines  = ['12', '13']

# Get the orientation parameters from LEDA
ort = EdgeTable('edge_leda.csv', cols=['Name', 'ledaRA', 'ledaDE', 'ledaPA', 'ledaAxIncl'])
ort.add_index('Name')

for imsk, msk in enumerate(msktyp):
    gallist = [os.path.basename(file).split('.')[0] for file in 
                sorted(glob.glob('fitsdata/*.co.'+seq+'_dil.snrpk.fits.gz'))] 
    tablelist=[]
    if msk == 'str':
        dotypes = ['mom0', 'emom0']
        unit    = ['K km/s', 'K km/s']
    if msk == 'dil':
        dotypes = ['snrpk', 'mom0', 'emom0', 'mom1', 'emom1', 'mom2', 'emom2']
        unit    = ['', 'K km/s', 'K km/s', 'km/s', 'km/s', 'km/s', 'km/s']
    if msk == 'smo':
        dotypes = ['mom0', 'emom0', 'mom1', 'emom1', 'mom2', 'emom2']
        unit    = ['K km/s', 'K km/s', 'km/s', 'km/s', 'km/s', 'km/s']
    for gal in gallist:
        file0 = 'fitsdata/'+gal+'.co.'+seq+'_'+msk+'.'+dotypes[0]+'.fits.gz'
        print(file0)
        if not os.path.exists(file0):
            continue
        adopt_incl = ort.loc[gal]['ledaAxIncl']
        print('Adopted inclination is {} deg'.format(adopt_incl))
        # Read the first image
        for iline, line in enumerate(lines):
            for j, type in enumerate(dotypes):
                if iline == 0 and j == 0:
                    print('Reading',file0)
                    galtab = fitsextract(file0, bunit=unit[0], 
                            col_lbl=dotypes[0]+'_'+line,
                            keepnan=True, stride=[3,3,1],
                            ra_gc=15*ort.loc[gal]['ledaRA'],
                            dec_gc=ort.loc[gal]['ledaDE'],
                            pa=ort.loc[gal]['ledaPA'],
                            inc=adopt_incl,
                            ortlabel='LEDA', first=True)
                    gname = Column([np.string_(gal)]*len(galtab), name='Name', description='Galaxy Name')
                    galtab.add_column(gname, index=0)
                    print(galtab[20:50])
                else:
                    if line == '12' or msk == 'str' or type == 'snrpk':
                        getfile = 'fitsdata/'+gal+'.co.'+seq+'_'+msk+'.'+type+'.fits.gz'
                    else:
                        getfile = 'fitsdata/'+gal+'.13co.'+seq+'_mk12_'+msk+'.'+type+'.fits.gz'
                    if os.path.exists(getfile):
                        print('Reading',getfile)
                        addtb = fitsextract(getfile, bunit=unit[j], col_lbl=type+'_'+line, 
                                        keepnan=True, stride=[3,3,1])
                        jointb = join(galtab, addtb, keys=['ix','iy'])
                        galtab = jointb
                    else:
                        newcol = Column(data=[np.nan]*len(galtab), name=type+'_'+line, 
                                        unit=unit[j], dtype='f4')
                        galtab.add_column(newcol)
            # Add the H2 column density, with and without deprojection
            if line == '12':
                sigmol = msd_co(galtab['mom0_12'], name='sigmol')
                e_sigmol = msd_co(galtab['emom0_12'], name='e_sigmol')
                sigmol_fo = msd_co(galtab['mom0_12']*np.cos(np.radians(adopt_incl)), name='sigmol_fo')
                e_sigmol_fo = msd_co(galtab['emom0_12']*np.cos(np.radians(adopt_incl)), name='e_sigmol_fo')
                galtab.add_columns([sigmol, e_sigmol, sigmol_fo, e_sigmol_fo])
        tablelist.append(galtab)

    if len(tablelist) > 0:
        t_merge = vstack(tablelist)
        for line in lines:
            if 'snrpk' in dotypes:
                t_merge['snrpk_'+line].description = line+'CO peak signal to noise ratio'
            t_merge['mom0_'+line].description = line+'CO integrated intensity using {} mask'.format(msk)
            t_merge['emom0_'+line].description = line+'CO error in mom0 assuming {} mask'.format(msk)
            if msk == 'dil':
                t_merge['mom1_'+line].description = line+'CO intensity wgtd mean velocity using {} mask'.format(msk)
                t_merge['emom1_'+line].description = line+'CO error in mom1 assuming {} mask'.format(msk)
                t_merge['mom2_'+line].description = line+'CO intensity wgtd vel disp using {} mask'.format(msk)
                t_merge['emom2_'+line].description = line+'CO error in mom2 assuming {} mask'.format(msk)
            if line == '12':
                t_merge['sigmol'].description = 'apparent H2+He surf density not deprojected'
                t_merge['e_sigmol'].description = 'error in sigmol not deprojected'
                t_merge['sigmol_fo'].description = 'H2+He surf density deprojected to face-on using ledaAxIncl'
                t_merge['e_sigmol_fo'].description = 'error in sigmol deprojected to face-on'
        t_merge.meta['date'] = datetime.today().strftime('%Y-%m-%d')   
        print(t_merge[20:50])

    if (len(gallist) > 1):
        outname = 'edge'
    else:
        outname = gal
    if imsk == 0:
        t_merge.write(outname+'.comom_'+seq+'.hdf5', path=msk, overwrite=True, 
                serialize_meta=True, compression=True)
    else:
        t_merge.write(outname+'.comom_'+seq+'.hdf5', path=msk, append=True, 
                serialize_meta=True, compression=True)
