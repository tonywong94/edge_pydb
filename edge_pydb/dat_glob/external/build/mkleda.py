#!/usr/bin/env python

from astropy.table import Table, vstack, join
from astropy import units as u
import numpy as np
from datetime import datetime
from astroquery.ipac.ned import Ned

import warnings
warnings.filterwarnings('ignore')

#all_list = ['edge', 'edge_aca']
all_list = ['dr3']

def phot_inc(axratio, ttype):
    q = axratio  # minor/major axis ratio
    q0 = 10**(-0.60 - 0.045*ttype)  # Bottinelli+ 83 (1983A&A...118....4B)
    # For 10 galaxies where the observed q < q_0, we replace q_0 by the smaller 
    #   of 0.135 and q to avoid an undefined value for cos i. 
    #   q_0=0.135 corresponds to ttype=6, close to the latest types in EDGE.
    q0[q<q0] = np.minimum(0.135,q[q<q0])
    inc = np.degrees(np.arccos( np.sqrt((q**2-q0**2)/(1-q0**2)) ))  # Hubble 1926ApJ....64..321H
    return inc

# Retrieve PA from NED if no LEDA value
def getNedPA(galname, defvalue=-25):
    ned_table = Ned.get_table(galname, table='diameters')
    ned_table.add_index('Frequency targeted')
    for keysel in ['r (SDSS Isophotal)','K_s (2MASS isophotal)',
                   'K_s (LGA/2MASS isophotal)']:
        try:
            idx = ned_table.loc_indices[keysel]
            print('Found',keysel,'for',galname)
            if not np.isscalar(idx):
                idx = idx[0]
            return ned_table[idx]['NED Position Angle']
        except KeyError:
            #print('Did not find',keysel,'for',galname)
            pass
    return defvalue

    
# Dict for renaming the LEDA columns
colnames = {
            'name' : 'Name',
            'al2000' : 'ledaRA',
            'de2000' : 'ledaDE',
            'type' : 'ledaMorph',
            'bar' : 'ledaBar',
            'ring' : 'ledaRing',
            'multiple': 'ledaMultiple',
            't' : 'ledaType',
            'pa' : 'ledaPA',
            'bt' : 'ledaBt',
            'it' : 'ledaIt',
            'vmaxg' : 'ledaVmaxg',
            'mfir' : 'ledaFIR',
            'v' : 'ledaVrad',
            'ag' : 'ledaA_Bgal',
            'incl' : 'ledaIncl',
            'vrot' : 'ledaVrot',
            'vvir' : 'ledaVvir',
            'modz' : 'ledaModz'
           }

tablist = []

# Get CALIFA ID numbers
pipe3d = Table.read('../../califa/build/CALIFA_3_joint_classnum_new_update.csv',
                    format='ascii.no_header')
idtab = Table(pipe3d['col1','col3'])
idtab.rename_column('col1','ID')
idtab.rename_column('col3','Name')

for list in all_list:
    print('Working on',list)
    # Process the output from HyperLEDA
    if list == 'edge_aca':
        t = Table.read(list+'_ledaout.txt',format='ascii.fixed_width',header_start=2)
        t = Table(t, masked=True, copy=False)
        t.remove_columns(['name'])
        t['objname'] = t['objname'].astype('<U15')
        t.rename_column('objname','name')
    else:
        t = Table.read(list+'_ledaout.txt',format='ascii.fixed_width',header_start=2)
        t['name'] = t['name'].astype('|S24')
        t.add_column(t['name'], index=1, name='ledaName')
        t.remove_columns(['objname'])

    for key in colnames.keys():
        t.rename_column(key, colnames[key])
    t['ledaName'].description = 'Galaxy name as known by LEDA' 

    # Some CALIFA names were unrecognized by LEDA and had to be manually substituted
    if list == 'edge':
        renames = {
                   # LEDAname : CALIFAname
                   'NGC4211A'  : 'NGC4211NED02',
                   'NGC6027E'  : 'NGC6027',
                   'PGC029708' : 'UGC05498NED01',
                   'UGC04299'  : 'IC2247'
                  }
    if list == 'alma':
        renames = {
                   'PGC000312' : 'IC1528',
                   'PGC072803' : 'NGC7783NED01',
                   'PGC002029' : 'UGC00335NED02',
                   'PGC066150' : 'UGC11680NED02',
                   'NGC7237'   : 'UGC11958',
                   'PGC070084' : 'VV488NED02'
                  }
    if list == 'edge_aca':
        renames = {
                   'PGC069310' : 'CGCG429-012',
                   'PGC000312' : 'IC1528',
                   'PGC073143' : 'MCG-01-01-012',
                   'PGC011767' : 'MCG-01-09-006',
                   'PGC013426' : 'MCG-01-10-015',
                   'PGC013535' : 'MCG-01-10-019',
                   'PGC065095' : 'MCG-01-52-012',
                   'PGC001841' : 'MCG-02-02-030',
                   'PGC064373' : 'MCG-02-51-004',
                   'UGC12463'  : 'NGC7559B',
                   'PGC072803' : 'NGC7783NED01',
                   'PGC002029' : 'UGC00335NED02',
                   'PGC066150' : 'UGC11680NED02',
                   'NGC7237'   : 'UGC11958',
                   'PGC070084' : 'VV488NED02'
                  }
    if list == 'dr3':
        renames = {
                   # LEDAname : CALIFAname
                   'PGC3089961': '2MASXJ09065870',
                   'NGC1143'   : 'ARP118',
                   'IC0225'    : 'BKD2008WR346',
                   'IC0701'    : 'IC701',
                   'UGC04299'  : 'IC2247',
                   'UGC09713'  : 'IC4534',
                   'PGC021757' : 'LSBCF560-04',
                   'NGC3406'   : 'NGC3406NED01',
                   'NGC4211A'  : 'NGC4211NED02',
                   'PGC000312' : 'IC1528',
                   'PGC049949' : 'NGC5421NED02',
                   'UGC08967'  : 'NGC5434B',
                   'NGC6027E'  : 'NGC6027',
                   'PGC058100' : 'NGC6150B',
                   'NGC6166'   : 'NGC6166NED01',
                   'PGC072803' : 'NGC7783NED01',
                   'UGC01370'  : 'SDSSJ015424',
                   'SDSSJ100142.11+371443.2' : 'SDSSJ100141.02+371447.4',
                   'NGC3655'   : 'SN2002ji',
                   'PGC001914' : 'UGC00312NED01',
                   'PGC002029' : 'UGC00335NED02',
                   'PGC029708' : 'UGC05498NED01',
                   'PGC066146' : 'UGC11680NED01',
                   'PGC066150' : 'UGC11680NED02',
                   'NGC7237'   : 'UGC11958',
                   'PGC071033' : 'UGC12494NOTES01',
                   'PGC059943' : 'VIIZw700',
                   'PGC070084' : 'VV488NED02'
                  }
    t.add_index('Name')
    for key in renames.keys():
        t.loc[key]['Name']=renames[key]

    # Calculated columns
    t['ledaD25'] = 0.1*10**t['logd25']
    t['ledaD25'].unit = 'arcmin'
    t['ledaD25'].description = 'Apparent B diameter from LEDA /logd25/ linearized'
    t['ledaD25'].info.format='.2f'

    t['ledaAxrat'] = 10**(-t['logr25'])
    t['ledaAxrat'].description = 'Minor to major axis ratio from LEDA /logr25/ linearized' 
    t['ledaAxrat'].info.format='.4f'

    # Override the bad axis ratio of logr25=0.46 for NGC 3687 in LEDA 
    # (Dmitry Makarov, priv comm, 29may2020).
    # Use values from SDSS-III DR12 (2015ApJS..219...12A)
    if list == 'edge' or list == 'dr3':
        t.loc['NGC3687']['ledaAxrat'] = 10**(0.815-0.916)
        t.loc['NGC3687']['ledaPA'] = 167.1

    # Missing redshift for PGC071033 in LEDA, use modz for UGC12494
    if list == 'dr3':
        t.loc['UGC12494NOTES01']['ledaModz'] = 33.95

    t['ledaAxIncl'] = phot_inc(t['ledaAxrat'],t['ledaType'])
    t['ledaAxIncl'].unit = 'deg'
    t['ledaAxIncl'].description = 'Inclination estimated from LEDA axratio using Bottinelli+83' 
    t['ledaAxIncl'].info.format='.1f'

    t['ledaDistMpc'] = 10**(0.2*(t['ledaModz']-25))
    t['ledaDistMpc'].unit = 'Mpc'
    t['ledaDistMpc'].description = 'Luminosity distance in Mpc corresponding to ledaModz' 
    t['ledaDistMpc'].info.format='.2f'

    t['ledaHIflux'] = 10**(-0.4*(t['m21']-17.4))
    t['ledaHIflux'].unit = 'Jy km / s'
    t['ledaHIflux'].description = '21cm flux from LEDA /m21/ linearized' 
    t['ledaHIflux'].info.format='.2f'

    # Descriptions
    t['Name'].description = 'Galaxy Name'

    t['ledaRA'].unit = 'hourangle'
    t['ledaRA'].convert_unit_to('deg')
    t['ledaRA'] = np.round(t['ledaRA'], 7) * u.deg
    t['ledaRA'].description = 'RA J2000 from LEDA /al2000/'

    t['ledaDE'].unit = 'deg'
    t['ledaDE'].description = 'DEC J2000 from LEDA /de2000/' 

    t['ledaA_Bgal'].unit = 'mag'
    t['ledaA_Bgal'].description = 'Galactic A_B from LEDA /ag/'

    t['ledaType'].description = 'Morphological type from LEDA /t/'

    t['ledaPA'].unit = 'deg'
    t['ledaPA'].description = 'PA from LEDA /pa/, N to E' 

    t['ledaIncl'].unit = 'deg'
    t['ledaIncl'].description = 'Morph inclination from LEDA /incl/' 

    t['ledaVrad'].unit = 'km / s'
    t['ledaVrad'].description = 'cz from mean data from LEDA /v/' 

    t['ledaVmaxg'].unit = 'km / s'
    t['ledaVmaxg'].description = 'HI max v_rot uncorrected for incl from LEDA /vmaxg/' 

    t['ledaVrot'].unit = 'km / s'
    t['ledaVrot'].description = 'HI max v_rot corrected for incl from LEDA /vrot/' 

    t['ledaMorph'].description = 'Hubble type from LEDA /type/' 
    t['ledaBar'].description = 'B = bar present from LEDA /bar/' 
    t['ledaRing'].description = 'R = ring present from LEDA /ring/' 
    t['ledaMultiple'].description = 'M = multiple system from LEDA /multiple/' 

    t['ledaBt'].unit = 'mag'
    t['ledaBt'].description = 'Apparent B total magnitude from LEDA /bt/' 

    t['ledaIt'].unit = 'mag'
    t['ledaIt'].description = 'Apparent I total magnitude from LEDA /it/' 

    t['ledaFIR'].unit = 'mag'
    t['ledaFIR'].description = 'FIR flux as magnitude from LEDA /mfir/' 

    t['ledaVvir'].unit = 'km / s'
    t['ledaVvir'].description = 'Virgo infall corrected cz from LEDA /vvir/' 

    t['ledaModz'].unit = 'mag'
    t['ledaModz'].description = 'Dist modulus from LEDA /modz/ based on /vvir/' 

    # Fill in missing PA values using NED
    badind = t['ledaPA'].mask.nonzero()[0]
    print('Indices with missing PA:',badind)
    for i, idx in enumerate(badind):
        newpa = getNedPA(t['Name'][idx])
        print('Setting PA for {} to {}'.format(t['Name'][idx],newpa))
        t['ledaPA'][idx] = newpa

    t.remove_columns(['logd25','logr25','e_logr25','m21'])
    t.pprint()
    tablist.append(t)

if len(tablist) > 0:
    alltab = vstack(tablist)
if len(idtab) > 0:
    newtab = join(alltab, idtab, keys='Name', join_type='left')
    newtab['ID'].description = 'CALIFA ID'
    alltab = newtab
alltab.add_index('Name')
alltab.sort('Name')
alltab.meta['date'] = datetime.today().strftime('%Y-%m-%d')
alltab.meta['comments'] = ('Galaxy Parameters from HyperLEDA')
print(alltab.meta)
alltab.write('edge_leda.csv', format='ascii.ecsv', delimiter=',', overwrite=True)


