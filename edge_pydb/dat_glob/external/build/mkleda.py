#!/usr/bin/env python

from astropy.table import Table, join
import numpy as np
from datetime import datetime
#import numpy.ma as ma

#list = 'alma'
list = 'edge'

def phot_inc(axratio, ttype):
    q = axratio  # minor/major axis ratio
    q0 = 10**(-0.60 - 0.045*ttype)  # Bottinelli+ 83 (1983A&A...118....4B)
    # For 10 galaxies where the observed q < q_0, we replace q_0 by the smaller 
    #   of 0.135 and q to avoid an undefined value for cos i. 
    #   q_0=0.135 corresponds to ttype=6, close to the latest types in EDGE.
    q0[q<q0] = np.minimum(0.135,q[q<q0])
    inc = np.degrees(np.arccos( np.sqrt((q**2-q0**2)/(1-q0**2)) ))  # Hubble 1926ApJ....64..321H
    return inc
    
# Process the output from HyperLEDA
t = Table.read(list+'_ledaout.txt',format='ascii.fixed_width',header_start=2)
t['name'] = t['name'].astype('|S15')

# Rename the LEDA columns
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
for key in colnames.keys():
    t.rename_column(key, colnames[key])

# Some CALIFA names were unrecognized by LEDA and had to be manually substituted
if list == 'edge':
    renames = {
               'NGC4211A' : 'NGC4211NED02',
               'NGC6027E' : 'NGC6027',
               'PGC029708' : 'UGC05498NED01',
               'UGC04299' : 'IC2247'
              }
if list == 'alma':
    renames = {
               'PGC000312' : 'IC1528',
               'PGC072803' : 'NGC7783NED01',
               'PGC002029' : 'UGC00335NED02',
               'PGC066150' : 'UGC11680NED02',
               'NGC7237' : 'UGC11958',
               'PGC070084' : 'VV488NED02'
              }
t.add_index('Name')
for key in renames.keys():
    t.loc[key]['Name']=renames[key]


# Calculated columns
t['ledaD25'] = 0.1*10**t['logd25']
t['ledaD25'].unit = 'arcmin'
t['ledaD25'].description = 'Apparent B diameter from LEDA /logd25/ linearized'
t['ledaD25'].format='.2f'

t['ledaAxrat'] = 10**(-t['logr25'])
t['ledaAxrat'].description = 'Minor to major axis ratio from LEDA /logr25/ linearized' 
t['ledaAxrat'].format='.4f'

# Override the bad axis ratio for NGC 3687 in LEDA (Dmitry Makarov, priv comm, 29may2020).
# Use values from SDSS-III DR12 (2015ApJS..219...12A)
t.loc['NGC3687']['ledaAxrat'] = 10**(0.815-0.916)
t.loc['NGC3687']['ledaPA'] = 167.1

t['ledaAxIncl'] = phot_inc(t['ledaAxrat'],t['ledaType'])
t['ledaAxIncl'].unit = 'deg'
t['ledaAxIncl'].description = 'Inclination estimated from LEDA axratio using Bottinelli+83' 
t['ledaAxIncl'].format='.1f'

t['ledaDistMpc'] = 10**(0.2*(t['ledaModz']-25))
t['ledaDistMpc'].unit = 'Mpc'
t['ledaDistMpc'].description = 'Luminosity distance in Mpc corresponding to ledaModz' 
t['ledaDistMpc'].format='.2f'

t['ledaHIflux'] = 10**(-0.4*(t['m21']-17.4))
t['ledaHIflux'].unit = 'Jy km / s'
t['ledaHIflux'].description = '21cm flux from LEDA /m21/ linearized' 
t['ledaHIflux'].format='.2f'

# Descriptions
t['Name'].description = 'Galaxy Name'

t['ledaRA'].unit = 'hourangle'
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

# Fill in missing PA values using 2MASS XSC
badind = t['ledaPA'].mask.nonzero()[0]
print(badind)
# Should be NGC447, NGC2487, NGC5394, NGC 5784, NGC5953, NGC6125
newpa = [15, 50, 30, 40, 65, 50]
for i, idx in enumerate(badind):
    print('Setting PA for {} to {}'.format(t['Name'][idx],newpa[i]))
    t['ledaPA'][idx] = newpa[i]

t.remove_columns(['objname', 'logd25', 'logr25', 'e_logr25', 'm21'])
t.meta['date'] = datetime.today().strftime('%Y-%m-%d')
print(t.meta)
t.write(list+'_leda.csv', format='ascii.ecsv', delimiter=',', overwrite=True)


