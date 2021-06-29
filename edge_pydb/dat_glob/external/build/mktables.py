#!/usr/bin/env python

# Build external tables from Alberto's master table and NED output

from astropy.coordinates import SkyCoord
from astropy.table import Table, join
import numpy as np
from datetime import datetime
#import numpy.ma as ma

# Which tables to rewrite with this execution of the script
rewrite=['ned','wise','nsa']

# Alberto's master table, last updated April 14, 2017
t = Table.read('DETableFinal.csv', format='ascii.csv')

# Write the NED table:
if 'ned' in rewrite:
    nedt = Table.read('ned_diam.txt', format='ascii.csv', delimiter='\t')
    sc = SkyCoord(nedt['RA_J2000'],nedt['Dec_J2000']) #convert to degrees
    nedRA=[]
    nedDE=[]
    for obj in sc:
        nedRA.append(float(obj.to_string(precision=5).split(' ')[0]))
        nedDE.append(float(obj.to_string(precision=5).split(' ')[1]))
    newt=Table()
    newt['Name'] = nedt['InputName']
    newt['Name'].description = 'Galaxy Name'
    newt['nedRA'] = nedRA
    newt['nedRA'].unit = 'deg'
    newt['nedRA'].description = 'J2000 RA from NED'
    newt['nedDE'] = nedDE
    newt['nedDE'].unit = 'deg'
    newt['nedDE'].description = 'J2000 DEC from NED'
    newt['nedVopt']=nedt['cz']
    newt['nedVopt'].unit='km / s'
    newt['nedVopt'].description='cz from NED'
    newt['nedz']=nedt['Redshift']
    newt['nedz'].description='z from NED'
    newt['nedMajDia']=nedt['maj_diam_am']
    newt['nedMajDia'].unit='arcmin'
    newt['nedMajDia'].description='NED Basic Data major axis diameter'
    newt['nedMinDia']=nedt['min_diam_am']
    newt['nedMinDia'].unit='arcmin'
    newt['nedMinDia'].description='NED Basic Data minor axis diameter'
    newt.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    newt.meta['comments'] = ('Galaxy coordinates and diameters from NED')
    print(newt.meta)
    newt.write('edge_ned.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# Write the WISE table:
t['W1'].unit = 'mag'
t['W1'].description = '3.4 um Vega magnitude from image photometry from table by Bitsakis. NGC0598 and NGC4676A are from w1mpro in allWISE catalog'
t['eW1'].unit = 'mag'
t['eW1'].description = 'Error in W1 in Vega magnitudes from table by Bitsakis. NGC0598 and NGC4676A are from w1mpro+2.5*log(w1snr) in allWISE catalog'
t['W2'].unit = 'mag'
t['W2'].description = '4.6 um Vega magnitude from Bitsakis'
t['eW2'].unit = 'mag'
t['eW2'].description = 'Error in W2 in Vega magnitudes from Bitsakis'
t['W3'].unit = 'mag'
t['W3'].description = '12 um Vega magnitude from Bitsakis'
t['eW3'].unit = 'mag'
t['eW3'].description = 'Error in W3 in Vega magnitudes from Bitsakis'
t['W4'].unit = 'mag'
t['W4'].description = '22 um Vega magnitude from Bitsakis'
t['eW4'].unit = 'mag'
t['eW4'].description = 'Error in W4 in Vega magnitudes from Bitsakis'
t['W1lum'].unit = 'dex(erg / s)'
t['W1lum'].description = 'Luminosity in W1 from caDistMpc and W1 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W2lum'].unit = 'dex(erg / s)'
t['W2lum'].description = 'Luminosity in W2 from caDistMpc and W2 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W3lum'].unit = 'dex(erg / s)'
t['W3lum'].description = 'Luminosity in W3 from caDistMpc and W3 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W4lum'].unit = 'dex(erg / s)'
t['W4lum'].description = 'Luminosity in W4 from caDistMpc and W4 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W4SFR'].unit = 'solMass / yr'
t['W4SFR'].description = 'SFR from W4lum using Catalan-Torrecilla+14 calibration'
if 'wise' in rewrite:
    outcols=['Name']
    clist=['W1', 'eW1', 'W2', 'eW2', 'W3', 'eW3', 'W4', 'eW4', 'W1lum', 'W2lum', 'W3lum', 'W4lum', 'W4SFR']
    for cname in clist:
        outcols.append(cname)
    for cname in ['W1lum', 'W2lum', 'W3lum', 'W4lum']:
        logged = np.around(np.log10(t[cname]),decimals=4)
        t[cname][:] = logged
    newt = t[outcols]
    newt.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    newt.meta['comments'] = ('WISE photometry from Bitsakis')
    print(newt.meta)
    newt.write('edge_wise.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# Write the NSA table:
t['nsaZdist'].description = 'Redshift from NASA-Sloan Atlas using pecular velocity model of Willick+ 97'
t['nsaAu'].unit = 'mag'
t['nsaAu'].description = 'Galactic extinction in u from Schlegel+ 97'
t['nsaAg'].unit = 'mag'
t['nsaAg'].description = 'Galactic extinction in g from Schlegel+ 97'
t['nsaAr'].unit = 'mag'
t['nsaAr'].description = 'Galactic extinction in r from Schlegel+ 97'
t['nsaAi'].unit = 'mag'
t['nsaAi'].description = 'Galactic extinction in i from Schlegel+ 97'
t['nsaAz'].unit = 'mag'
t['nsaAz'].description = 'Galactic extinction in z from Schlegel+ 97'
if 'nsa' in rewrite:
    outcols=['Name']
    for cname in t.colnames:
        if cname.startswith("nsa"):
            outcols.append(cname)
    newt = t[outcols]
    newt.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    newt.meta['comments'] = ('Galaxy Parameters from NASA-Sloan Atlas')
    print(newt.meta)
    newt.write('edge_nsa.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

