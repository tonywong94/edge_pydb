#!/usr/bin/env python

# Radial profile parameters from Bolatto et al. (2017)

from astropy.table import Table
from datetime import datetime

# Alberto's master table, last updated April 14, 2017
t = Table.read('../../external/build/DETableFinal.csv', format='ascii.csv')

# Write the radial distributions table (edge_rdist.csv):
t['coScaleMol'].unit = 'kpc'
t['coScaleMol'].description = 'Exponential scale length for CO disk derived by filling in undetected values in annulii with 1-sigma'
t['coeScaleMol'].unit = 'kpc'
t['coeScaleMol'].description = 'Statistical error from fit to exponential scale length'
t['coScaleMolHi'].unit = 'kpc'
t['coScaleMolHi'].description = 'Upper limit to exponential scale length by filling in annuli with 2-sigma values'
t['coScaleMolLo'].unit = 'kpc'
t['coScaleMolLo'].description = 'Lower limit to exponential scale length by filling in annulii with zeros'
t['coNormMol'].unit = 'solMass / pc2'
t['coNormMol'].description = 'Normalization of CO exponential disk profile, i.e. density at R=0'
t['coeNormMol'].unit = 'solMass / pc2'
t['coeNormMol'].description = 'Error in normalization of CO exponential disk profile'
t['coScaleSt'].unit = 'kpc'
t['coScaleSt'].description = 'Exponential scale length for the mass of the stellar disk'
t['coeScaleSt'].unit = 'kpc'
t['coeScaleSt'].description = 'Formal error in stellar scale length fit' 
t['coNormSt'].unit = 'solMass / pc2'
t['coNormSt'].description = 'Normalization of stellar exponential disk profile, i.e. density at R=0'
t['coeNormSt'].unit = 'solMass / pc2'
t['coeNormSt'].description = 'Formal error in stellar disk normalization'
t['coR50Mol'].unit = 'kpc'
t['coR50Mol'].description = 'Radius enclosing 50% of the molecular mass'
t['coeR50Mol'].unit = 'kpc'
t['coeR50Mol'].description = 'Error in radius enclosing 50% of the molecular mass, including beam size'
t['coR50St'].unit = 'kpc'
t['coR50St'].description = 'Radius enclosing 50% of the stellar mass'
t['coeR50St'].unit = 'kpc'
t['coeR50St'].description = 'Error in radius enclosing 50% of the stellar mass, including beam size'
t['coScaleSFR'].unit = 'kpc'
t['coScaleSFR'].description = 'Exponential scale length for SFR from extinction corrected Ha'
t['coeScaleSFR'].unit = 'kpc'
t['coeScaleSFR'].description = 'Error in exponential scale length for SFR from extinction corrected Ha'

outcols=['Name']
clist = ['coScaleMol', 'coeScaleMol', 'coScaleMolHi', 'coScaleMolLo', 'coNormMol', 
         'coeNormMol', 'coScaleSt', 'coeScaleSt', 'coNormSt', 'coeNormSt', 'coR50Mol', 
         'coeR50Mol', 'coR50St', 'coeR50St', 'coScaleSFR', 'coeScaleSFR']
newt = t[outcols+clist]
for cname in clist:
    newt[cname].name = t[cname].name.replace('co','rd')
newt.meta['date'] = datetime.today().strftime('%Y-%m-%d')
newt.meta['comments'] = ('Radial distribution parameters from Bolatto+ 2017ApJ...846..159B')
print(newt.meta)
newt.write('edge_rdist17.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

