#!/usr/bin/env python

from astropy.table import Table
import numpy as np
from datetime import datetime

infile = "jam.csv"

# Get data
t = Table.read(infile, format='ascii.csv')
t.pprint()

t['Name'].description = 'Galaxy Name'
t['jaType'].description = 'Morphological Type'
t['jaDistMpc'].description = 'Distance in Mpc'
t['jaDistMpc'].unit = 'Mpc'
t['jaIncl'].description = 'Inclination from outer isophotes'
t['jaIncl'].unit = 'deg'
t['jaMstar'].description = 'Total Stellar mass'
t['jaMstar'].unit = 'dex(solMass)'
t['jaMorphPA'].description = 'PA from outer isophotes'
t['jaMorphPA'].unit = 'deg'
t['jaKinPA'].description = 'PA from fitting CO kinematics'
t['jaKinPA'].unit = 'deg'
t['jaReff'].description = 'Effective radius'
t['jaReff'].unit = 'arcsec'
t['jaVstar'].description = 'Stellar velocity dispersion at Reff'
t['jaVstar'].unit = 'km / s'
t['jaVsys'].description = 'Systemic velocity referenced to optical convention cz'
t['jaVsys'].unit = 'km / s'

t.meta['date'] = datetime.today().strftime('%Y-%m-%d')
t.meta['comments'] = ('Galaxy Parameters from Table 1 of 2018MNRAS.477..254L')
print(t.meta)
t.write('edge_jampars.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

