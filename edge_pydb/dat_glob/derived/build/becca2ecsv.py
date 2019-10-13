#!/usr/bin/env python

from astropy.table import Table
import numpy as np
from datetime import datetime

infile = "EDGE_COparameters_20191004.csv"
head_start = 21  # previously 18
dat_start = 22   # previously 19

# Get header names from top of file
headnames = {}
with open(infile) as f:
    for line in f:
        if ":" in line:
            str=line.split(':',1)
            headnames[str[0]] = str[1].strip().strip(",")

print('Information found for {} columns'.format(len(headnames)))
for key in headnames:
    print(key)

# Get data, add 'rf' string before column names
t = Table.read(infile, format='ascii.csv', header_start=head_start, data_start=dat_start)
t.pprint()
for cname in t.colnames:
    t[cname].description = headnames[cname]

t['rfVsys'].unit = 'km / s'
t['rfVsys'] = t['rfVsys'].astype(float)
t['rfLSRK2helio'].unit = 'km / s'
t['rfRflat'].unit = 'arcsec'
t['rfVrotMaxL'].unit = 'km / s'
t['rfeVrotMaxL'].unit = 'km / s'
t['rfVrotMax'].unit = 'km / s'
t['rfeVrotMax'].unit = 'km / s'
t['rfPA'].unit = 'deg'
t['rfInc'].unit = 'deg'
t['rfKinXoff'].unit = 'arcsec'
t['rfKinYoff'].unit = 'arcsec'
t['rfKinRA'].unit = 'hourangle'
t['rfKinDecl'].unit = 'deg'

for extracol in ['ledaRA','ledaDE']:
    if extracol in t.colnames:
        t.remove_column(extracol)

t.meta['date'] = datetime.today().strftime('%Y-%m-%d')
print(t.meta)
t.write('edge_rfpars.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

