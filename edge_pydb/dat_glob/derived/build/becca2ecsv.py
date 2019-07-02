#!/usr/bin/env python

from astropy.table import Table
import numpy as np

infile = "EDGE_COparameters_20180604.csv"

# Get header names from top of file
headnames = {}
with open(infile) as f:
    for line in f:
        if ":" in line:
            str=line.split(':',1)
            headnames[str[0]] = str[1].strip()

for key in headnames:
    print(key)

# Get data, add 'rf' string before column names
t = Table.read(infile, format='ascii.csv', header_start=18, data_start=19)
t.pprint()
for cname in t.colnames:
    t[cname].description = headnames[cname]
    #t[cname].name = t[cname].name

t.write('edge_rfpars.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

