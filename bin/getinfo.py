import os, glob
import h5py
import astropy
from astropy.table import Table

def getPath(file):
    h5f = h5py.File(file, 'r')
    return [key for key in h5f.keys() if "__table_column_meta__" not in key]

h5files = glob.glob('*.hdf5')

for h5f in h5files:
    output = h5f+'.info.txt'
    h5info = open(output, 'w')
    for path in getPath(h5f):
        h5info.write("filename: {}\npath: {}\n".format(h5f, path))
        tab = Table.read(h5f, path=path)
        astropy.table.info.table_info(tab, out=h5info)
        h5info.write("\n")
    h5info.close()

