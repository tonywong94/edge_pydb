from edge_pydb import conversion, fitsextract, xy2hist, util
from astropy.table import Table as _Table
from astropy.table import join as _join


class EdgeTable:
    def __init__(self, file, path):
        self.table = _Table.read(util.fetch(file), path)
        self.path = path
        self.srcfile = file

    def join(self, file, join_type='inner', path=''):
        if 'csv' in file:
            target = _Table.read(util.fetch(file), format='ascii.ecsv')
        elif 'hdf5' in file:
            if not path:
                # using the same path as the source file
                path = self.path
            target = _Table.read(util.fetch(file), path)
        self.table = _join(self.table, target, join_type=join_type, keys='Name')
        
