from edge_pydb import conversion, fitsextract, xy2hist, util
from astropy.table import Table as _Table
from astropy.table import join as _join


class EdgeTable(_Table):
    def __init__(self, file='', path=''):
        super().__init__()
        if file and path:
            self.read(file, path)
        self.srcfile = file
        self.path = path
        self.joined = []
        
    def read(self, file, path):
        self.table = _Table.read(util.fetch(file), path=path)
        self.__dict__.update(self.table.__dict__)
        
    def join(self, file, join_type='inner', path=''):
        if 'csv' in file:
            target = _Table.read(util.fetch(file), format='ascii.ecsv')
        elif 'hdf5' in file:
            if not path:
                # using the same path as the source file
                path = self.path
            target = _Table.read(util.fetch(file), path)
        self.table = _join(self.table, target, join_type=join_type, keys='Name')
        # update the data
        self.__dict__.update(self.table.__dict__)
        self.joined.append((file, path, join_type))