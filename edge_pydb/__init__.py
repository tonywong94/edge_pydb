from edge_pydb import conversion, fitsextract, xy2hist, util
from astropy.table import Table as _Table
from astropy.table import join as _join
import h5py as _h5py

class EdgeTable(_Table):
    def __init__(self, file='', path='', cols=None, data=None, masked=None, names=None, dtype=None,
                 meta=None, copy=True, rows=None, copy_indices=True,
                 **kwargs):
        super().__init__(masked=masked)
        if file:
            if file == 'list':
                print("Choose from the following files to read:")
                util.listfiles(printing=True)
            elif file.endswith('csv'):
                self.read(file)
            elif path:
                self.read(file, path)
            else:
                # no path specified with hdf5 file
                # f = _h5py.File(util.fetch(file), 'r')
                print('Paths in',file,':\n', util.getPath(file))
#         else:
#             print("Choose from the following files to read:")
#             util.listfiles(printing=True)
        if cols:
            data = []
            for i in cols:
                data.append(self.table[i])
            self.table = _Table(data=data)
            self.__dict__.update(self.table.__dict__)
        self.srcfile = file
        self.path = path
        self.joined = []
        
    def read(self, file, path=''):
        if 'csv' in file:
            try:
                self.table = _Table.read(util.fetch(file), format='ascii.ecsv')
            except ValueError:
                self.table = _Table.read(util.fetch(file), format='ascii.csv')
        elif path:
            self.table = _Table.read(util.fetch(file), path=path)
        self.__dict__.update(self.table.__dict__)
        
    def join(self, table, join_type='inner', keys='Name'):
        # if table:
        #     self.table = _join(self.table, table.table)
        #     # update the data
        #     self.__dict__.update(self.table.__dict__)
        #     return
        # if 'csv' in file:
        #     target = _Table.read(util.fetch(file), format='ascii.ecsv')
        # elif 'hdf5' in file:
        #     if not path:
        #         # using the same path as the source file
        #         path = self.path
        #     if file.endswith('csv'):
        #         target = _Table.read(util.fetch(file), path=path)
        #     elif path:
        #         target = _Table.read(util.fetch(file), path=path)
        # else:
        #     # raise the error
        #     target = None
        if isinstance(table, _Table):
            self.table = _join(self.table, table, join_type=join_type, keys=keys)
        elif isinstance(table, self.__class__):
            self.table = _join(self.table, table.table, join_type=join_type, keys=keys)
            self.joined.append((table.srcfile, join_type))
        # update the data
        self.__dict__.update(self.table.__dict__)




