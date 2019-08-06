from edge_pydb import conversion, fitsextract, xy2hist, util
from astropy.table import Table as _Table
from astropy.table import join as _join


class EdgeTable(_Table):
    def __init__(self, file='', path='', cols=None):
        super().__init__()
        if file:
            if file.endswith('csv'):
                self.read(file)
            elif path:
                self.read(file, path)
        else:
            print("Choose from the following files to read:")
            util.listfiles(printing=True)
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
            self.table = _Table.read(util.fetch(file), format='ascii.ecsv')
        elif path:
            self.table = _Table.read(util.fetch(file), path=path)
        self.__dict__.update(self.table.__dict__)
        
    def join(self, file='', table=None, join_type='inner', path=''):
        if table:
            self.table = _join(self.table, table.table)
            # update the data
            self.__dict__.update(self.table.__dict__)
            return
        if 'csv' in file:
            target = _Table.read(util.fetch(file), format='ascii.ecsv')
        elif 'hdf5' in file:
            if not path:
                # using the same path as the source file
                path = self.path
            if file.endswith('csv'):
                target = _Table.read(util.fetch(file), path=path)
            elif path:
                target = _Table.read(util.fetch(file), path=path)
        else:
            # raise the error
            target = None
        self.table = _join(self.table, target, join_type=join_type, keys='Name')
        # update the data
        self.__dict__.update(self.table.__dict__)
        if file.endswith('csv'):
            self.joined.append((file, join_type))
        else:
            self.joined.append((file, path, join_type))