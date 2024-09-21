from ._astropy_init import __version__, test
from . import beam_sample, conversion, fitsextract, plotting, util
from astropy.table import Table as _Table
from astropy.table import join as _join


'''
Definition of the EdgeTable class.  An EdgeTable is a regular AstroPy table in 
either ECSV or HDF5 format.  For an HDF5 table, the Path must be specified, otherwise 
a list of available Paths is provided.
Keyword 'cols' allows a subset of available columns to be read in.
'''


class EdgeTable(_Table):
    """
    Definition of the EdgeTable class.  An EdgeTable is a regular Astropy table in 
    either ECSV or HDF5 format.

    Parameters
    ----------
    file : str
        Name of the file.  This is generally NOT a full path, as the package
        maintains a list of available tables.  The name should end with .csv
        or .hdf5.  Use 'list' to provide a listing of available tables.
    path : str
        Name of the path within the file, for HDF5 files only.  If not given,
        a list of available paths is output.
    cols : list of str
        A list of columns to read in.  Other columns are discarded.
    """

    def __init__(self, file='', path='', cols=None, data=None, masked=None, names=None, 
                 dtype=None, meta=None, copy=True, rows=None, copy_indices=True,
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
                print('Paths in',file,':\n', util.getPath(file))
        if cols:
            self.table.keep_columns(cols)
#             data = []
#             for i in cols:
#                 data.append(self.table[i])
#             self.table = _Table(data=data)
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
        
    def join(self, table, join_type='inner', keys=None):
        # check for ix or iy in both tables 
        if keys is None:
            keys = ['Name']
        join_keys = [i for i in keys]
        if 'ix' in self.colnames and 'ix' in table.colnames:
            join_keys.append('ix')
            if 'iy' in self.colnames and 'iy' in table.colnames:
                join_keys.append('iy')
        if isinstance(table, self.__class__):
            self.table = _join(self.table, table.table, join_type=join_type, keys=join_keys)
            self.joined.append((table.srcfile, join_type))
        elif isinstance(table, _Table):
            self.table = _join(self.table, table, join_type=join_type, keys=join_keys)
            self.joined.append((table.srcfile, join_type))
        else:
            raise Exception('cannot merge the two tables, \
                the second table data type is neither EdgeTable nor astropy table.')
        # update the data
        self.__dict__.update(self.table.__dict__)


