from edge_pydb.plotting import gridplot as _gridplot
from matplotlib.colors import Normalize as _Normalize
import numpy as _np
from astropy.table import Table as _Table
import astropy as _astropy
import sys
import os as _os
import json as _json
import shutil as _shutil
import requests as _requests
import h5py as _h5py


# Initial setup script to read the file location from config file
_ROOT = _os.path.abspath(_os.path.dirname(__file__))

_filepath = _os.path.join(_ROOT, '_config.json')
_runtime = False
_config = {}

try:
    _fp = open(_filepath, 'r')
    if _os.stat(_filepath).st_size > 0:
        _config = _json.load(_fp)

except FileNotFoundError:
    _fp = open(_filepath, 'w')
    # _config = {}

except OSError as _err:
    if _err.errno == 30:
        print("WARNING! Read-only file system, cannot save the package data file location.\n" +
              "For better longterm performance, consider providing a config file by using the extConfig() function.")
        print("If you need to change the files in the package data, please consider running as root, \
        the manipulation of files requires the sudo priority.")
    _runtime = True
    # _config = {}


def _walkthrough(dir=_ROOT, max_depth=2):
    retval = {}
    dir = _os.path.abspath(dir)
    base_depth = dir.rstrip(_os.path.sep).count(_os.path.sep)
    for _root, _dirs, _files in _os.walk(dir):
        cur_depth = _root.count(_os.path.sep)
        if max_depth >= 0 and cur_depth > max_depth + base_depth:
            continue
        for _file in _files:
            if _file.endswith('.csv') or _file.endswith('.hdf5'):
                if _file in retval:
                    print("{} redundant file detected\n--Current location: {}\n++New location: {}\n".format(
                        _file, retval[_file], _os.path.join(_root, _file)))
                retval[_file] = _os.path.join(_root, _file)
    return retval


if not _config:
    # print(os.listdir(_ROOT))
    _config = _walkthrough()
    if not _runtime:
        _json.dump(_config, _fp)

if not _runtime:
    _fp.close()


# update the files
def updatefiles(dir=_ROOT, max_depth=-1):
    tmp = _walkthrough(dir, max_depth)
    for k, v in tmp.items():
        if k not in _config.keys():
            print("Add file %s" % k)
        elif v != _config[k]:
            print("Update file %s" % k)
    _config.update(tmp)
    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)


def download(file, url, loc='', user='', password=''):
    if not loc:
        if not _os.path.exists(_ROOT + '/data'):
            _os.mkdir(_ROOT + '/data')
        loc = _ROOT + '/data/'
    data = _requests.get(url + file, verify=False, auth=(user, password))
    if data.status_code != 200:
        data.raise_for_status()
    if not loc.endswith('/'):
        loc += '/'
    with open(loc + file, "wb") as _fp:
        _fp.write(data.content)


def save_config(src):
    '''
    This function will write _config back to a file.
    '''
    global _config
    with open(src, 'w') as _fp:
        _json.dump(_config, _fp)


def load_config(src, readonly=False):
    '''
    read the config from a file
    '''
    _fp = open(src, 'r')
    global _config
    if readonly:
        _config = _json.load(_fp)
    else:
        _config.update(_json.load(_fp))
        _filepath = src
        _runtime = False
    _fp.close()


def listfiles(contain='', values=False, printing=False):
    '''
    List the current available files in the package data directory
    If values=True, give full directory paths
    If contain='hdf', list only file names with the 'hdf' substring

    Parameters:
        contain: the target substring to find in the file name to list, if not provided, then will print all files
    '''
    files = []
    if values:
        for val in _config.values():
            if contain in val:
                if printing:
                    print(val)
                files.append(val)
    else:
        for key in _config.keys():
            if contain in key:
                if printing:
                    print(key)
                files.append(key)
    return files


def fetch(names):
    '''
    Get all the files by its file name, can either be a single file or a list of files

    Parameters:
        names: a list of files or a single file name in package data
    '''
    if isinstance(names, list):
        retval = []
        for name in names:
            if name not in _config.keys():
                raise FileNotFoundError(
                    "Cannot find the specified file: %s" % name)
            # if dir:
            #     dirpath = _os.path.abspath(_os.path.dirname(name))
            #     if dirpath not in retval:
            #         retval.append(path)
            else:
                retval.append(_config[name])
        return retval
    else:
        if names not in _config.keys():
            raise FileNotFoundError(
                "Cannot find the specified file: %s" % names)
        # if dir:
        #     return _os.path.abspath(_os.path.dirname(names))
        else:
            return _config[names]


def addfile(src, dest='', copy=True, overwrite=False):
    if _runtime:
        print("WARNING! No sudo permission, take care, will break")
    name = _os.path.basename(src)
    src = _os.path.abspath(src)
    if copy:
        if dest:
            if name in _config.keys():
                if overwrite:
                    _os.remove(_config[name])
                else:
                    raise FileExistsError(_config[name])
            _shutil.copyfile(src, dest)
            _config[name] = dest
        else:
            if name in _config.keys():
                if overwrite:
                    _shutil.copyfile(src, _config[name])
                else:
                    raise FileExistsError('%s exists: ' % name + _config[name])

            if not _os.path.exists(_ROOT + '/data'):
                _os.mkdir(_ROOT + '/data')
            _shutil.copyfile(src, _ROOT + '/data/' + name)
            _config[name] = _ROOT + '/data/' + name
    else:
        if name in _config.keys():
            if overwrite:
                _os.remove(_config[name])
            else:
                raise FileExistsError('%s exists: ' % name + _config[name])
        _config[name] = src

    print("Update file %s" % name)
    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)
    else:
        print("WARNING! The location of this file will be saved runtime only")


def add_from_dir(src, dest='', copy=True, overwrite=False, max_depth=-1):
    '''
    This function copy files from a src directory to dest directory,
    if dest is empty, then it creates a directory just under the edge_pydb
    package directory. The copy always assumes a topdown copy (i.e. from parent
    directory to child directories)

    Parameters
    ----------
    src: source directory to copy from
    dest: destination directory to copy to
    copy: if copy is false, files will not be copied, and instead the path of
            these files will be recorded in the _config.json
    overwrite: if file is at the destination, overwrite the file if true, else create
            and copy into a subdirectory data/ under the edge_pydb package directory
    max_depth: specify the depth the copy should perform.
            -1 means copy all directories
            0 means just under the root directory and do not go into subdirectories
    '''
    # this function assume a topdown copy
    if _runtime:
        print("WARNING! No sudo permission, take care, will break")
    dirname = _os.path.basename(src)
    # dirname = _os.path.abspath(_os.path.dirname(src))
    if not dest:
        dest = _ROOT + '/' + dirname

    if copy:
        if max_depth >= 0:
            src = _os.path.normpath(src)
            tmp = _walkthrough(src, max_depth)
            for k, v in tmp.items():
                # copy2 will not raise FileExist, but overwrite directly
                _shutil.copy2(v, dest)
                # print("here ", k, v)
        else:
            try:
                _shutil.copytree(src, dest)
            except FileExistsError:
                if overwrite:
                    _shutil.rmtree(dest)
                    _shutil.copytree(src, dest)
                else:
                    dest = _ROOT + '/data/' + dirname
                    _shutil.copytree(src, dest)
        updatefiles(dest)
    else:
        updatefiles(src, max_depth)

    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)
    else:
        print("WARNING! The location of this file will be saved runtime only")


def getPath(file):
    '''get the hdf5 path'''
    h5f = _h5py.File(fetch(file), 'r')
    return [key for key in h5f.keys() if "__table_column_meta__" not in key]


def md_generate(csv_output, h5_output):
    """Generate markdown file for csv and txt for hdf5"""
    csvfiles = open(csv_output, 'w')
    h5files = open(h5_output, 'w')
    files = listfiles()
    title = ""
    for file in files:
        if file.endswith(".csv"):
            title += "- [" + file + "]" + "(#" + file.replace('.', '') + ")\n"
    csvfiles.write(title + "\n\n")
    for file in files:
        other_info = ""
        if file.endswith(".csv"):
            # print(file)
            # check the ecsv valid
            with open(fetch(file), 'r') as fp:
                lines = fp.readlines()
                if "ECSV" in lines[0]:
                    print("Working on {}".format(file))
                    name = "## {}\n\n".format(file)
                    comment = ""
                    header = "| name | unit | datatype | format | description |\n|---|---|---|---|---|\n"
                    to_print = header
                    aux = []
                    for i in range(len(lines)):
                        line = lines[i]
                        second_line = ""
                        if "# - {" not in line:
                            continue
                        if line[-2] != "}" and line[-1] == "\n":
                            line = line + " "
                            second_line = lines[i+1][5:-2]
                            i += 1
                        params = line_proc(line, aux, 4)
                        output = [" " for j in range(5)]
                        for param in params:
                            if param[0] == "name":
                                output[0] = param[1]
                            elif param[0] == "unit":
                                output[1] = param[1]
                            elif param[0] == "datatype":
                                output[2] = param[1]
                            elif param[0] == "format":
                                output[3] = param[1]
                            elif param[0] == "description":
                                output[4] = param[1] + second_line
                            elif param[0] == "comments":
                                comment = param[1]
                            else:
                                other_info = param[0] + ": " + param[1]
                        if output[0] != " ":
                            to_print += "| {} | {} | {} | {} | {} |\n"\
                                .format(output[0], output[1], output[2], output[3], output[4])
                    if other_info:
                        other_info += "\n\n"
                    csvfiles.write(name + comment + "\n\n" +
                                   other_info + to_print + "\n")
        elif file.endswith(".hdf5"):
            # h5files.write("{}\n".format(file))
            for path in getPath(file):
                h5files.write("filename: {}\npath: {}\n".format(file, path))
                tab = _Table.read(fetch(file), path=path)
                _astropy.table.info.table_info(tab, out=h5files)
                h5files.write("\n")
            h5files.write("\n")
    csvfiles.close()
    h5files.close()


def line_proc(line, other, val):
    cp = line[5:-2]
    newline = cp.split(", ", val)
    params = []
    for i in range(len(newline)):
        seg = newline[i]
        if "description: " in seg:
            i += 1
            while i < len(newline):
                seg += ", " + newline[i]
                i += 1

        tmp = seg.split(": ", 2)
        if len(tmp) == 2:
            params.append((tmp[0], tmp[1]))
        else:
            other += tmp
    return params


def add_url(file, root_url="https://github.com/tonywong94/edge_pydb/blob/master"):
    with open(file, 'r') as fp:
        lines = fp.readlines()
    for i in range(len(lines)):
        if "##" in lines[i]:
            fname = lines[i].rstrip().split(' ')[1]
            url = root_url
            substr = fetch(fname).split('/')
            flag = False
            for j in range(len(substr)):
                if substr[j] == 'edge_pydb' and ('edge_pydb' not in substr[j+1:-1]):
                    flag = True
                if flag:
                    url += '/' + substr[j]
            lines[i] = "## [{}]({})\n".format(fname, url.rstrip())
    with open(file, 'w') as fp:
        fp.writelines(lines)


def to_markdown(csv_out='index_csv.md', h5_out='index_hdf.txt', add_url=True):
    md_generate(csv_out, h5_out)
    if add_url:
        add_url(csv_out)
    return


def plotgallery(hdf_files=None, scale='auto', nx=7, ny=6, pad=8, 
                minperc=1, maxperc=99, paths=None, basedir='.'):
    '''
    Make multi-page gridplots for all galaxies in all available HDF5 files.

    === Parameters ===
    hdf_files : list of str
        Names of HDF5 files, should be available via EdgeTable.  Default is
        to process all available HDF5 files except with 'cocube' in the filename.
    scale : str
        'auto' (default) uses a jet (rainbow) colormap normalized to each galaxy
            individually.
        'perc' uses a gist_ncar_r colormap ranging from 1st to 99th percentile
            over all galaxies in the plotted column.
    nx : int
        number of subplots in horizontal direction
    ny : int
        number of subplots in vertical direction
    pad : int
        Padding in pixels around edges of bounding box
    minperc : float
        Minimum percentile for scale='perc'.  Default is 1%.
    maxperc : float
        Maximum percentile for scale='perc'.  Default is 99%.
    paths: list of str
        Names of paths (subtables) to plot.  Default is to plot all.
    basedir : str
        The directory into which to write the files.
    '''

    # Get the list of available files
    if hdf_files is None:
        hdf_files = listfiles(contain='hdf')
        print('Files to be processed:\n{}'.format(hdf_files))
    elif isinstance(hdf_files, str):
        hdf_files = [hdf_files]

    # Loop over files
    for dofile in hdf_files:
        if 'cocube' in dofile:
            continue
        # Loop over paths within each file
        if paths is None:
            paths = getPath(dofile)
            print('\nPaths in {}:\n{}'.format(dofile, paths))
        for dopath in paths:
            tab = _Table.read(fetch(dofile), path=dopath)
            print('\nWorking on {}'.format(dopath))
            for j in range(9, len(tab.colnames)):
                if tab.colnames[j] == 'cosi':
                    continue
                if scale == 'perc':
                    vmin = _np.nanpercentile(
                        tab[tab.colnames[j]], minperc, interpolation='nearest')
                    vmax = _np.nanpercentile(
                        tab[tab.colnames[j]], maxperc, interpolation='nearest')
                    if vmax == vmin:
                        vmax = vmin + 1
                    print('\n{} has vmin={} and vmax={}'.format(
                        tab.colnames[j], vmin, vmax))
                    norm = _Normalize(vmin=vmin, vmax=vmax)
                    cm = 'nipy_spectral'
                    outfile = _os.path.join(
                        basedir, dofile, dopath, tab.colnames[j]+'_perc.pdf')
                else:
                    print('')
                    norm = None
                    cm = 'jet'
                    outfile = _os.path.join(
                        basedir, dofile, dopath, tab.colnames[j]+'_auto.pdf')
                if not _os.path.isdir(_os.path.join(basedir, dofile, dopath)):
                    _os.makedirs(_os.path.join(basedir, dofile, dopath))
                if 'hex' in dofile:
                    _gridplot(edgetab=tab, columnlist=tab.colnames[j], vshow=True,
                              plotstyle='dot', clipedge=True, pad=pad, nx=nx, ny=ny,
                              cmap=cm, norm=norm, pdfname=outfile)
                else:
                    _gridplot(edgetab=tab, columnlist=tab.colnames[j], vshow=True,
                              plotstyle='image', clipedge=True, pad=pad, nx=nx, ny=ny,
                              cmap=cm, norm=norm, pdfname=outfile)
    return

