import os as _os
import json as _json
import shutil as _shutil



_ROOT = _os.path.abspath(_os.path.dirname(__file__))

_filepath = _os.path.join(_ROOT, '_config.json')
_runtime = False

try:
    _fp = open(_filepath, 'r+')
    if _os.stat(_filepath).st_size == 0:
        _config = {}
    else:
        _config = _json.load(_fp)

except FileNotFoundError:
    _fp = open(_filepath, 'w') 
    _config = {}

except OSError as _err:
    if _err.errno == 30:
        print("WARNING! Read-only file system, cannot record the package data file location.\n" + \
            "Please consider change the mode of this file for reliability.")
        print("If you need to change the files in the package data, please consider run as root, the manipulation of files requires the sudo priority.")
    _runtime = True
    _config = {}

def _walkthrough(dir=_ROOT):
    retval = {}
    for _root, _dirs, _files in _os.walk(dir):
        for _file in _files:
            if _file.endswith('.csv') or _file.endswith('.hdf5'):
                if _file in retval:
                    print("{} location will be changed\nCurrent location: {}\nNew location: {}".format(_file, retval[_file], _os.path.join(_root, _file)))
                retval[_file] = _os.path.join(_root, _file)
    return retval


if not _config:
    # print(os.listdir(_ROOT))
    _config = _walkthrough()

if not _runtime:
    _json.dump(_config, _fp)
    _fp.close()


'''
List the current available files in the package data directory
list all the specified file type such as hdf5 or csv
list all the files otherwise
'''

def listfiles(file_type=None):
    suffix = ''
    if file_type:
        suffix = file_type

    _files = []
    for key in _config:
        if key.endswith(suffix):
            print(key)
            _files.append(key)
    return _files
    

'''
Get all the files by its file name, can either be a single file or a list of files
'''

def getfiles(names, dir=False):
    if isinstance(names, list):
        retval = []
        for name in names:
            if name not in _config:
                raise Exception("Cannot find the specified file: %s" % name)
            if dir:
                dirpath = _os.path.abspath(_os.path.dirname(name))
                if dirpath not in retval:
                    retval.append(path)
            else:
                retval.append(name)
        return retval
    else:
        if names not in _config:
            raise Exception("Cannot find the specified file: %s" % names)
        if dir:
            return _os.path.abspath(_os.path.dirname(names))
        else:
            return _config[names]


# update the files
def Updatefiles(dir=_ROOT):
    tmp = _walkthrough(dir)
    for k, v in tmp:
        if k not in _config:
            print("Add file %s" % k)
        elif v != _config[k]:
            print("Update file %s" % k)
    _config.update(tmp)
    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)


def Addfile(src, dest='', copy=True):
    if _runtime:
        print("WARNING! No sudo permission, take care, will break")
    name = _os.path.basename(src)
    src = _os.path.abspath(src)
    if copy:
        if dest:
            if name in _config:
                os.remove(_config[name])
            _shutil.copyfile(src, dest)
            _config[name] = dest
        else:
            if not _os.path.exists(_ROOT + '/data'):
                _os.mkdir(_ROOT + '/data')

            _shutil.copyfile(src, _ROOT + '/data/' + name)
            _config[name] = _ROOT + '/data/' + name
    else:
        _config[name] = src

    print("Updated file %s" % name)    
    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)
    else:
        print("WARNING! The location of this file will be recorded runtime only")


def AddDir(src, dest='', copy=True, overwrite=False):
    if _runtime:
        print("WARNING! No sudo permission, take care, will break")
    dirname = _os.path.basename(src)
    # dirname = _os.path.abspath(_os.path.dirname(src))
    if not dest:
        dest = _ROOT + '/' + dirname
        
    if copy:
        try:
            _shutil.copytree(src, dest)
        except FileExistsError:
            if overwrite:
                _shutil.rmtree(dest)
                _shutil.copytree(src, dest)
            else:
                _shutil.copy(src, dest)
        Updatefiles(dest)
    else:
        Updatefiles(src)

    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)
    else:
        print("WARNING! The location of this file will be recorded runtime only")


from edge_pydb import conversion, fitsextract, xy2hist
