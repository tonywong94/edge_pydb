# find _configuration _file
import os as _os
import json as _json
import shutil as _shutil



_ROOT = _os.path.abspath(_os.path.dirname(__file__))

_filepath = _os.path.join(_ROOT, '_config.json')
# print(_ROOT)
# print(_filepath)
try:
    _fp = open(_filepath, 'r+')
    _config = _json.load(_fp)

except:
    _fp = open(_filepath, 'w') 
    _config = {}


def _walkthrough(dir=''):
    retval = {}
    if not dir:
        dir = _ROOT
    for _root, _dirs, _files in _os.walk(dir):
        # print(_files)
        for _file in _files:
            if _file.endswith('.csv') or _file.endswith('.hdf5') and _file not in retval:
                # print(_file)
                # print(os.path.join(_ROOT, _file))
                # print(_dirs)
                retval[_file] = _os.path.join(_root, _file)
    return retval


if not _config:
    # print(os.listdir(_ROOT))
    _config = _walkthrough()

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
def Updatefiles():
    tmp = _walkthrough()
    for k, v in tmp:
        if k not in _config:
            print("Added file %s" % k)
        elif v != _config[k]:
            print("Updated file %s" % k)
    _config = tmp
    with open(_filepath, 'w') as _fp:
        _json.dump(_config, _fp)


def Addfile(src, dest='', copy=True):
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
    with open(_filepath, 'w') as _fp:
        _json.dump(_config, _fp)


def AddDir(src, dest='', copy=True, overwrite=False):
    dirname = _os.path.basename(src)
    # dirname = _os.path.abspath(_os.path.dirname(src))
    if not dest:
        dest = _ROOT + '/' + dirname
        
    if copy:
        try:
            _shutil.copytree(src, dest)
        except:
            if overwrite:
                _shutil.rmtree(dest)
                _shutil.copytree(src, dest)
            else:
                _shutil.copy(src, dest)
        _config.update(_walkthrough(dest))
    else:
        _config.update(_walkthrough(src))

    with open(_filepath, 'w') as _fp:
        _json.dump(_config, _fp)


from edge_pydb import conversion, fitsextract
