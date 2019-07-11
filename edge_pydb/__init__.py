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
            "For better longterm performance, consider provide a config file by using the extConfig() function.")
        print("If you need to change the files in the package data, please consider run as root, the manipulation of files requires the sudo priority.")
    _runtime = True
    _config = {}

def _walkthrough(dir=_ROOT):
    retval = {}
    for _root, _dirs, _files in _os.walk(_os.path.abspath(dir)):
        for _file in _files:
            if _file.endswith('.csv') or _file.endswith('.hdf5'):
                if _file in retval:
                    print("{} redundant file detected\n--Current location: {}\n++New location: {}".format(_file, retval[_file], _os.path.join(_root, _file)))
                retval[_file] = _os.path.join(_root, _file)
    return retval


if not _config:
    # print(os.listdir(_ROOT))
    _config = _walkthrough()
    if not _runtime:
        _json.dump(_config, _fp)

if not _runtime:
    _fp.close()


'''
This function will read the config from a file or write back to a file.
If the persistent is true, then will directly use the provided file rather than the _config in the package directory.
'''
def extConfig(src, mode, persistent=False):
    valid = {'w', 'r', 'u'}
    if mode not in valid:
        print("The valid mode is\n" + \
            "'w': write\n'r': read\n'u':update, update the current config with the new file")
        return
    if mode == 'w':
        _fp = open(src, 'w')
        _json.dump(_config, _fp)
    else:
        _fp = open(src, 'r')
        if mode == 'r':
            _config = _json.load(_fp)
        elif mode == 'u':
            _config.update(_json.load(_fp))
    _fp.close()

    if persistent:
        _filepath = src
        _runtime = False



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
    for key in _config.keys():
        if key.endswith(suffix):
            print(key)
            _files.append(key)
    return _files
    

'''
Get all the files by its file name, can either be a single file or a list of files
'''

def getfiles(names):
    if isinstance(names, list):
        retval = []
        for name in names:
            if name not in _config.keys():
                raise Exception("Cannot find the specified file: %s" % name)
            # if dir:
            #     dirpath = _os.path.abspath(_os.path.dirname(name))
            #     if dirpath not in retval:
            #         retval.append(path)
            else:
                retval.append(_config[name])
        return retval
    else:
        if names not in _config.keys():
            raise Exception("Cannot find the specified file: %s" % names)
        # if dir:
        #     return _os.path.abspath(_os.path.dirname(names))
        else:
            return _config[names]


# update the files
def Updatefiles(dir=_ROOT):
    tmp = _walkthrough(dir)
    for k, v in tmp.items():
        if k not in _config.keys():
            print("Add file %s" % k)
        elif v != _config[k]:
            print("Update file %s" % k)
    _config.update(tmp)
    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)


def Addfile(src, dest='', copy=True, overwrite=False):
    if _runtime:
        print("WARNING! No sudo permission, take care, will break")
    name = _os.path.basename(src)
    src = _os.path.abspath(src)
    if copy:
        if dest:
            if name in _config.keys():
                if overwrite:
                    os.remove(_config[name])
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
                dest = _ROOT + '/data/' + dirname
                _shutil.copytree(src, dest)
        Updatefiles(dest)
    else:
        Updatefiles(src)

    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)
    else:
        print("WARNING! The location of this file will be recorded runtime only")


from edge_pydb import conversion, fitsextract, xy2hist
