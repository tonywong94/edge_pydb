import os as _os
import json as _json
import shutil as _shutil
import requests as _requests


# Initial setup script to read the file location from config file
_ROOT = _os.path.abspath(_os.path.dirname(__file__))

_filepath = _os.path.join(_ROOT, '_config.json')
_runtime = False
_config = {}

try:
    _fp = open(_filepath, 'r+')
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


def _walkthrough(dir=_ROOT):
    retval = {}
    for _root, _dirs, _files in _os.walk(_os.path.abspath(dir)):
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
def updatefiles(dir=_ROOT):
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
    if readonly:
        _config = _json.load(_fp)
    else:
        _config.update(_json.load(_fp))
        _filepath = src
        _runtime = False
    _fp.close()


def listfiles(contain='', printing=False):
    '''
    List the current available files in the package data directory
    list all the file contain the specified substring
    list all the files otherwise

    Parameters:
        contain: the target substring to find in the file name to list, if not provided, then will print all files
    '''
    files = []
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
                raise FileNotFoundError("Cannot find the specified file: %s" % name)
            # if dir:
            #     dirpath = _os.path.abspath(_os.path.dirname(name))
            #     if dirpath not in retval:
            #         retval.append(path)
            else:
                retval.append(_config[name])
        return retval
    else:
        if names not in _config.keys():
            raise FileNotFoundError("Cannot find the specified file: %s" % names)
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


def add_from_dir(src, dest='', copy=True, overwrite=False):
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
        updatefiles(dest)
    else:
        updatefiles(src)

    if not _runtime:
        with open(_filepath, 'w') as _fp:
            _json.dump(_config, _fp)
    else:
        print("WARNING! The location of this file will be saved runtime only")
