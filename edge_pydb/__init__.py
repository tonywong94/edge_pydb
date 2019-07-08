# find _configuration _file
import os as _os
import json as _json


_ROOT = _os.path.abspath(_os.path.dirname(__file__))

_filepath = _os.path.join(_ROOT, 'config.json')
# print(_ROOT)
# print(_filepath)
try:
    _fp = open(_filepath, 'r+')
    _config = _json.load(_fp)

except:
    _fp = open(_filepath, 'w') 
    _config = {}

if not _config:
    # print(os.listdir(_ROOT))
    for _root, _dirs, _files in _os.walk(_ROOT):
        # print(_files)
        for _file in _files:
            if _file.endswith('.csv') or _file.endswith('.hdf5'):
                # print(_file)
                # print(os.path.join(_ROOT, _file))
                # print(_dirs)
                _config[_file] = _os.path.join(_root, _file)

_json.dump(_config, _fp)
_fp.close()

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

def getfiles(names):
    if isinstance(names, list):
        retval = []
        for name in names:
            if name not in _config:
                raise Exception("Cannot find the specified file: %s" % name)
            retval.append(_config[name])
        return retval
    else:
        if names not in _config:
            raise Exception("Cannot find the specified file: %s" % names)
        return _config[names]


from edge_pydb import edge_conv, fitsextract

# update the _files