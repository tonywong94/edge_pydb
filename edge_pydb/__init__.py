# find configuration file
import os
import json


_ROOT = os.path.abspath(os.path.dirname(__file__))

filepath = os.path.join(_ROOT, 'config.json')
# print(_ROOT)
# print(filepath)
try:
    fp = open(filepath, 'r+')
    config = json.load(fp)

except:
    fp = open(filepath, 'w') 
    config = {}

if not config:
    # print(os.listdir(_ROOT))
    for root, dirs, files in os.walk(_ROOT):
        # print(files)
        for file in files:
            if file.endswith('.csv') or file.endswith('.hdf5'):
                # print(file)
                # print(os.path.join(_ROOT, file))
                config[file] = os.path.join(_ROOT, file)

json.dump(config, fp)
fp.close()

def list_file(file_type=None):
    suffix = ''
    if file_type:
        suffix = file_type

    for key in config:
        if key.endswith(suffix):
            print(key)
    return

def get_file(name):
    if name not in config:
        raise Exception("Does not find the specified file")

    return config[name]


from edge_pydb import edge_conv, fitsextract

# update the files