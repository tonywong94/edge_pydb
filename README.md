# edge_pydb
Python-based database for CARMA EDGE

The EDGE database has several components:
* Zero-dimensional tables (one value per galaxy), found in `dat_glob`.
* One-dimensional tables (e.g. radial profiles or spectra), found in `dat_prof` and `dat_spec`.
* Downsampled 2D or 3D images, saved as HDF5 binary tables in the `img_` directories.

Install the package directly from PyPI (https://pypi.org/project/edge-pydb/) using

`pip install --user edge_pydb`

or, if you prefer to keep the latest source code handy, by cloning this Github repository and running

`pip install --user .`

in the directory containing `setup.py`.

The `--user` flag ensures the package is not installed in your system-wide Python directories, which you probably don't have write access to.  The package tries to save and update a configuration file `_config.json`, so installation in your user area is recommended.

If you have access to the protected archive (not on Github), you should also obtain the `team_files.py` script from Tony which will download and install the additional tables into the appropriate locations.  You can run this script in a shell or iPython environment (but should avoid doing so in the same directory as `setup.py`).  Currently the Github repository only includes HDF5 files for one galaxy, for testing purposes.

A `demo_notebk` folder provides examples of accessing database values from a Jupyter notebook.
