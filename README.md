# edge_pydb
Python-based database for CARMA EDGE

The EDGE database has several components:
* Zero-dimensional tables (one value per galaxy), found in `dat_glob`.
* One-dimensional tables (e.g. radial profiles or spectra), found in `dat_prof`.
* Downsampled 2D or 3D images, saved as HDF5 binary tables in the `img_` directories.

Sub-directories provide a space for building the database from the original text or FITS files.  Sample HDF5 files are provided for one galaxy for testing purposes; the entire set is too large to host on a Github repo.

A `demo_notebk` folder provides examples of accessing database values from a Jupyter notebook.
