# edge_pydb
Python-based database for CARMA EDGE.  This package requires Python 3.

The EDGE database has several components (see `index_csv.md` and `index_hdf.txt` for details):

* Zero-dimensional CSV tables (one value per galaxy), found in `dat_glob`.

* One-dimensional CSV tables (e.g. radial profiles or spectra), found in `dat_prof` and `dat_spec`.

* Downsampled 2D or 3D images, saved as HDF5 binary tables in the `img_` directories.

## Installation

Install the package directly from PyPI (https://pypi.org/project/edge-pydb/) using

    pip install --user edge_pydb

or, if you prefer to keep the latest source code handy, by cloning this Github repository and running

    pip install --user .

in the directory containing `setup.py`.

The `--user` flag ensures the package is not installed in your system-wide Python directories, which you probably don't have write access to.  The package tries to save and update a configuration file `_config.json`, which provides the paths to all the database tables, so installation in your user area is recommended.  Alternatively, you may save the configuration file in a different location using (e.g.):

    import edge_pydb.util as edgeutil
    edgeutil.save_config('config.json')

You will then need to load this file whenever you start the package:

    edgeutil.load_config('config.json')

The Github package only contains data for a single galaxy (NGC 4047), for demonstration and testing purposes.  Larger data files can be downloaded as ZIP archives from Box.  It is recommended that you unpack additional files into a single directory that is easily accessible on your file system, and not inside your Python libraries (`site-packages` area).  Here is the suggested way to incorporate these into your runtime environment:

Leave the git directory in which this README is located.  (This ensures that the package runs from your `site-packages` area and not the current directory.)  Open an iPython shell and type:

    edgeutil.listfiles(values=True)

This should show only the Github data installed in `site-packages`.  Now suppose the additional data files are in a folder called `pybase`.  Then

    edgeutil.add_from_dir('/path/to/pybase/', max_depth=0, copy=False)

will add the additional files to your environment.  This only needs to be done once after package installation, unless you add new files to `pybase`.  Here the `max_depth=0` parameter prevents files in subdirectories from being added.  Use the `listfiles` command above to verify that the expected tables are available.

## Basic Usage

    from edge_pydb import EdgeTable
    EdgeTable('list')

makes a listing of the available files.

    ctrpos = EdgeTable('edge_coflux_smo7.csv')

loads a CSV file.  `ctrpos` can now be treated like an astropy table, for example `ctrpos.info()` will summarize the contents and `ctrpos.pprint()` will print some of the data.

    ctrpos = EdgeTable('edge_coflux_smo7.csv', cols=['Name', 'coRactr_smo7', 'coDectr_smo7'])

loads the three specified columns only from the CSV file.

    leda  = EdgeTable('edge_leda.csv', cols=['Name', 'ledaD25', 'ledaPA', 'ledaIncl'])
    ctrpos.join(leda)

will merge a sub-table from `edge_leda.csv` into `ctrpos`.  We must select the `Name` column from both tables for the join to work.

    comom = EdgeTable('NGC4047.comom_smo7.hdf5', path='smo')

loads an HDF5 file.  The path must be given, otherwise a listing of available paths is provided.

A `demo_notebk` folder provides examples of accessing and plotting database values in a Jupyter notebook.

## HDF5 File Contents

Detailed listings of the HDF5 files are provided in [index_hdf.txt](https://github.com/tonywong94/edge_pydb/blob/master/index_hdf.txt) at the top level.  Note that each HDF5 file can bundle several tables or "paths," and only one path can be read into `EdgeTable` at a time.

- **[label].pipe3d.hdf5**: These are CALIFA data products from Pipe3D.  As described in [Sanchez et al. (2016a)](http://adsabs.harvard.edu/abs/2016RMxAA..52..171S), there are five collections of images, bundled as `ELINES`, `SFH`, `SSP`, `flux_elines`, and `indices`.  Each collection can be found in two separate Pipe3D runs, one performed on the native resolution CALIFA data, regridded to match the CARMA spatial grid (with 1" pixels), and another performed on the smoothed CALIFA data, matched to the 7" CARMA resolution (and also on the same spatial grid).  Thus you will find an `ELINES_rg` table as well as an `ELINES_sm` table with identical structure and data size; this yields 10 paths in total.  Being resolution matched, the `_sm` products are recommended when making comparisons with the CO data.

- **[label].comom_smo7.hdf5**, **[label].cocube_smo7.hdf5**: These are the CARMA CO moment maps and data cubes.  All are at a resolution of 7 arcsec (FWHM Gaussian beam).  Moment maps were generated using three different methods (`str`, `dil`, `smo`), with each method being a separate table (path) within the HDF5 file.  The straight (`str`) moment maps are generated without masking and have very poor signal-to-noise.  To reject noise, the dilated (`dil`) moment maps use a dilated mask that starts at a high significance contour (3.5&sigma; or greater in two consecutive channels) and expands to a surrounding 2&sigma; contour.  The smoothed (`smo`) moment maps use a mask that is obtained by smoothing the cube spatially (to 14") before constructing a dilated mask.  For most purposes the dilated masks produce the best results.

Available datasets will expand over time, but current values for [label] include 'edge' (which uses a square sampling grid spaced by 3", sufficient for Nyquist sampling the CARMA beam) and 'edge_allpix' (all pixel values saved, resulting in much longer tables).

## References

For more information about the CALIFA Data Release 3:

- [Sanchez et al. (2016b)](http://adsabs.harvard.edu/abs/2016A&A...594A..36S), "CALIFA, the Calar Alto Legacy Integral Field Area survey. IV. Third public data release."

For more information about EDGE:

- [Bolatto et al. (2017)](http://adsabs.harvard.edu/abs/2017ApJ...846..159B), "The EDGE-CALIFA Survey: Interferometric Observations of 126 Galaxies with CARMA."

If you use this package in a publication please also cite:

- Wong et al. (2021), "The EDGE-CALIFA Survey: An Extragalactic Database for Galaxy Evolution Studies," in prep.
