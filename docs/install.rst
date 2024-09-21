Installation
------------

Required packages are:

-  `matplotlib <https://matplotlib.org/>`__
-  `numpy <https://numpy.org/>`__
-  `scipy <https://scipy.org/>`__
-  `astropy <https://www.astropy.org/>`__
-  `pandas <https://pandas.pydata.org/>`__
-  `radio_beam <https://radio-beam.readthedocs.io/>`__
-  `requests <https://requests.readthedocs.io/>`__
-  `h5py <https://www.h5py.org/>`__
-  `uncertainties <https://uncertainties.readthedocs.io/>`__

Also needed for building the database (but not for general use) are:

-  `pyFIT3D <https://ifs.astroscu.unam.mx/pyPipe3D/>`__
-  `reproject <https://reproject.readthedocs.io/>`__
-  `CO_conversion_factor <https://github.com/astrojysun/COConversionFactor/>`__

Install the package directly from PyPI
(https://pypi.org/project/edge-pydb/) using

::

   pip install --user edge_pydb

or, if you prefer to keep the latest source code handy, by cloning the
`Github repository <https://github.com/tonywong94/edge_pydb>`_ and running

::

   pip install --user .

in the directory containing ``setup.py``.

The ``--user`` flag ensures the package is not installed in your
system-wide Python directories, which you may not have write access to.
The package tries to save and update a configuration file
``_config.json``, which provides the paths to all the database tables,
so installation in your user area is recommended. Alternatively, you may
save the configuration file in a different location using (e.g.):

::

   import edge_pydb.util as edgeutil
   edgeutil.save_config('/path/to/config.json')

You will then need to load this file whenever you start the package:

::

   edgeutil.load_config('/path/to/config.json')

To uninstall you may use

::

   pip uninstall edge_pydb

but note that the ``_config.json`` will not be removed, so to fully
uninstall you will need to delete the ``edge_pydb`` directory manually
using (e.g.):

::

   rm -r ~/Library/Python/3.9/lib/python/site-packages/edge_pydb/

Merging in large data sets
--------------------------

The Github package only contains data for a single galaxy (NGC 4047),
for demonstration and testing purposes. Larger data files can be
downloaded from Zenodo:

-  `HDF5 Files for CARMA EDGE <https://zenodo.org/records/10256732>`__

It is recommended that you unpack additional files into a single
directory that is easily accessible on your file system, and not
embedded within your Python libraries (``site-packages`` area). Here is
the suggested way to incorporate these into your runtime environment.

Ensure you are *not* in the directory in which ``setup.py`` is located, since 
you want to run the package from your ``site-packages`` area and not the
current directory. Open an iPython shell and type:

::

   import edge_pydb.util as edgeutil
   edgeutil.listfiles(values=True)

This should show only the Github data installed in ``site-packages``.
Now suppose the additional data files are in a folder called ``pybase``.
Then type:

::

   edgeutil.add_from_dir('/path/to/pybase/', max_depth=0, copy=False)

to add the additional files to your environment. This only needs to be
done once after package installation, unless you add new files to
``pybase``. Here the ``max_depth=0`` parameter prevents files in
subdirectories from being added. Use the ``listfiles`` command above to
verify that the expected tables are available.
