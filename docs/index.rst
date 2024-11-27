edge-pydb
=========

Python-based database for the `Extragalactic Database for Galaxy Evolution (EDGE) 
<https://www.astro.umd.edu/~bolatto/EDGE/>`_.

The EDGE database has several components:

-  Zero-dimensional CSV tables (one value per galaxy), found in ``dat_glob``
   of the ``site-packages`` directory.

-  One-dimensional CSV tables (e.g. radial profiles or spectra), found 
   in ``dat_prof`` and ``dat_spec``.

-  Downsampled 2D or 3D images, saved as HDF5 binary tables in the
   ``img_`` directories, or as large HDF5 files in a user-specified area.

For the latest version matching this documentation, please install from the
`Github repository <https://github.com/tonywong94/edge_pydb>`_.

If you use this package in a publication please cite:

-  `Wong et
   al. (2024) <http://adsabs.harvard.edu/abs/2024ApJS..271...35W>`__,
   “The EDGE-CALIFA Survey: An Extragalactic Database for Galaxy
   Evolution Studies.”

Basic Usage
^^^^^^^^^^^

::

    from edge_pydb import EdgeTable
    EdgeTable('list')

makes a listing of the available files.

::

    ctrpos = EdgeTable('edge_coflux_smo7.csv')

loads a CSV file.  ``ctrpos`` can now be treated like an `astropy table <https://docs.astropy.org/en/stable/table/>`_, 
for example ``ctrpos.info()`` will summarize the contents and ``ctrpos.pprint()`` 
will print some of the data.

::

    ctrpos = EdgeTable('edge_coflux_smo7.csv', cols=['Name', 'coRactr_smo7', 'coDectr_smo7'])

loads the three specified columns only from the CSV file.

::

    leda = EdgeTable('edge_leda.csv', cols=['Name', 'ledaD25', 'ledaPA', 'ledaIncl'])
    ctrpos.join(leda)

will merge a sub-table from ``edge_leda.csv`` into ``ctrpos``.  We must select the 
``Name`` column from both tables for the join to work.  For pixel tables the ``ix`` and
``iy`` columns must also be selected.

::

    comom = EdgeTable('NGC4047.2d_smo7.hdf5', path='comom_smo')

loads an HDF5 file.  The path must be given, otherwise a listing of available paths is provided.

A `demo_notebk <https://github.com/tonywong94/edge_pydb/tree/main/demo_notebk>`_ folder 
and various subfolders provide examples of accessing and plotting 
database values in a Jupyter notebook.

More Information
^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   install.rst
   Tutorial
   global_props
   SFLaw-edge
   Alambda_demo
   api.rst

