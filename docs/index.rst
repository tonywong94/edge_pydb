edge-pydb
=========

Python-based database for the Extragalactic Database for Galaxy Evolution (EDGE).

The EDGE database has several components:

-  Zero-dimensional CSV tables (one value per galaxy), found in ``dat_glob``.

-  One-dimensional CSV tables (e.g. radial profiles or spectra), found 
   in ``dat_prof`` and ``dat_spec``.

-  Downsampled 2D or 3D images, saved as HDF5 binary tables in the
   ``img_`` directories, or as large HDF5 files in a user-specified area.

Getting started
^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   install.rst
   api.rst

